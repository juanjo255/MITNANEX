"""Reference-guided pipeline: mitochondrial reads, consensus and variant calls.

Port of ``mitnanex_reference.sh`` and ``scripts/*.sh``:

1. quality / length filter                       (chopper)
2. map to the reference mitogenome               (minimap2, samtools)
3. Flye --meta assembly of the mapped reads to discard NUMT/low-quality reads,
   keep the reads that map to the mitogenome contig (consensus)
4. re-map those reads to the reference           (minimap2 with read group for GATK)
5. variant calling in mitochondrial mode         (GATK Mutect2 + FilterMutectCalls)
6. figures: depth along the genome, allele fractions, spectrum (pandas/matplotlib)
7. optional: annotate variants with rCRS regions (bcftools), lift genes over
   to the consensus (liftoff), haplogroups (haplogrep3), MITOMAP query.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

import pandas as pd

from .. import vcf as vcf_utils
from ..external import check_tools, run, run_pipeline
from ..io import iter_fastq, read_fasta, read_flye_assembly_info, write_fasta
from .base import Pipeline, default_prefix

log = logging.getLogger("mitnanex")

# Bundled human mitochondrial reference (rCRS, NC_012920.1) and its annotation.
REFSEQ_MT = resources.files("mitnanex") / "data" / "refseqMT"

# BED annotations of rCRS regions: file -> (INFO tag, VCF header line).
RCRS_ANNOTATIONS = {
    "HV.bed": ("Hypervariable", '##INFO=<ID=Hypervariable,Number=1,Type=String,Description="Hypervariable region">'),
    "HP.bed": ("Homopolymer", '##INFO=<ID=Homopolymer,Number=0,Type=Flag,Description="Homopolymer">'),
    "HS.bed": ("Hotspot", '##INFO=<ID=Hotspot,Number=0,Type=Flag,Description="Mutational hotspot">'),
    "CDS.bed": ("CDS", '##INFO=<ID=CDS,Number=1,Type=String,Description="Protein-coding gene">'),
    "RNR.bed": ("RNR", '##INFO=<ID=RNR,Number=1,Type=String,Description="rRNA gene">'),
    "TRN.bed": ("TRN", '##INFO=<ID=TRN,Number=1,Type=String,Description="tRNA gene">'),
    "DLOOP.bed": ("DLOOP", '##INFO=<ID=DLOOP,Number=0,Type=Flag,Description="D-loop">'),
}


@dataclass
class ReferenceConfig:
    """Options of the reference pipeline (CLI flags in parentheses).

    Attributes:
        reference: Reference mitogenome FASTA, single contig (``-r``).
        reads: Input FASTQ (``-i``).
        gff: Reference GFF3 for liftoff (``-g``); default bundled rCRS annotation.
        outdir / prefix / threads: Output folder, file prefix and threads.
        minimap2_preset: minimap2 ``-x`` preset (``--mm2-preset``).
        min_mapq: Minimum mapping quality (``--mapq``).
        filter_reads: Run chopper (``--no-filter`` disables it).
        min_length / max_length / min_mean_quality: chopper options.
        keep_percent: Filtlong ``--keep_percent`` before the NUMT-removal assembly.
        min_pruning / kmer_sizes / max_assembly_region_size: GATK Mutect2 options.
        flye_preset / other_flye_opts: Flye options.
        annotate_rcrs: Annotate variants with rCRS regions (only for rCRS references).
        liftoff: Lift genes from the reference to the consensus.
        haplogrep: Path to the haplogrep3 executable (classification is skipped if ``None``).
        haplogrep_trees / top_hits: haplogrep3 options.
        mitomap: Send PASS variants to the MITOMAP web service (opt-in).
    """

    reference: Path
    reads: Path
    gff: Path | None = None
    outdir: Path = Path("mitnanex_results")
    prefix: str | None = None
    threads: int = 8
    minimap2_preset: str = "map-ont"
    min_mapq: int = 30
    filter_reads: bool = True
    min_length: int = 500
    max_length: int | None = None
    min_mean_quality: float = 10
    keep_percent: float = 80
    min_pruning: int = 3
    kmer_sizes: list[int] = field(default_factory=lambda: [15, 25])
    max_assembly_region_size: int | None = None
    flye_preset: str = "--nano-hq"
    other_flye_opts: list[str] = field(default_factory=list)
    annotate_rcrs: bool = False
    liftoff: bool = False
    haplogrep: Path | None = None
    haplogrep_trees: list[str] = field(default_factory=lambda: ["phylotree-rcrs@17.2"])
    top_hits: int = 3
    mitomap: bool = False


class ReferencePipeline(Pipeline):
    """Reference-guided extraction of mitochondrial reads and variant calling."""

    name = "reference"
    required_tools = ("minimap2", "samtools", "filtlong", "flye", "gatk")

    def __init__(self, config: ReferenceConfig, timestamped: bool = False, resume: bool = False) -> None:
        super().__init__(config.outdir, config.prefix or default_prefix(config.reads), timestamped, resume)
        self.cfg = config
        self.cfg.reads = Path(config.reads).resolve()
        if config.filter_reads:
            self.required_tools = (*self.required_tools, "chopper")

        # The reference must be a single contig (the mitogenome).
        reference = read_fasta(config.reference)
        if len(reference) != 1:
            raise ValueError(f"The reference must contain exactly one contig (found {len(reference)}).")
        self.contig_id = reference.loc[0, "id"]
        self.contig_len = int(reference.loc[0, "length"])
        # Reads longer than the mitogenome are likely chimeric/NUMT-containing (old default).
        self.cfg.max_length = config.max_length or self.contig_len

        # File names in one place.
        self.ref_dir = self.outdir / "reference"
        self.ref = self.ref_dir / f"{self.contig_id}.fasta"
        self.filtered_reads = self.path(f".filtQ{config.min_mean_quality:g}.fastq")
        self.flye_dir = self.outdir / "flye_for_numts"
        self.consensus = self.outdir / "MT_genome.fasta"
        self.mt_reads = self.path(".mt_reads.fastq")
        self.bam = self.path(f".{self.contig_id}.sorted.bam")
        self.gatk_dir = self.outdir / "VariantCall" / "gatk_mutect2"
        self.vcf_raw = self.gatk_dir / f"{self.prefix}.{self.contig_id}.gatk.vcf"
        self.vcf = self.gatk_dir / f"{self.prefix}.{self.contig_id}.gatk.filt.vcf"

    def steps(self):
        cfg = self.cfg
        steps = [("Prepare reference (copy, faidx, GATK dictionary)", self.prepare_reference, [self.ref.with_suffix(".dict")])]
        if cfg.filter_reads:
            steps.append((f"Filter reads (Q>={cfg.min_mean_quality:g}, {cfg.min_length}-{cfg.max_length} bp; chopper)",
                          self.filter_reads, [self.filtered_reads]))
        steps += [
            ("Map to reference, assemble with Flye to remove NUMTs, re-map", self.map_reads, [self.bam]),
            ("Variant calling (GATK Mutect2, mitochondria mode)", self.variant_calling, [self.vcf]),
            ("Variant figures", self.plot_variants, []),
        ]
        if cfg.annotate_rcrs:
            steps.append(("Annotate variants with rCRS regions (bcftools)", self.annotate_vcf, []))
        if cfg.liftoff:
            steps.append(("Lift annotation to the consensus (liftoff)", self.liftoff, []))
        if cfg.haplogrep:
            steps.append(("Haplogroup classification (haplogrep3)", self.haplogroups, []))
        if cfg.mitomap:
            steps.append(("Annotate variants with MITOMAP (web service)", self.mitomap, []))
        return steps

    @property
    def reads(self) -> Path:
        """Reads entering the mapping step (filtered or raw)."""
        return self.filtered_reads if self.cfg.filter_reads else self.cfg.reads

    # -- steps -------------------------------------------------------------------
    def prepare_reference(self) -> None:
        """Copy the reference into the output folder and index it there.

        The original scripts expected ``.fai``/``.dict`` files next to the user's
        reference (the indexing commands were commented out), so GATK failed
        for any reference other than the bundled one.
        """
        self.ref_dir.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(self.cfg.reference, self.ref)
        run(["samtools", "faidx", self.ref], log_file=self.tool_log)
        self.ref.with_suffix(".dict").unlink(missing_ok=True)
        run(["gatk", "CreateSequenceDictionary", "-R", self.ref, "-O", self.ref.with_suffix(".dict")], log_file=self.tool_log)

    def filter_reads(self) -> None:
        cfg = self.cfg
        run(["chopper", "--threads", cfg.threads, "-q", cfg.min_mean_quality, "--minlength", cfg.min_length,
             "--maxlength", cfg.max_length, "--input", cfg.reads], stdout=self.filtered_reads, log_file=self.tool_log)

    def _map(self, reference: Path, reads: Path, bam: Path, extra: list | None = None, min_mapq: int | None = None) -> None:
        """minimap2 | samtools view (mapped primary reads) | samtools sort, then index."""
        cfg = self.cfg
        mapq = cfg.min_mapq if min_mapq is None else min_mapq
        run_pipeline(
            [["minimap2", "-t", cfg.threads, "--secondary=no", "-ax", cfg.minimap2_preset, *(extra or []), reference, reads],
             ["samtools", "view", "-@", cfg.threads, "-b", "-q", mapq, "-F", "0x804", "-"],
             ["samtools", "sort", "-@", cfg.threads, "-o", bam, "-"]],
            log_file=self.tool_log,
        )
        run(["samtools", "index", bam], log_file=self.tool_log)

    def map_reads(self) -> None:
        """Mapping + Flye-based NUMT removal (port of ``scripts/map_reads.sh``)."""
        cfg = self.cfg
        # 1. reads mapped to the reference mitogenome
        first_bam = self.path(".first_map.bam")
        self._map(self.ref, self.reads, first_bam, extra=["-g", "1k"])
        first_reads = self.path(".first_map.fastq")
        run(["samtools", "fastq", "-@", cfg.threads, first_bam], stdout=first_reads, log_file=self.tool_log)
        filtlong_reads = self.path(".first_map.filtlong.fastq")
        run(["filtlong", "--keep_percent", cfg.keep_percent, "--mean_q_weight", "10", first_reads],
            stdout=filtlong_reads, log_file=self.tool_log)

        # 2. assemble them: NUMT and low-quality reads end up in other contigs
        run(["flye", "-t", cfg.threads, "--meta", cfg.flye_preset, filtlong_reads, *cfg.other_flye_opts,
             "-o", self.flye_dir], log_file=self.tool_log)
        info = read_flye_assembly_info(self.flye_dir / "assembly_info.txt")
        # Mitogenome contig = the longest one (as in the original script).
        contig = info.sort_values("length", ascending=False).iloc[0]
        write_fasta(read_fasta(self.flye_dir / "assembly.fasta", ids=[contig["seq_name"]]), self.consensus, line_width=80)
        log.info("Consensus mitogenome: %s (%s bp, %sx, circular=%s)", contig["seq_name"], contig["length"],
                 contig["coverage"], contig["circular"])

        # 3. keep only the reads that map to the mitogenome contig
        flye_bam = self.flye_dir / f"aln_{self.prefix}.sorted.bam"
        self._map(self.flye_dir / "assembly.fasta", filtlong_reads, flye_bam, extra=["-k", "25"])
        contig_bam = self.flye_dir / f"aln_{self.prefix}_{contig['seq_name']}.bam"
        run(["samtools", "view", "-@", cfg.threads, "-b", "-o", contig_bam, flye_bam, contig["seq_name"]], log_file=self.tool_log)
        run(["samtools", "fastq", "-@", cfg.threads, contig_bam], stdout=self.mt_reads, log_file=self.tool_log)

        # 4. final alignment to the reference, with a read group (required by GATK)
        self._map(self.ref, self.mt_reads, self.bam, extra=["-R", r"@RG\tID:sample\tSM:" + self.prefix], min_mapq=0)
        for tmp in (first_bam, first_bam.with_suffix(".bam.bai"), first_reads, filtlong_reads):
            tmp.unlink(missing_ok=True)

        self.summary.update(
            consensus=str(self.consensus), consensus_length=int(contig["length"]),
            consensus_circular=contig["circular"] == "Y", mt_reads=sum(1 for _ in iter_fastq(self.mt_reads)),
        )
        log.info("Mitochondrial reads: %s (%d reads)", self.mt_reads, self.summary["mt_reads"])

    def variant_calling(self) -> None:
        """GATK Mutect2 in mitochondria mode followed by FilterMutectCalls."""
        cfg = self.cfg
        self.gatk_dir.mkdir(parents=True, exist_ok=True)
        mutect = ["gatk", "Mutect2", "-R", self.ref, "-L", self.contig_id, "--mitochondria-mode",
                  "--annotation", "DepthPerAlleleBySample", "--min-pruning", cfg.min_pruning,
                  "-I", self.bam, "-O", self.vcf_raw]
        for k in cfg.kmer_sizes:
            mutect += ["--kmer-size", k]
        if cfg.max_assembly_region_size:
            mutect += ["--max-assembly-region-size", cfg.max_assembly_region_size]
        run(mutect, log_file=self.tool_log)
        run(["gatk", "FilterMutectCalls", "--mitochondria-mode", "-R", self.ref, "-V", self.vcf_raw, "-O", self.vcf],
            log_file=self.tool_log)
        variants = vcf_utils.read_vcf(self.vcf)
        variants.drop(columns=["INFO"]).to_csv(self.vcf.with_suffix(".tsv"), sep="\t", index=False)
        self.summary.update(variants_total=len(variants), variants_pass=int(variants["pass"].sum()), vcf=str(self.vcf))
        log.info("Variants: %d called, %d PASS (%s)", len(variants), variants["pass"].sum(), self.vcf)

    def plot_variants(self) -> None:
        variants = vcf_utils.read_vcf(self.vcf)
        depth = vcf_utils.read_depth(self.bam)
        plot_dir = self.outdir / "plots"
        vcf_utils.plot_variants_along_genome(variants, depth, plot_dir / f"{self.prefix}_variants_along_genome.png")
        vcf_utils.plot_variant_summary(variants, plot_dir / f"{self.prefix}_variant_summary.png")
        self.summary["median_depth"] = float(depth["depth"].median())

    def annotate_vcf(self) -> None:
        """Annotate variants with rCRS regions (port of ``scripts/annotate_vcf.sh``).

        Fixes two bugs of the original: the bundled BED files name the contig
        ``chrM`` while the reference/VCF use ``NC_012920.1`` (nothing matched),
        and bcftools needs the annotation files bgzipped and tabix-indexed.
        """
        check_tools(["bcftools", "bgzip", "tabix"])
        if self.contig_len != 16569:
            log.warning("annotate_rcrs: the reference is not rCRS-sized (%d bp); annotations may be wrong.", self.contig_len)
        annot_dir = self.outdir / "annotation" / "rcrs_regions"
        annot_dir.mkdir(parents=True, exist_ok=True)
        current = self.vcf
        for bed_name, (tag, header) in RCRS_ANNOTATIONS.items():
            bed = pd.read_csv(REFSEQ_MT / bed_name, sep="\t", header=None).iloc[:, :4]
            bed[0] = self.contig_id  # rename chrM -> the VCF contig name
            bed_path = annot_dir / f"{bed_name}"
            bed.sort_values([1, 2]).to_csv(bed_path, sep="\t", header=False, index=False)
            run(["bgzip", "-f", bed_path], log_file=self.tool_log)
            run(["tabix", "-f", "-p", "bed", f"{bed_path}.gz"], log_file=self.tool_log)
            header_file = annot_dir / f"{tag}.hdr"
            header_file.write_text(header + "\n")
            annotated = annot_dir / f"{self.prefix}.{tag}.vcf"
            run(["bcftools", "annotate", "-a", f"{bed_path}.gz", "-c", f"CHROM,FROM,TO,{tag}", "-h", header_file,
                 "-o", annotated, current], log_file=self.tool_log)
            current = annotated
        final = self.vcf.with_name(self.vcf.name.replace(".vcf", ".annotated.vcf"))
        shutil.move(current, final)
        self.summary["annotated_vcf"] = str(final)
        log.info("Annotated VCF: %s", final)

    def liftoff(self) -> None:
        """Lift genes from the reference annotation onto the consensus (``scripts/annotate_mito.sh``)."""
        check_tools(["liftoff"])
        gff = self.cfg.gff or REFSEQ_MT / "chrMT_NC_012920.gff3"
        out = self.outdir / "annotation" / f"{self.prefix}.liftoff.gff"
        out.parent.mkdir(parents=True, exist_ok=True)
        run(["liftoff", "-g", gff, "-mm2_options=-a --end-bonus 0 --eqx -N 5 -p 0.8 -k10 -w3",
             "-o", out, self.consensus, self.ref], log_file=self.tool_log, cwd=out.parent)
        log.info("Consensus annotation: %s", out)

    def haplogroups(self) -> None:
        """Classify the haplogroup with haplogrep3 (``scripts/haplogroup_class.sh``)."""
        out_dir = self.outdir / "haplogroup"
        out_dir.mkdir(parents=True, exist_ok=True)
        for tree in self.cfg.haplogrep_trees:
            run([self.cfg.haplogrep, "install-tree", tree], log_file=self.tool_log)
            run([self.cfg.haplogrep, "classify", f"--tree={tree}", "--in", self.vcf, "--hits", self.cfg.top_hits,
                 "--extend-report", "--out", out_dir / f"haplogrep3.{tree}.txt"], log_file=self.tool_log)

    def mitomap(self) -> None:
        """Send PASS variants to MITOMAP's MITOMASTER web service (``scripts/mitomap.py``).

        Opt-in only: this uploads the sample's variant list to mitomap.org.
        """
        import requests  # optional dependency, only needed here

        variants = vcf_utils.read_vcf(self.vcf)
        snv_list = self.outdir / "variants.mitomaster_input.tsv"
        variants.loc[variants["pass"], ["CHROM", "POS", "REF", "ALT"]].to_csv(snv_list, sep="\t", header=False, index=False)
        with open(snv_list, "rb") as handle:
            response = requests.post(
                "https://mitomap.org/mitomaster/websrvc.cgi",
                files={"file": handle, "fileType": ("", "snvlist"), "output": ("", "detail")},
                timeout=300,
            )
        response.raise_for_status()
        out = self.outdir / "variant_annot.mitomap.txt"
        out.write_text(response.text)
        log.info("MITOMAP annotation: %s", out)
