"""Command-line interface: ``mitnanex {denovo,reference,extract} ...``.

The short flags of the original shell scripts are kept for backwards
compatibility (e.g. ``mitnanex denovo -i reads.fastq -t 8 -p 0.4 -s 0.6``).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .external import ToolError
from .extract import ExtractionConfig, MitoReadExtractor
from .pipelines import DenovoConfig, DenovoPipeline, ReferenceConfig, ReferencePipeline
from .pipelines.base import setup_logging

log = logging.getLogger("mitnanex")


def _proportion(value: str) -> float:
    """``-p`` accepts a fraction (0-1] or a number of reads (> 1)."""
    number = float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be > 0")
    return int(number) if number > 1 else number


def _int_list(value: str) -> list[int]:
    return [int(v) for v in value.split(",") if v]


def _add_extraction_args(parser: argparse.ArgumentParser) -> None:
    """Options of the read clustering / selection, shared by `denovo` and `extract`."""
    group = parser.add_argument_group("read clustering and selection")
    defaults = ExtractionConfig()
    group.add_argument("-s", "--min-identity", type=float, default=defaults.min_identity,
                       help="minimum overlap identity (matches/alignment length) to cluster two reads. "
                            "Depends on read accuracy: ~0.6 for R10.4.1, lower for R9")
    group.add_argument("-c", "--min-coverage", type=int, default=None,
                       help="minimum reads per cluster; default: automatic (largest coverage gap)")
    group.add_argument("-x", "--min-clusters", type=int, default=defaults.min_clusters,
                       help="automatic cut-off must keep more than this many clusters")
    group.add_argument("--keep-internal", action="store_true",
                       help="do not discard internal-match overlaps (shared repeats)")
    group.add_argument("--kmer", type=int, default=defaults.kmer_size, help="k for composition profiles")
    group.add_argument("--profile", choices=("members", "longest"), default=defaults.profile,
                       help="cluster composition from all member reads pooled (members) or from the longest "
                            "read only (longest, original MITNANEX)")
    group.add_argument("--group-rule", choices=("max", "median"), default=defaults.group_rule,
                       help="mitochondrial K-means group: the one holding the highest-coverage cluster (max) "
                            "or the one with the highest median coverage (median, original MITNANEX)")
    group.add_argument("-y", "--haplotype", action="store_true",
                       help="one composition group per cluster; keeps only the best cluster (haplotype separation)")
    group.add_argument("--plot-max-reads", type=int, default=defaults.plot_max_reads,
                       help="reads used in read-level composition plots")
    group.add_argument("--no-plots", action="store_true", help="skip the diagnostic figures")


def _extraction_config(args: argparse.Namespace, min_overlap_len: int | None = None) -> ExtractionConfig:
    return ExtractionConfig(
        min_identity=args.min_identity,
        min_overlap_len=min_overlap_len,
        drop_internal=not args.keep_internal,
        min_coverage=args.min_coverage,
        min_clusters=args.min_clusters,
        kmer_size=args.kmer,
        profile=args.profile,
        group_rule=args.group_rule,
        haplotype_mode=args.haplotype,
        plot_max_reads=args.plot_max_reads,
        make_plots=not args.no_plots,
    )


def _add_common_args(parser: argparse.ArgumentParser, default_threads: int) -> None:
    parser.add_argument("-i", "--reads", type=Path, required=True, help="input reads (FASTQ, optionally .gz)")
    parser.add_argument("-o", "-w", "--outdir", type=Path, default=Path("mitnanex_results"),
                        help="output folder")
    parser.add_argument("-r", "--prefix", help="prefix for output files [input file name]")
    parser.add_argument("-t", "--threads", type=int, default=default_threads, help="threads")
    parser.add_argument("-d", "--timestamp", action="store_true", help="append date and time to the output folder")
    parser.add_argument("--resume", action="store_true", help="skip steps whose outputs already exist")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mitnanex",
        description="MITNANEX - MITochondrial NANopore reads EXtractor",
        epilog="Authors: Juan Jose Picon Cossio, Javier Correa Alvarez, Gustavo Gamez de las Armas",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # --- denovo -----------------------------------------------------------------
    denovo = sub.add_parser("denovo", help="reference-free extraction + mitogenome assembly",
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    _add_common_args(denovo, default_threads=4)
    denovo.add_argument("-p", "--proportion", type=_proportion, default=0.4,
                        help="subsample for the overlap step: fraction or number of reads")
    denovo.add_argument("-m", "--min-len", type=int, default=None, help="minimum read length")
    denovo.add_argument("-M", "--max-len", type=int, default=None, help="maximum read length")
    denovo.add_argument("-K", "--minimap-kmer", type=int, default=15, help="minimap2 k for all-vs-all overlaps")
    denovo.add_argument("-S", "--miniasm-stage", type=int, default=7,
                        help="miniasm hidden -S option: stop graph cleaning after this stage (0 = miniasm default)")
    denovo.add_argument("-q", "--min-mapq", type=int, default=0, help="minimum MAPQ when recruiting reads")
    denovo.add_argument("-k", "--keep-percent", type=float, default=90, help="filtlong --keep_percent")
    denovo.add_argument("-f", "--flye-mode", default="--nano-hq", help="Flye read type flag")
    denovo.add_argument("--seed", type=int, default=11, help="subsampling seed")
    _add_extraction_args(denovo)

    # --- reference ----------------------------------------------------------------
    ref = sub.add_parser("reference", help="reference-guided extraction + variant calling",
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    _add_common_args(ref, default_threads=8)
    ref.add_argument("-R", "--reference", type=Path, required=True, help="reference mitogenome (single-contig FASTA)")
    ref.add_argument("-g", "--gff", type=Path, default=None, help="reference GFF3 for liftoff [bundled rCRS]")
    ref.add_argument("--mm2-preset", default="map-ont", help="minimap2 -x preset")
    ref.add_argument("--mapq", type=int, default=30, help="minimum mapping quality")
    ref.add_argument("--min-length", type=int, default=500, help="chopper minimum read length")
    ref.add_argument("--max-length", type=int, default=None, help="chopper maximum read length [reference length]")
    ref.add_argument("--min-mean-quality", type=float, default=10, help="chopper minimum mean read quality")
    ref.add_argument("--keep-percent", type=float, default=80, help="filtlong --keep_percent")
    ref.add_argument("--no-filter", action="store_true", help="do not filter reads with chopper")
    ref.add_argument("--min-pruning", type=int, default=3, help="GATK --min-pruning")
    ref.add_argument("--kmer-size", type=_int_list, default=[15, 25], help="GATK k-mer sizes, comma separated")
    ref.add_argument("--max-assembly-region-size", type=int, default=None, help="GATK --max-assembly-region-size")
    ref.add_argument("--flye-preset", default="--nano-hq", help="Flye read type flag")
    ref.add_argument("--other-flye-opts", default="", help="extra Flye options, quoted")
    ref.add_argument("--annotate-rcrs", action="store_true", help="annotate variants with rCRS regions (human)")
    ref.add_argument("--liftoff", action="store_true", help="lift the reference GFF onto the consensus")
    ref.add_argument("--haplogrep", type=Path, default=None, help="path to haplogrep3 to classify haplogroups")
    ref.add_argument("--trees", default="phylotree-rcrs@17.2", help="haplogrep3 trees, comma separated")
    ref.add_argument("--top-hits", type=int, default=3, help="haplogrep3 top hits")
    ref.add_argument("--mitomap", action="store_true", help="send PASS variants to the MITOMAP web service")

    # --- extract (core step only) -------------------------------------------------
    extract = sub.add_parser("extract", help="run only the read clustering/selection on existing reads + PAF",
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    extract.add_argument("-i", "--reads", type=Path, required=True, help="reads used to compute the overlaps (FASTQ)")
    extract.add_argument("-a", "--paf", type=Path, required=True, help="all-vs-all overlaps (minimap2 -x ava-ont)")
    extract.add_argument("-o", "--outdir", type=Path, default=Path("mitnanex_extract"), help="output folder")
    extract.add_argument("-r", "--prefix", default=None, help="prefix for output files [reads file name]")
    extract.add_argument("--min-overlap-len", type=int, default=None, help="minimum overlap length")
    extract.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    _add_extraction_args(extract)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "denovo":
            config = DenovoConfig(
                reads=args.reads, outdir=args.outdir, prefix=args.prefix, threads=args.threads,
                proportion=args.proportion, min_len=args.min_len, max_len=args.max_len,
                minimap_kmer=args.minimap_kmer, extraction=_extraction_config(args),
                miniasm_stage=args.miniasm_stage or None, min_mapq=args.min_mapq,
                keep_percent=args.keep_percent, flye_mode=args.flye_mode, seed=args.seed,
            )
            pipeline = DenovoPipeline(config, timestamped=args.timestamp, resume=args.resume)
        elif args.command == "reference":
            config = ReferenceConfig(
                reference=args.reference, reads=args.reads, gff=args.gff, outdir=args.outdir, prefix=args.prefix,
                threads=args.threads, minimap2_preset=args.mm2_preset, min_mapq=args.mapq,
                filter_reads=not args.no_filter, min_length=args.min_length, max_length=args.max_length,
                min_mean_quality=args.min_mean_quality, keep_percent=args.keep_percent,
                min_pruning=args.min_pruning, kmer_sizes=args.kmer_size,
                max_assembly_region_size=args.max_assembly_region_size, flye_preset=args.flye_preset,
                other_flye_opts=args.other_flye_opts.split(), annotate_rcrs=args.annotate_rcrs,
                liftoff=args.liftoff, haplogrep=args.haplogrep, haplogrep_trees=args.trees.split(","),
                top_hits=args.top_hits, mitomap=args.mitomap,
            )
            pipeline = ReferencePipeline(config, timestamped=args.timestamp, resume=args.resume)
        else:  # extract
            setup_logging(args.outdir / "mitnanex.log", args.verbose)
            prefix = args.prefix or args.reads.name.split(".")[0]
            MitoReadExtractor(_extraction_config(args, args.min_overlap_len)).run(args.reads, args.paf, args.outdir, prefix)
            return 0

        setup_logging(pipeline.outdir / "mitnanex.log", args.verbose)
        pipeline.run()
        return 0
    except (ToolError, RuntimeError, ValueError, FileNotFoundError) as error:
        log.error("%s", error)
        return 1


if __name__ == "__main__":
    sys.exit(main())
