"""Common machinery for MITNANEX pipelines: output folder, logging, steps."""

from __future__ import annotations

import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Callable, Sequence

from ..external import check_tools

log = logging.getLogger("mitnanex")

BANNER = r"""
  __  __ ___ _____ _   _    _    _   _ _______  __
 |  \/  |_ _|_   _| \ | |  / \  | \ | | ____\ \/ /
 | |\/| || |  | | |  \| | / _ \ |  \| |  _|  \  /
 | |  | || |  | | | |\  |/ ___ \| |\  | |___ /  \
 |_|  |_|___| |_| |_| \_/_/   \_\_| \_|_____/_/\_\
"""


def setup_logging(log_file: Path | None = None, verbose: bool = False) -> None:
    """Log to the terminal and, optionally, to a file in the output folder.

    Every message carries its own timestamp (the bash version computed the
    timestamp once, so every message showed the start time).
    """
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )


def default_prefix(reads: Path) -> str:
    """File name up to the first dot (``sample.fastq.gz`` -> ``sample``)."""
    return reads.name.split(".")[0]


class Pipeline:
    """Base class: sets up the output folder and runs named steps in order.

    Subclasses define :attr:`required_tools` and :meth:`steps`.

    Args:
        outdir: Output folder (created if needed).
        prefix: Prefix of every output file.
        timestamped: Append date and time to ``outdir`` (old ``-d`` option).
        resume: Skip steps whose declared outputs already exist.
    """

    name = "pipeline"
    required_tools: Sequence[str] = ()

    def __init__(self, outdir: Path, prefix: str, timestamped: bool = False, resume: bool = False) -> None:
        if timestamped:
            outdir = outdir.with_name(f"{outdir.name}_{datetime.now():%Y-%m-%d_%H-%M-%S}")
        self.outdir = Path(outdir)
        self.prefix = prefix
        self.resume = resume
        self.tool_log = self.outdir / "logs" / "tools.log"  # stderr of every external tool
        self.summary: dict[str, object] = {}

    def path(self, suffix: str) -> Path:
        """Output path ``<outdir>/<prefix><suffix>``."""
        return self.outdir / f"{self.prefix}{suffix}"

    def steps(self) -> list[tuple[str, Callable[[], None], list[Path]]]:
        """Return ``(description, function, outputs)`` for each step, in order."""
        raise NotImplementedError

    def run(self) -> dict[str, object]:
        """Run every step; returns the summary that is also written as JSON."""
        self.outdir.mkdir(parents=True, exist_ok=True)
        log.info(BANNER)
        log.info("%s pipeline | prefix: %s | output: %s", self.name, self.prefix, self.outdir.resolve())
        check_tools(self.required_tools)

        start = time.perf_counter()
        for i, (description, function, outputs) in enumerate(self.steps(), start=1):
            if self.resume and outputs and all(p.exists() for p in outputs):
                log.info("Step %d: %s -> outputs exist, skipped (--resume)", i, description)
                continue
            log.info("Step %d: %s", i, description)
            step_start = time.perf_counter()
            function()
            log.info("Step %d done in %.1f s", i, time.perf_counter() - step_start)

        self.summary["elapsed_seconds"] = round(time.perf_counter() - start, 1)
        summary_file = self.path("_summary.json")
        summary_file.write_text(json.dumps(self.summary, indent=2, default=str))
        log.info("### MITNANEX finished in %.1f s ###", self.summary["elapsed_seconds"])
        log.info("Summary: %s", summary_file)
        return self.summary
