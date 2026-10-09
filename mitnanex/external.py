"""Thin, logged wrappers around the external command-line tools.

Commands are passed as argument lists (never through a shell string), so
paths with spaces are safe and there is no shell-injection risk. Pipes such as
``minimap2 ... | samtools view ...`` are built with :func:`run_pipeline`,
which connects the processes directly and fails if *any* step fails (the
original bash scripts silently ignored errors in the middle of a pipe).
"""

from __future__ import annotations

import logging
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Sequence

log = logging.getLogger(__name__)

Command = Sequence[str | Path | int | float]


class ToolError(RuntimeError):
    """Raised when an external tool is missing or exits with an error."""


def check_tools(tools: Sequence[str]) -> None:
    """Fail early, with a helpful message, if any required executable is missing."""
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        raise ToolError(
            f"Missing executables: {', '.join(missing)}. "
            "Install them with: mamba env create -f environment.yml && mamba activate mitnanex"
        )


def _as_args(command: Command) -> list[str]:
    """Convert every element of a command to ``str``."""
    return [str(arg) for arg in command]


def _log_command(commands: Sequence[Command], stdout: Path | None) -> None:
    text = " | ".join(shlex.join(_as_args(c)) for c in commands)
    log.info("$ %s%s", text, f" > {stdout}" if stdout else "")


def _tail(path: Path, n_bytes: int = 3000) -> str:
    """Last ``n_bytes`` of a (possibly large) text file."""
    with open(path, "rb") as handle:
        handle.seek(max(path.stat().st_size - n_bytes, 0))
        return handle.read().decode(errors="replace")


def run(command: Command, stdout: str | Path | None = None, log_file: str | Path | None = None,
        cwd: str | Path | None = None) -> None:
    """Run one command.

    Args:
        command: Program and arguments.
        stdout: Write the standard output to this file.
        log_file: Append the standard error to this file (otherwise it is
            only shown if the command fails).
        cwd: Working directory for the command.
    """
    run_pipeline([command], stdout=stdout, log_file=log_file, cwd=cwd)


def run_pipeline(commands: Sequence[Command], stdout: str | Path | None = None,
                 log_file: str | Path | None = None, cwd: str | Path | None = None) -> None:
    """Run ``cmd1 | cmd2 | ... > stdout`` and check every exit status.

    Standard error never goes to a pipe: with pipes, a verbose downstream
    tool could fill its stderr buffer while we wait for an upstream one,
    dead-locking the pipeline. With ``log_file`` every process appends its
    stderr to that file *as it runs* (so progress of long tools such as Flye
    or GATK can be followed with ``tail -f``); otherwise temporary files are used.
    """
    stdout = Path(stdout) if stdout else None
    _log_command(commands, stdout)
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, "a") as handle:
            handle.write(f"\n### {' | '.join(shlex.join(_as_args(c)) for c in commands)}\n")
        # One append-mode handle per process: O_APPEND keeps concurrent writes intact.
        stderr_files = [open(log_file, "ab") for _ in commands]
    else:
        stderr_files = [tempfile.TemporaryFile() for _ in commands]

    processes: list[subprocess.Popen] = []
    try:
        with open(stdout, "wb") if stdout else open(os.devnull, "wb") as out_handle:
            previous_stdout = None
            for i, command in enumerate(commands):
                last = i == len(commands) - 1
                process = subprocess.Popen(
                    _as_args(command),
                    stdin=previous_stdout,
                    stdout=out_handle if last else subprocess.PIPE,
                    stderr=stderr_files[i],
                    cwd=cwd,
                )
                # Close our copy of the pipe so the upstream process gets SIGPIPE if needed.
                if previous_stdout is not None:
                    previous_stdout.close()
                previous_stdout = process.stdout
                processes.append(process)
            for process in processes:
                process.wait()
    finally:
        # Collect stderr text for error messages before closing the handles.
        stderr_texts = []
        for handle in stderr_files:
            if log_file:
                stderr_texts.append(None)
            else:
                handle.seek(0)
                stderr_texts.append(handle.read().decode(errors="replace"))
            handle.close()

    failed = [(c, p, t) for c, p, t in zip(commands, processes, stderr_texts) if p.returncode != 0]
    if failed:
        # With a log file the stderr of all steps is interleaved there: show its tail.
        log_tail = _tail(log_file) if log_file else ""
        messages = [f"{shlex.join(_as_args(c))} exited with {p.returncode}" + (f":\n{t[-2000:]}" if t else "")
                    for c, p, t in failed]
        raise ToolError("\n".join(messages) + (f"\n--- end of {log_file} ---\n{log_tail}" if log_tail else ""))


def capture(command: Command) -> str:
    """Run a command and return its standard output as text."""
    log.debug("$ %s", shlex.join(_as_args(command)))
    result = subprocess.run(_as_args(command), capture_output=True, text=True)
    if result.returncode != 0:
        raise ToolError(f"{shlex.join(_as_args(command))} exited with {result.returncode}:\n{result.stderr[-2000:]}")
    return result.stdout
