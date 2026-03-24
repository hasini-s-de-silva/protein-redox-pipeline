from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable


def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def safe_name(path: Path) -> str:
    return path.stem.replace(" ", "_").replace("/", "_")


def discover_inputs(input_dir: Path, extensions: tuple[str, ...]) -> list[Path]:
    files: list[Path] = []
    for ext in extensions:
        files.extend(sorted(input_dir.glob(f"*{ext}")))
    return sorted(set(files))


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def copy_template_input(template_input: Path, destination: Path) -> None:
    if not template_input.exists():
        raise FileNotFoundError(f"Template input file not found: {template_input}")
    shutil.copy2(template_input, destination)


def prepare_run_directory(
    base_run_dir: Path,
    input_file: Path,
    template_input: Path,
) -> Path:
    run_dir = base_run_dir / safe_name(input_file)
    ensure_dir(run_dir)

    # Copy the protein/input file into the run directory
    copied_input = run_dir / input_file.name
    shutil.copy2(input_file, copied_input)

    # Copy BioDC batch input template into the run directory
    copied_template = run_dir / "input.txt"
    copy_template_input(template_input, copied_template)

    return run_dir


def run_biodc(
    biodc_script: Path,
    run_dir: Path,
    python_executable: str = sys.executable,
    timeout: int = 3600,
) -> subprocess.CompletedProcess[str]:
    """
    Runs BioDC in the run directory.

    You may need to adjust the command below depending on how your local BioDC
    setup expects batch execution to be triggered.
    """
    command = [python_executable, str(biodc_script)]

    result = subprocess.run(
        command,
        cwd=run_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    return result


def write_run_log(
    log_path: Path,
    sample_name: str,
    command_status: int,
    stdout: str,
    stderr: str,
) -> None:
    with log_path.open("w", encoding="utf-8") as f:
        f.write(f"sample: {sample_name}\n")
        f.write(f"return_code: {command_status}\n")
        f.write("\n=== STDOUT ===\n")
        f.write(stdout or "")
        f.write("\n\n=== STDERR ===\n")
        f.write(stderr or "")


def append_summary_row(
    summary_csv: Path,
    row: dict[str, str | int],
) -> None:
    file_exists = summary_csv.exists()
    fieldnames = ["sample", "input_file", "run_dir", "return_code", "status"]

    with summary_csv.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


def pipeline(
    biodc_script: Path,
    input_dir: Path,
    template_input: Path,
    output_root: Path,
    extensions: tuple[str, ...],
    timeout: int,
) -> None:
    inputs = discover_inputs(input_dir, extensions)
    if not inputs:
        raise FileNotFoundError(
            f"No input files found in {input_dir} with extensions: {extensions}"
        )

    run_root = output_root / f"runs_{timestamp()}"
    ensure_dir(run_root)

    summary_csv = run_root / "pipeline_summary.csv"

    print(f"Found {len(inputs)} input files.")
    print(f"Run root: {run_root}")

    for input_file in inputs:
        sample = safe_name(input_file)
        print(f"\n[INFO] Running sample: {sample}")

        run_dir = prepare_run_directory(
            base_run_dir=run_root,
            input_file=input_file,
            template_input=template_input,
        )

        try:
            result = run_biodc(
                biodc_script=biodc_script,
                run_dir=run_dir,
                timeout=timeout,
            )
            status = "success" if result.returncode == 0 else "failed"

            write_run_log(
                log_path=run_dir / "run.log",
                sample_name=sample,
                command_status=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
            )

            append_summary_row(
                summary_csv=summary_csv,
                row={
                    "sample": sample,
                    "input_file": input_file.name,
                    "run_dir": str(run_dir),
                    "return_code": result.returncode,
                    "status": status,
                },
            )

            print(f"[INFO] Completed {sample} with status: {status}")

        except Exception as exc:
            write_run_log(
                log_path=run_dir / "run.log",
                sample_name=sample,
                command_status=-1,
                stdout="",
                stderr=str(exc),
            )

            append_summary_row(
                summary_csv=summary_csv,
                row={
                    "sample": sample,
                    "input_file": input_file.name,
                    "run_dir": str(run_dir),
                    "return_code": -1,
                    "status": f"error: {exc}",
                },
            )

            print(f"[ERROR] {sample}: {exc}")

    print("\n[INFO] Pipeline finished.")
    print(f"[INFO] Summary written to: {summary_csv}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run BioDC across multiple protein inputs."
    )
    parser.add_argument(
        "--biodc-script",
        type=Path,
        default=Path("biodc/V2.2/BioDCv2.py"),
        help="Path to BioDCv2.py",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("examples"),
        help="Directory containing input structures/files",
    )
    parser.add_argument(
        "--template-input",
        type=Path,
        default=Path("configs/input.txt"),
        help="Template BioDC input.txt file for batch runs",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("results"),
        help="Directory where run folders will be created",
    )
    parser.add_argument(
        "--extensions",
        nargs="+",
        default=[".pdb", ".inp", ".txt"],
        help="Input file extensions to include",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=3600,
        help="Timeout in seconds for each BioDC run",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pipeline(
        biodc_script=args.biodc_script.resolve(),
        input_dir=args.input_dir.resolve(),
        template_input=args.template_input.resolve(),
        output_root=args.output_root.resolve(),
        extensions=tuple(args.extensions),
        timeout=args.timeout,
    )


if __name__ == "__main__":
    main()