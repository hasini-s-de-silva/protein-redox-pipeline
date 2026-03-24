from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable


REDOX_PATTERNS = {
    "redox_potential": [
        re.compile(r"redox\s+potential\s*[:=]\s*([-+]?\d*\.?\d+)", re.IGNORECASE),
        re.compile(r"E0\s*[:=]\s*([-+]?\d*\.?\d+)", re.IGNORECASE),
        re.compile(r"E°\s*[:=]\s*([-+]?\d*\.?\d+)", re.IGNORECASE),
    ],
    "reorganization_energy": [
        re.compile(
            r"reorganization\s+free\s+energy\s*[:=]\s*([-+]?\d*\.?\d+)",
            re.IGNORECASE,
        ),
        re.compile(r"lambda\s*[:=]\s*([-+]?\d*\.?\d+)", re.IGNORECASE),
    ],
    "electronic_coupling": [
        re.compile(
            r"electronic\s+coupling\s*[:=]\s*([-+]?\d*\.?\d+)",
            re.IGNORECASE,
        ),
        re.compile(r"Hab\s*[:=]\s*([-+]?\d*\.?\d+)", re.IGNORECASE),
    ],
}


def iter_text_files(run_dir: Path) -> Iterable[Path]:
    for ext in ("*.txt", "*.log", "*.out", "*.dat", "*.csv"):
        yield from run_dir.glob(ext)


def extract_first_match(text: str, patterns: list[re.Pattern[str]]) -> str | None:
    for pattern in patterns:
        match = pattern.search(text)
        if match:
            return match.group(1)
    return None


def parse_run(run_dir: Path) -> dict[str, str]:
    record = {
        "sample": run_dir.name,
        "redox_potential": "",
        "reorganization_energy": "",
        "electronic_coupling": "",
        "source_file": "",
    }

    for file_path in iter_text_files(run_dir):
        try:
            text = file_path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue

        found_any = False
        for field, patterns in REDOX_PATTERNS.items():
            if not record[field]:
                value = extract_first_match(text, patterns)
                if value is not None:
                    record[field] = value
                    found_any = True

        if found_any and not record["source_file"]:
            record["source_file"] = str(file_path)

        if all(record[field] for field in ("redox_potential", "reorganization_energy", "electronic_coupling")):
            break

    return record


def analyse_results(results_root: Path, output_csv: Path) -> None:
    run_dirs = [p for p in results_root.iterdir() if p.is_dir()]
    if not run_dirs:
        raise FileNotFoundError(f"No run directories found in: {results_root}")

    rows = [parse_run(run_dir) for run_dir in sorted(run_dirs)]

    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "sample",
                "redox_potential",
                "reorganization_energy",
                "electronic_coupling",
                "source_file",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[INFO] Analysis complete. Wrote: {output_csv}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Parse BioDC result folders into a single CSV."
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        required=True,
        help="Path to a single runs_YYYYMMDD_HHMMSS folder",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/aggregated_results.csv"),
        help="Path to output CSV file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    analyse_results(
        results_root=args.results_root.resolve(),
        output_csv=args.output_csv.resolve(),
    )


if __name__ == "__main__":
    main()