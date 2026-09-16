#!/usr/bin/env python3
"""Convert semicolon-separated residue labels to '+'-joined residue numbers."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


RESIDUE_NUMBER = re.compile(r"-?\d+[A-Za-z]?")


def residue_number(label: str, line_number: int) -> str:
    """Return the numeric residue identifier from a label such as 'I115'."""
    label = label.strip()
    match = RESIDUE_NUMBER.search(label)
    if match is None:
        raise ValueError(
            f"Line {line_number}: cannot find a residue number in {label!r}"
        )
    return match.group(0)


def convert(input_path: Path, output_path: Path) -> None:
    with input_path.open("r", encoding="utf-8-sig", newline="") as source:
        lines = source.readlines()

    if not lines:
        raise ValueError(f"Input file is empty: {input_path}")

    converted = [lines[0]]  # Keep the first line exactly as it is.
    for line_number, line in enumerate(lines[1:], start=2):
        newline = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
        content = line.rstrip("\r\n")

        if not content.strip():
            converted.append(line)
            continue

        fields = content.split(";")
        first_column = fields[0]  # Keep the first column exactly as it is.
        residue_labels = [field for field in fields[1:] if field.strip()]
        numbers = "+".join(
            residue_number(label, line_number) for label in residue_labels
        )
        converted.append(f"{first_column};{numbers}{newline}")

    with output_path.open("w", encoding="utf-8", newline="") as destination:
        destination.writelines(converted)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Keep the header and first column intact, and replace semicolon-separated "
            "residue labels in later columns with residue numbers joined by '+'."
        )
    )
    parser.add_argument("input_csv", type=Path, help="CSV file to convert")
    parser.add_argument(
        "-o", "--output", type=Path, help="Output file (default: INPUT_residue_numbers.csv)"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = args.input_csv
    output_path = args.output or input_path.with_name(
        f"{input_path.stem}_residue_numbers{input_path.suffix}"
    )

    if input_path.resolve() == output_path.resolve():
        raise ValueError("Input and output paths must be different")

    convert(input_path, output_path)
    print(f"Created: {output_path}")


if __name__ == "__main__":
    main()
