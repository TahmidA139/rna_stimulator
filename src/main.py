#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py - Entry point for the RNA Sequence Simulator.

Parses command-line arguments, validates them, creates a Simulator instance,
and saves the generated sequences to a FASTA file.
"""

import sys
import argparse

from simulator import Simulator


def create_parser() -> argparse.ArgumentParser:
    """
    Create and configure the command-line argument parser.

    Returns:
        argparse.ArgumentParser: Configured parser with all required arguments.
    """
    parser = argparse.ArgumentParser(
        prog="rna-simulator",
        description="Generate synthetic RNA sequences in FASTA format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--num-sequences", "-n",
        type=int,
        default=10,
        help="Number of sequences to generate.",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default="sequences.fasta",
        help="Output FASTA file path.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=100,
        help="Minimum ORF length in nucleotides.",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=1000,
        help="Maximum ORF length in nucleotides.",
    )
    parser.add_argument(
        "--flanking-prob",
        type=float,
        default=0.5,
        help="Probability (0.0–1.0) that a sequence has flanking regions.",
    )
    parser.add_argument(
        "--flanking-length",
        type=int,
        default=50,
        help="Length of flanking regions in nucleotides.",
    )
    parser.add_argument(
        "--completeness",
        type=float,
        default=0.7,
        help="Ratio (0.0–1.0) of complete ORFs vs partial ORFs.",
    )
    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    """
    Validate parsed command-line arguments for logical consistency.

    Args:
        args (argparse.Namespace): Parsed arguments from argparse.

    Returns:
        None

    Raises:
        ValueError: If any argument fails a validation check.
    """
    if args.num_sequences < 1:
        raise ValueError("--num-sequences must be >= 1.")
    if args.min_length < 3:
        raise ValueError("--min-length must be >= 3.")
    if args.max_length < args.min_length:
        raise ValueError("--max-length must be >= --min-length.")
    if not (0.0 <= args.flanking_prob <= 1.0):
        raise ValueError("--flanking-prob must be between 0.0 and 1.0.")
    if args.flanking_length < 1:
        raise ValueError("--flanking-length must be >= 1.")
    if not (0.0 <= args.completeness <= 1.0):
        raise ValueError("--completeness must be between 0.0 and 1.0.")


def main() -> int:
    """
    Main entry point for the RNA sequence simulator.

    Parses and validates arguments, creates a Simulator, generates sequences,
    and saves them to the specified FASTA output file.

    Returns:
        int: 0 on success, 1 on error.
    """
    parser = create_parser()
    args = parser.parse_args()

    try:
        validate_arguments(args)
    except ValueError as e:
        sys.stderr.write(f"ERROR: {e}\n")
        return 1

    try:
        simulator = Simulator(
            num_sequences=args.num_sequences,
            min_orf_length=args.min_length,
            max_orf_length=args.max_length,
            flanking_probability=args.flanking_prob,
            flanking_length=args.flanking_length,
            completeness_ratio=args.completeness,
        )
        simulator.save_fasta(args.output)
    except (ValueError, IOError) as e:
        sys.stderr.write(f"ERROR: {e}\n")
        return 1

    sys.stdout.write(f"✓ Successfully generated {args.num_sequences} sequences\n")
    sys.stdout.write(f"✓ Output saved to: {args.output}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
