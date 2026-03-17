# -*- coding: utf-8 -*-
"""
sequence_lib.py - Reusable library functions for RNA sequence operations.

Provides utilities for GC content calculation, ambiguity detection,
random codon/sequence generation, and FASTA file writing.
"""

import random
from typing import List, Tuple

# Valid RNA nucleotides and IUPAC ambiguity codes
RNA_NUCLEOTIDES: List[str] = ['A', 'U', 'G', 'C']
IUPAC_AMBIGUITY_CODES: List[str] = [
    'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'
]
START_CODON: str = 'AUG'
STOP_CODONS: List[str] = ['UAA', 'UAG', 'UGA']


def get_gc_content(sequence: str) -> float:
    """
    Calculate the GC content of an RNA sequence.

    Args:
        sequence (str): RNA sequence string (uppercase, U for uracil).

    Returns:
        float: GC content as a percentage (0.0 to 100.0).

    Raises:
        ValueError: If sequence is empty.

    Example:
        >>> get_gc_content("AUGC")
        50.0
    """
    if not sequence:
        raise ValueError("Sequence must not be empty.")
    gc_count = sum(1 for nt in sequence.upper() if nt in ('G', 'C'))
    return round((gc_count / len(sequence)) * 100, 2)


def get_ambiguity_content(sequence: str) -> float:
    """
    Calculate the percentage of IUPAC ambiguity codes in an RNA sequence.

    Args:
        sequence (str): RNA sequence string (uppercase).

    Returns:
        float: Ambiguity content as a percentage (0.0 to 100.0).

    Raises:
        ValueError: If sequence is empty.

    Example:
        >>> get_ambiguity_content("AUGNNN")
        50.0
    """
    if not sequence:
        raise ValueError("Sequence must not be empty.")
    ambig_count = sum(1 for nt in sequence.upper() if nt in IUPAC_AMBIGUITY_CODES)
    return round((ambig_count / len(sequence)) * 100, 2)


def generate_random_codon() -> str:
    """
    Generate a random 3-nucleotide RNA codon.

    Returns:
        str: A random codon consisting of A, U, G, C nucleotides.

    Example:
        >>> codon = generate_random_codon()
        >>> len(codon) == 3
        True
    """
    return ''.join(random.choices(RNA_NUCLEOTIDES, k=3))


def is_start_codon(codon: str) -> bool:
    """
    Check whether a given codon is a start codon (AUG).

    Args:
        codon (str): A 3-nucleotide RNA codon string.

    Returns:
        bool: True if the codon is AUG, False otherwise.

    Example:
        >>> is_start_codon("AUG")
        True
        >>> is_start_codon("UAA")
        False
    """
    return codon.upper() == START_CODON


def is_stop_codon(codon: str) -> bool:
    """
    Check whether a given codon is a stop codon (UAA, UAG, or UGA).

    Args:
        codon (str): A 3-nucleotide RNA codon string.

    Returns:
        bool: True if the codon is UAA, UAG, or UGA, False otherwise.

    Example:
        >>> is_stop_codon("UAA")
        True
        >>> is_stop_codon("AUG")
        False
    """
    return codon.upper() in STOP_CODONS


def generate_random_sequence(length: int) -> str:
    """
    Generate a random RNA sequence of the specified length.

    Args:
        length (int): Desired length of the random sequence (must be >= 1).

    Returns:
        str: A random RNA sequence of the given length using A, U, G, C.

    Raises:
        ValueError: If length is less than 1.

    Example:
        >>> seq = generate_random_sequence(10)
        >>> len(seq) == 10
        True
    """
    if length < 1:
        raise ValueError(f"Sequence length must be >= 1, got {length}.")
    return ''.join(random.choices(RNA_NUCLEOTIDES, k=length))


def write_fasta(sequences: List[Tuple[str, str, str]], output_file: str) -> None:
    """
    Write a list of sequences to a file in FASTA format.

    Args:
        sequences (List[Tuple[str, str, str]]): List of tuples, each containing
            (sequence_id, description, sequence).
        output_file (str): Path to the output FASTA file.

    Returns:
        None

    Raises:
        ValueError: If sequences list is empty.
        IOError: If the file cannot be written.

    Example:
        >>> write_fasta([("seq1", "length=10 gc_content=50.0", "AUGCAUGCAU")], "out.fasta")
    """
    if not sequences:
        raise ValueError("Sequences list must not be empty.")
    with open(output_file, 'w', encoding='utf-8') as fh:
        for seq_id, description, sequence in sequences:
            fh.write(f">{seq_id} {description}\n")
            # Write sequence in 70-character lines (standard FASTA)
            for i in range(0, len(sequence), 70):
                fh.write(sequence[i:i + 70] + '\n')
