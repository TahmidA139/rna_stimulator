# -*- coding: utf-8 -*-
"""
simulator.py - Main Simulator class for generating synthetic RNA sequences.

Orchestrates the generation of complete and partial ORFs with optional
flanking regions, and produces FASTA output with per-sequence metadata.
"""

import random
from typing import List, Tuple

from sequence_lib import (
    generate_random_codon,
    generate_random_sequence,
    get_gc_content,
    get_ambiguity_content,
    write_fasta,
    is_stop_codon,
    STOP_CODONS,
    START_CODON,
)


class Simulator:
    """
    Generates synthetic RNA sequences with configurable characteristics.

    Produces a mix of complete ORFs (AUG...stop codon) and partial ORFs,
    optionally flanked by non-coding regions. Each sequence is annotated
    with metadata including length, GC content, ambiguity content, type,
    and flanked status.

    Attributes:
        num_sequences (int): Number of sequences to generate.
        min_orf_length (int): Minimum ORF length in nucleotides.
        max_orf_length (int): Maximum ORF length in nucleotides.
        flanking_probability (float): Probability that a sequence has flanking regions.
        flanking_length (int): Length of each flanking region in nucleotides.
        completeness_ratio (float): Proportion of sequences that are complete ORFs.
    """

    def __init__(
        self,
        num_sequences: int = 10,
        min_orf_length: int = 100,
        max_orf_length: int = 1000,
        flanking_probability: float = 0.5,
        flanking_length: int = 50,
        completeness_ratio: float = 0.7,
    ) -> None:
        """
        Initialize the Simulator with sequence generation parameters.

        Args:
            num_sequences (int): Number of sequences to generate (>= 1).
            min_orf_length (int): Minimum ORF length in nucleotides (>= 3).
            max_orf_length (int): Maximum ORF length (>= min_orf_length).
            flanking_probability (float): Probability of flanking regions (0.0–1.0).
            flanking_length (int): Length of flanking regions in nucleotides (>= 1).
            completeness_ratio (float): Fraction of complete ORFs (0.0–1.0).

        Raises:
            ValueError: If any parameter fails validation.
        """
        if num_sequences < 1:
            raise ValueError(f"num_sequences must be >= 1, got {num_sequences}.")
        if min_orf_length < 3:
            raise ValueError(f"min_orf_length must be >= 3, got {min_orf_length}.")
        if max_orf_length < min_orf_length:
            raise ValueError(
                f"max_orf_length ({max_orf_length}) must be >= min_orf_length ({min_orf_length})."
            )
        if not (0.0 <= flanking_probability <= 1.0):
            raise ValueError(f"flanking_probability must be 0.0–1.0, got {flanking_probability}.")
        if flanking_length < 1:
            raise ValueError(f"flanking_length must be >= 1, got {flanking_length}.")
        if not (0.0 <= completeness_ratio <= 1.0):
            raise ValueError(f"completeness_ratio must be 0.0–1.0, got {completeness_ratio}.")

        self.num_sequences = num_sequences
        self.min_orf_length = min_orf_length
        self.max_orf_length = max_orf_length
        self.flanking_probability = flanking_probability
        self.flanking_length = flanking_length
        self.completeness_ratio = completeness_ratio

    def generate_orf(self, complete: bool = True) -> str:
        """
        Generate a single ORF sequence.

        If complete=True, generates AUG + random non-stop codons + stop codon.
        If complete=False, generates a random RNA sequence of random length
        within [min_orf_length, max_orf_length] without guaranteed start/stop.

        Args:
            complete (bool): If True, produce a complete ORF; otherwise partial.

        Returns:
            str: The generated ORF sequence as an uppercase RNA string.
        """
        target_length = random.randint(self.min_orf_length, self.max_orf_length)

        if not complete:
            return generate_random_sequence(target_length)

        # Complete ORF: AUG + body codons + stop codon
        # Reserve 6 nt for AUG + stop; remainder is body (rounded to codon boundary)
        body_length = max(0, target_length - 6)
        num_body_codons = body_length // 3

        body_codons: List[str] = []
        for _ in range(num_body_codons):
            codon = generate_random_codon()
            # Regenerate if a stop codon appears in the body
            while is_stop_codon(codon):
                codon = generate_random_codon()
            body_codons.append(codon)

        stop_codon = random.choice(STOP_CODONS)
        return START_CODON + ''.join(body_codons) + stop_codon

    def generate_sequence(self) -> Tuple[str, bool, bool]:
        """
        Generate a complete sequence with randomized features.

        Decides whether this sequence is a complete or partial ORF based on
        completeness_ratio, then optionally adds 5' and/or 3' flanking regions.

        Returns:
            Tuple[str, bool, bool]: (sequence, is_complete, is_flanked)
                - sequence: the full RNA sequence string
                - is_complete: True if this is a complete ORF
                - is_flanked: True if flanking regions were added
        """
        is_complete = random.random() < self.completeness_ratio
        orf = self.generate_orf(complete=is_complete)

        is_flanked = random.random() < self.flanking_probability
        if is_flanked:
            five_prime = generate_random_sequence(self.flanking_length)
            three_prime = generate_random_sequence(self.flanking_length)
            sequence = five_prime + orf + three_prime
        else:
            sequence = orf

        return sequence, is_complete, is_flanked

    def generate_sequences(self) -> List[Tuple[str, str, str]]:
        """
        Generate all sequences with associated metadata.

        Creates num_sequences sequences, each with a unique identifier and
        a description string containing length, GC content, ambiguity
        content, type (complete/partial), and flanked status.

        Returns:
            List[Tuple[str, str, str]]: List of (sequence_id, description, sequence)
                tuples ready for FASTA output.
        """
        results: List[Tuple[str, str, str]] = []

        for i in range(1, self.num_sequences + 1):
            seq_id = f"seq_{i:03d}"
            sequence, is_complete, is_flanked = self.generate_sequence()

            gc = get_gc_content(sequence)
            ambiguity = get_ambiguity_content(sequence)
            seq_type = "complete" if is_complete else "partial"
            flanked = "yes" if is_flanked else "no"

            description = (
                f"length={len(sequence)} "
                f"gc_content={gc} "
                f"ambiguity={ambiguity} "
                f"type={seq_type} "
                f"flanked={flanked}"
            )
            results.append((seq_id, description, sequence))

        return results

    def save_fasta(self, output_file: str) -> None:
        """
        Generate all sequences and save them to a FASTA file.

        Calls generate_sequences() to produce the sequences and then
        writes them to the specified file using write_fasta() from
        sequence_lib.

        Args:
            output_file (str): Path to the output FASTA file.

        Returns:
            None

        Raises:
            IOError: If the file cannot be written.
        """
        sequences = self.generate_sequences()
        write_fasta(sequences, output_file)
