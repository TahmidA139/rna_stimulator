# RNA Sequence Simulator

A Python command-line tool that generates synthetic RNA sequences in FASTA format. The simulator produces Open Reading Frames (ORFs) with configurable characteristics — complete or partial, optionally flanked by non-coding regions — and annotates each sequence with metadata including length, GC content, ambiguity content, type, and flanked status.

---

## Features

- Generates any number of synthetic RNA sequences in valid FASTA format
- Produces **complete ORFs** (AUG start codon + random non-stop codons + stop codon)
- Produces **partial ORFs** (random sequences without guaranteed start/stop codons)
- Randomly adds **5' and 3' flanking non-coding regions** based on a configurable probability
- Annotates each sequence with: length, GC%, ambiguity%, type, and flanked status
- Fully configurable via command-line arguments
- Modular codebase split into reusable library (`sequence_lib.py`) and simulation logic (`simulator.py`)

---

## Requirements

- Python 3.11 or later
- NumPy
- Biopython
- Anaconda or Miniconda (for environment management)

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-username/rna-simulator.git
cd rna-simulator
```

### 2. Create the Conda environment

```bash
conda env create -f environment.yml
```

### 3. Activate the environment

```bash
conda activate rna-simulator
```

### 4. Verify the installation

```bash
python src/main.py --help
```

---

## Usage

Run the simulator from the project root directory:

```bash
python src/main.py [OPTIONS]
```

### Example 1: Default parameters

```bash
python src/main.py
```

Generates 10 sequences with default settings and saves to `sequences.fasta`.

```
✓ Successfully generated 10 sequences
✓ Output saved to: sequences.fasta
```

### Example 2: Custom sequence count and output file

```bash
python src/main.py -n 50 -o my_sequences.fasta
```

### Example 3: Custom length range, completeness, and flanking

```bash
python src/main.py -n 20 --min-length 200 --max-length 2000 --completeness 0.9 --flanking-prob 0.3 --flanking-length 100 -o custom.fasta
```

---

## Command-Line Arguments

| Argument              | Short | Type  | Default           | Description                                            |
|-----------------------|-------|-------|-------------------|--------------------------------------------------------|
| `--num-sequences`     | `-n`  | int   | `10`              | Number of sequences to generate                        |
| `--output`            | `-o`  | str   | `sequences.fasta` | Output FASTA file path                                 |
| `--min-length`        |       | int   | `100`             | Minimum ORF length in nucleotides                      |
| `--max-length`        |       | int   | `1000`            | Maximum ORF length in nucleotides                      |
| `--flanking-prob`     |       | float | `0.5`             | Probability (0.0–1.0) that a sequence has flanking regions |
| `--flanking-length`   |       | int   | `50`              | Length of each flanking region in nucleotides          |
| `--completeness`      |       | float | `0.7`             | Fraction (0.0–1.0) of sequences that are complete ORFs |

---

## Output Format

Sequences are written in standard FASTA format, with 70 characters per line:

```
>seq_001 length=526 gc_content=50.95 ambiguity=0.0 type=complete flanked=yes
GGGУUGGA...
>seq_002 length=420 gc_content=49.29 ambiguity=0.0 type=partial flanked=no
AAGGGCUA...
```

### Header fields

| Field         | Description                                                  |
|---------------|--------------------------------------------------------------|
| `length`      | Total length of the sequence in nucleotides                  |
| `gc_content`  | Percentage of G and C nucleotides (0.0–100.0)                |
| `ambiguity`   | Percentage of IUPAC ambiguity code nucleotides (0.0–100.0)   |
| `type`        | `complete` (has AUG + stop codon) or `partial` (random)      |
| `flanked`     | `yes` if non-coding flanking regions were added, else `no`   |

---

## Project Structure

```
rna-simulator/
├── README.md                    # This file
├── LICENSE                      # MIT License
├── pseudocode.txt               # Algorithm planning document
├── environment.yml              # Conda environment specification
├── src/
│   ├── __init__.py              # Package metadata (version, author)
│   ├── main.py                  # CLI entry point (argparse)
│   ├── simulator.py             # Simulator class
│   └── sequence_lib.py          # Reusable sequence utility functions
└── examples/
    ├── example_output.fasta     # Sample generated output
    └── example_run.txt          # Example commands and their output
```

---

## Algorithm Description

### Sequence Generation

1. For each of the `num_sequences` sequences:
   - Randomly decide if the sequence is a **complete** ORF (based on `completeness_ratio`)
   - **Complete ORF**: Build `AUG` + random body codons (stop codons excluded from body) + random stop codon (`UAA`, `UAG`, or `UGA`)
   - **Partial ORF**: Generate a fully random RNA sequence of random length in `[min_orf_length, max_orf_length]`
   - Randomly decide if flanking regions are added (based on `flanking_probability`)
   - If flanked: prepend a random 5' region and append a random 3' region, each of `flanking_length` nucleotides

### Metadata Calculation

- **GC Content**: Count of G and C nucleotides divided by total sequence length, multiplied by 100
- **Ambiguity Content**: Count of IUPAC ambiguity code characters (N, R, Y, S, W, K, M, B, D, H, V) divided by total length, multiplied by 100

---

## References

- IUPAC Nucleotide Codes: https://www.bioinformatics.org/sms/iupac.html
- FASTA Format: https://en.wikipedia.org/wiki/FASTA_format
- Open Reading Frames: https://en.wikipedia.org/wiki/Open_reading_frame
- Biopython Documentation: https://biopython.org/
- Python argparse: https://docs.python.org/3/library/argparse.html
- Python Type Hints: https://docs.python.org/3/library/typing.html
- PEP 8 Style Guide: https://www.python.org/dev/peps/pep-0008/

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Author

Your Name  
your.email@example.com
