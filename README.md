# Arginine Codon Usage Analysis

A Python tool for computing arginine codon frequencies in coronaviruses using NCBI GenBank annotations. This script correctly handles complex genomic features including:

- Separate ORFs with their own reading frames
- `codon_start` offsets
- Joins and ribosomal frameshifts (e.g., ORF1ab)
- Overlapping ORFs (e.g., in HCoV-OC43)

## Features

- Fetches genome data directly from NCBI GenBank
- Analyzes codon usage across all CDS features
- Generates formatted tables and visualizations
- Supports multiple coronavirus genomes:
  - SARS-CoV-2 (Wuhan-Hu-1)
  - Bat coronavirus RaTG13
  - Human coronavirus OC43

## Requirements

- Python 3.7+
- See `requirements.txt` for dependencies

## Installation

1. Clone or download this repository
2. Install dependencies:

```bash
pip install -r requirements.txt
```

## Usage

Run the script:

```bash
python arg_codon_table.py
```

The script will:
1. Fetch genome data from NCBI for each configured virus
2. Extract and analyze CDS sequences
3. Calculate arginine codon frequencies
4. Print a formatted table to the console
5. Generate visualization plots (PNG and PDF) in the `outputs/` directory

## Output

- **Console**: A formatted table showing arginine codon usage percentages for each virus
- **Plots**: Grouped bar charts saved as:
  - `outputs/arginine_codon_usage.png`
  - `outputs/arginine_codon_usage.pdf`

## Arginine Codons Analyzed

The script analyzes all six arginine codons:
- CGT
- CGC
- CGA
- CGG
- AGA
- AGG

## Configuration

To analyze different viruses or accessions, modify the `VIRUS_ACCESSIONS` dictionary in `arg_codon_table.py`:

```python
VIRUS_ACCESSIONS: Dict[str, str] = {
    "SARS-CoV-2": "MN908947.3",
    "RaTG13": "MN996532.1",
    "human CoV OC43": "AY585228.1",
}
```

## Notes

- The script requires an internet connection to fetch data from NCBI
- GenBank records are fetched using NCBI EFetch API
- CDS sequences are extracted in their correct reading frames, accounting for `codon_start` offsets
- The visualization highlights CGG codon usage with special formatting

## License

This project is provided as-is for research and educational purposes.

