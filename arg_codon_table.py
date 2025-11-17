#!/usr/bin/env python3
"""
Compute arginine codon frequencies for coronaviruses using NCBI GenBank
annotations (CDS features), correctly handling:

- Separate ORFs with their own reading frames
- codon_start offsets
- joins / ribosomal frameshifts (e.g. ORF1ab)
- overlapping ORFs (e.g. in HCoV-OC43)

Genomes / accessions:

- SARS-CoV-2 Wuhan-Hu-1:     MN908947.3
- Bat coronavirus RaTG13:    MN996532.1
- Human coronavirus OC43:    AY585228.1
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from io import StringIO
from typing import Dict, Iterable, List

import os
import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

VIRUS_ACCESSIONS: Dict[str, str] = {
    "SARS-CoV-2": "MN908947.3",      # Wuhan-Hu-1
    "RaTG13": "MN996532.1",
    "human CoV OC43": "AY585228.1",
}

ARG_CODONS: List[str] = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]


@dataclass
class VirusResult:
    name: str
    accession: str
    arg_counts: Dict[str, int]
    arg_freqs: Dict[str, float]  # fraction of arg codons (0–1)


# ---------------------------------------------------------------------
# Fetch and parse GenBank
# ---------------------------------------------------------------------


def fetch_genbank_record(accession: str) -> SeqRecord:
    """
    Fetch a GenBank record for the given accession using NCBI EFetch.
    """
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gb",
        "retmode": "text",
    }
    resp = requests.get(base, params=params, timeout=60)
    resp.raise_for_status()
    handle = StringIO(resp.text)
    record = SeqIO.read(handle, "genbank")
    return record


# ---------------------------------------------------------------------
# Extract CDS sequences in correct frame
# ---------------------------------------------------------------------


def iter_cds_nuc_seqs(record: SeqRecord) -> Iterable[str]:
    """
    Yield nucleotide sequences for each CDS in the correct reading frame.

    - Uses the feature location (handles join/complement/etc.)
    - Applies /codon_start (1, 2, or 3) if present
    - Trims any trailing bases not forming a complete codon
    """
    for feat in record.features:
        if feat.type != "CDS":
            continue

        cds_seq = feat.extract(record.seq)
        cds_seq = str(cds_seq).upper().replace("-", "")

        codon_start_str = feat.qualifiers.get("codon_start", ["1"])[0]
        codon_start = int(codon_start_str)

        if codon_start not in (1, 2, 3):
            raise ValueError(
                f"Unexpected codon_start={codon_start} in {record.id}"
            )

        if codon_start > 1:
            cds_seq = cds_seq[codon_start - 1 :]

        usable_len = (len(cds_seq) // 3) * 3
        cds_seq = cds_seq[:usable_len]

        yield str(cds_seq)


# ---------------------------------------------------------------------
# Codon counting
# ---------------------------------------------------------------------


def count_codons_from_cds(record: SeqRecord) -> Counter:
    """
    Count codons across all CDS features of a GenBank record, in
    their correct reading frames.
    """
    counts: Counter = Counter()
    for cds in iter_cds_nuc_seqs(record):
        for i in range(0, len(cds), 3):
            codon = cds[i : i + 3]
            if len(codon) == 3 and all(b in "ACGT" for b in codon):
                counts[codon] += 1
    return counts


def arginine_counts_and_freqs(counts: Counter) -> Dict[str, float]:
    """
    Compute arginine codon frequencies from codon counts.
    Returns dict mapping codon -> frequency (0–1) among arginine codons.
    """
    arg_counts = {c: int(counts.get(c, 0)) for c in ARG_CODONS}
    total_arg = sum(arg_counts.values())
    if total_arg == 0:
        raise ValueError("No arginine codons found in CDS set.")
    return {c: arg_counts[c] / total_arg for c in ARG_CODONS}


# ---------------------------------------------------------------------
# Table printing
# ---------------------------------------------------------------------


def format_and_print_table(results: List[VirusResult]) -> None:
    """
    Print a side-by-side table of codon frequencies.
    """
    name_col_width = max(len(r.name) for r in results) + 2
    codon_col_width = 6
    value_col_width = 8

    # Header row
    header = " " * codon_col_width
    for r in results:
        header += r.name.center(value_col_width)
    print(header)

    # Rows: one per arginine codon
    for codon in ARG_CODONS:
        row = codon.ljust(codon_col_width)
        for r in results:
            pct = r.arg_freqs[codon] * 100.0
            row += f"{pct:5.1f}%".rjust(value_col_width)
        print(row)


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------


def plot_arginine_codon_usage(
    results: List[VirusResult], outfile: str | None = None
) -> str | None:
    """
    Create a grouped bar chart of arginine codon usage (percent among arginine
    codons) across viruses.
    """
    try:
        import numpy as np  # type: ignore
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as import_exc:
        raise RuntimeError(
            "matplotlib/numpy not available; install to enable plotting."
        ) from import_exc

    viruses = [r.name for r in results]
    codons = ARG_CODONS

    # Matrix: rows=viruses, cols=codons, values=percentages
    data = np.array(
        [
            [results[i].arg_freqs[c] * 100.0 for c in codons]
            for i in range(len(results))
        ]
    )

    num_codons = len(codons)
    num_viruses = len(viruses)
    group_positions = np.arange(num_codons)
    total_group_width = 0.8
    bar_width = total_group_width / num_viruses

    fig, ax = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.7, alpha=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", which="both", direction="out")

    greys = np.linspace(0.25, 0.70, num_viruses) if num_viruses > 1 else [0.25]
    hatches = ["", "//", "xx", "\\\\", "++", ".."]

    # Highlight CGG
    try:
        cgg_idx = codons.index("CGG")
    except ValueError:
        cgg_idx = None
    if cgg_idx is not None:
        left = cgg_idx - (total_group_width / 2) - (bar_width * 0.15)
        right = cgg_idx + (total_group_width / 2) + (bar_width * 0.15)
        ax.axvspan(left, right, color="0.93", zorder=0)

    for vidx, virus in enumerate(viruses):
        offsets = group_positions - (total_group_width / 2) + (vidx + 0.5) * bar_width
        bars = ax.bar(
            offsets,
            data[vidx],
            width=bar_width,
            label=virus,
            color=str(greys[vidx]),
            edgecolor="black",
            linewidth=0.8,
            hatch=hatches[vidx % len(hatches)],
        )
        # Highlight CGG bars
        if cgg_idx is not None:
            try:
                bar_cgg = bars[cgg_idx]
                emph_grey = float(np.clip(greys[vidx] - 0.15, 0.05, 0.85))
                bar_cgg.set_facecolor(str(emph_grey))
                bar_cgg.set_edgecolor("black")
                bar_cgg.set_linewidth(1.4)
                height = bar_cgg.get_height()
                ax.text(
                    bar_cgg.get_x() + bar_cgg.get_width() / 2,
                    height + max(0.5, height * 0.015),
                    f"{height:.1f}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color="black",
                )
            except Exception:
                pass

    ax.set_xticks(group_positions)
    ax.set_xticklabels(codons, fontsize=10)
    ax.set_ylabel("Percent among arginine codons (%)", fontsize=11)
    ax.set_title("Arginine codon usage (CDS)", fontsize=12)
    ax.legend(
        title="Virus",
        fontsize=9,
        title_fontsize=10,
        frameon=False,
        ncols=min(3, num_viruses),
        handlelength=1.6,
        columnspacing=1.2,
        handletextpad=0.6,
        borderaxespad=0.2,
        loc="upper left",
    )

    # Save plot
    if outfile is None:
        out_dir = os.path.join(os.path.dirname(__file__), "outputs")
        os.makedirs(out_dir, exist_ok=True)
        outfile = os.path.join(out_dir, "arginine_codon_usage.png")

    fig.savefig(outfile, dpi=200)
    try:
        fig.savefig(os.path.splitext(outfile)[0] + ".pdf")
    except Exception:
        pass
    plt.close(fig)
    return outfile


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main() -> None:
    results: List[VirusResult] = []

    for name, acc in VIRUS_ACCESSIONS.items():
        print(f"Fetching {name} ({acc}) from NCBI...", flush=True)
        record = fetch_genbank_record(acc)
        codon_counts = count_codons_from_cds(record)
        arg_freqs = arginine_counts_and_freqs(codon_counts)
        results.append(
            VirusResult(
                name=name,
                accession=acc,
                arg_counts={c: codon_counts.get(c, 0) for c in ARG_CODONS},
                arg_freqs=arg_freqs,
            )
        )

    print()
    format_and_print_table(results)
    try:
        out_path = plot_arginine_codon_usage(results)
        if out_path:
            print(f"\nSaved plot to {out_path}")
    except Exception as exc:
        print(f"\n[plot skipped] {exc}")


if __name__ == "__main__":
    main()
