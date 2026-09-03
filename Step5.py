#Visualizing and Final reporting

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def run_reporting(input_csv="ranked_cftr_variants_final.csv"):
    print("--- Visualizations & Final Reporting ---")

    # 1. Load Step 4 Ranked Data
    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} prioritized variants for plotting.")

    # High-contrast plot styling
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 14,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "figure.titlesize": 16,
        }
    )

    # ==========================================
    # VISUALIZATION 1: BAR CHART OF MUTATION TYPES
    # ==========================================
    print("Generating Plot 1: Mutation Types Bar Chart...")
    plt.figure(figsize=(10, 6))

    def classify_effect(effect):
        eff = str(effect).lower() if pd.notna(effect) else ""
        if "missense" in eff:
            return "Missense"
        elif "stop_gained" in eff or "nonsense" in eff:
            return "Nonsense / Stop Gained"
        elif "disruption" in eff or "deletion" in eff or "insertion" in eff:
            return "Indels / In-Frame"
        elif "utr" in eff:
            return "UTR Variant"
        elif "synonymous" in eff:
            return "Synonymous"
        else:
            return "Other / Splice / Intronic"

    df["SIMPLE_EFFECT"] = (
        df["EFFECT"].apply(classify_effect)
        if "EFFECT" in df.columns
        else "Unknown"
    )
    effect_counts = df["SIMPLE_EFFECT"].value_counts()

    ax1 = sns.barplot(
        x=effect_counts.index,
        y=effect_counts.values,
        hue=effect_counts.index,
        palette="Set2",
        edgecolor="black",
        linewidth=1.2,
    )
    plt.title(
        "Distribution of CFTR Mutation Types", fontsize=14, fontweight="bold"
    )
    plt.xlabel("Mutation Category", fontsize=12, fontweight="bold")
    plt.ylabel("Number of Variants", fontsize=12, fontweight="bold")
    plt.xticks(rotation=15, ha="right")

    for p in ax1.patches:
        height = p.get_height()
        if height > 0:
            ax1.annotate(
                f"{int(height)}",
                (p.get_x() + p.get_width() / 2.0, height),
                ha="center",
                va="bottom",
                xytext=(0, 4),
                textcoords="offset points",
                fontweight="bold",
                color="#2B2B2B",
            )

    plt.tight_layout()
    plt.savefig("1_mutation_types_barchart.png", dpi=300)
    plt.close()

    # ==========================================
    # VISUALIZATION 2: LOLLIPOP / DOMAIN CLUSTER PLOT (GRCh38)
    # ==========================================
    print("Generating Plot 2: CFTR Protein Domain Lollipop Plot (GRCh38)...")
    fig, ax2 = plt.subplots(figsize=(12, 6))

    # CFTR Genomic Coordinates Aligned to GRCh38 (NC_000007.14:g.117480025-117668665)
    domains = [
        ("TMD1 (Transmembrane 1)", 117480000, 117510000, "#2b5c8f"),
        ("NBD1 (Nucleotide Binding 1)", 117530000, 117560000, "#2a9d8f"),
        ("R Domain (Regulatory)", 117570000, 117590000, "#e76f51"),
        ("TMD2 (Transmembrane 2)", 117600000, 117630000, "#9b5de5"),
        ("NBD2 (Nucleotide Binding 2)", 117640000, 117660000, "#f15bb5"),
    ]

    for name, start, end, color in domains:
        ax2.add_patch(
            plt.Rectangle(
                (start, -0.5),
                end - start,
                1.0,
                color=color,
                alpha=0.75,
                ec="black",
                lw=1,
                label=name,
            )
        )

    if "POS" in df.columns and "RANK_SCORE" in df.columns:
        stemlines = ax2.stem(
            df["POS"],
            df["RANK_SCORE"],
            linefmt="#4A4A4A",
            markerfmt="o",
            basefmt=" ",
        )
        plt.setp(
            stemlines.markerline,
            color="#D90429",
            markersize=8,
            markeredgecolor="black",
            markeredgewidth=1,
        )

    ax2.set_title(
        "CFTR Variant Clustering Across Genomic Coordinates & Domains (GRCh38)",
        fontsize=14,
        fontweight="bold",
    )
    ax2.set_xlabel(
        "Genomic Position (chr7: GRCh38)", fontsize=12, fontweight="bold"
    )
    ax2.set_ylabel(
        "Pathogenicity Rank Score", fontsize=12, fontweight="bold"
    )
    ax2.set_ylim(-1, max(df["RANK_SCORE"].max() + 2, 10))
    ax2.legend(
        loc="upper left",
        title="CFTR Structural Domains",
        frameon=True,
        facecolor="white",
        edgecolor="black",
    )

    plt.tight_layout()
    plt.savefig("2_cftr_domain_lollipop.png", dpi=300)
    plt.close()

    # ==========================================
    # VISUALIZATION 3: CADD SCORE DISTRIBUTION
    # ==========================================
    print("Generating Plot 3: CADD Score Distribution...")
    plt.figure(figsize=(9, 5))

    cadd_scores = (
        pd.to_numeric(df["CADD_PHRED"], errors="coerce").dropna()
        if "CADD_PHRED" in df.columns
        else pd.Series(dtype=float)
    )

    if len(cadd_scores) > 0:
        sns.histplot(
            cadd_scores,
            kde=True,
            color="#2a9d8f",
            bins=15,
            edgecolor="black",
            linewidth=1.2,
            alpha=0.6,
        )
        plt.axvline(
            20,
            color="#D90429",
            linestyle="--",
            linewidth=2.5,
            label="Pathogenic Cutoff (CADD >= 20)",
        )
        plt.title(
            "Distribution of CADD PHRED Pathogenicity Scores",
            fontsize=14,
            fontweight="bold",
        )
        plt.xlabel("CADD PHRED Score", fontsize=12, fontweight="bold")
        plt.ylabel("Variant Count", fontsize=12, fontweight="bold")
        plt.legend(
            loc="upper left",
            frameon=True,
            facecolor="white",
            edgecolor="black",
        )
    else:
        plt.text(
            0.5,
            0.5,
            "No numerical CADD PHRED scores present in dataset",
            ha="center",
            va="center",
            fontsize=12,
        )

    plt.tight_layout()
    plt.savefig("3_cadd_score_distribution.png", dpi=300)
    plt.close()

    # ==========================================
    # FINAL EXPORT: SUMMARY CSV
    # ==========================================
    print("Exporting Final Prioritized Summary CSV...")

    cols_mapping = {
        "POS": "Position",
        "REF": "Ref",
        "ALT": "Alt",
        "HGVS_P": "HGVS Protein",
        "SIMPLE_EFFECT": "Effect Type",
        "CLNSIG": "ClinVar Significance",
        "CADD_PHRED": "CADD PHRED",
        "RANK_SCORE": "Priority Rank Score",
        "RECOMMENDED_THERAPY": "Targeted Therapy",
    }

    available_cols = [c for c in cols_mapping.keys() if c in df.columns]
    df_export = df[available_cols].copy()
    df_export.rename(columns=cols_mapping, inplace=True)

    df_export.to_csv("final_ranked_cftr_summary.csv", index=False)

    print("\n==========================================")
    print("FINAL STEP COMPLETE! Clean outputs generated:")
    print(" - 1_mutation_types_barchart.png")
    print(" - 2_cftr_domain_lollipop.png")
    print(" - 3_cadd_score_distribution.png")
    print(" - final_ranked_cftr_summary.csv")
    print("==========================================\n")


if __name__ == "__main__":
    run_reporting()
