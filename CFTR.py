#STEP 1

#Extracting vcf  file
from cyvcf2 import VCF, Writer

input_vcf = "clinvar.vcf.gz"
output_vcf = "cftr_clinvar.vcf"
REGION = "7:117465784-117682387"

vcf = VCF(input_vcf)
out_vcf = Writer(output_vcf, vcf)

for variant in vcf(REGION):
    out_vcf.write_record(variant)

out_vcf.close()
vcf.close()
print("Extracted cftr_clinvar.vcf successfully")

#STEP 2 

#create a test vcf file
##zgrep "^#" clinvar.vcf.gz > cftr_test.vcf                    #copied header from original file to test file

#installation of snpEff
#install snpEff from https://pcingola.github.io/SnpEff/        #upload on codespace
##unzip snpEff_latest_core.zip
##grep "GRCh38" snpEff/snpEff.config                           # finds database
##java -jar snpEff/snpEff.jar download GRCh38.86                #downloads database

#Annotating test vcf file
##java -Xmx2g -jar snpEff/snpEff.jar -v -csvStats testresults GRCh38.86 cftr_test.vcf > cftr_testannotated.vcf

#Annotating original vcf file
##java -Xmx2g -jar snpEff/snpEff.jar -v -csvStats summary GRCh38.86 cftr_clinvar.vcf > cftr_clinvar_annotated.vcf

#STEP 3

#Data base integration (cross referncing) the original file
#parsing the data
import re
import numpy as np
import pandas as pd


def build_dataframe(vcf_path, output_csv="scored_annotations.csv"):
    print(f"building: Processing VCF from '{vcf_path}'...")

    # Known high-confidence CFTR2 variant mapping lookup table
    cftr2_known = {
        'Phe508del': 'CF-causing',
        'F508del': 'CF-causing',
        'Gly551Asp': 'CF-causing',
        'G551D': 'CF-causing',
        'Gly542Ter': 'CF-causing',
        'G542X': 'CF-causing',
        'Arg117H': 'Varying_Clinical_Consequence',
        'R117H': 'Varying_Clinical_Consequence',
        'Ile507del': 'CF-causing',
        'I507del': 'CF-causing',
    }

    records = []

    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]

            clnsig_match = re.search(r'CLNSIG=([^;]+)', info)
            clnsig_val = (
                clnsig_match.group(1) if clnsig_match else 'Not_Provided'
            )

            clndn_match = re.search(r'CLNDN=([^;]+)', info)
            clndn_val = clndn_match.group(1) if clndn_match else 'Not_Provided'

            af_match = re.search(r'AF(?:_EXAC|_TGP)?=([0-9.eE-]+)', info)
            af_val = float(af_match.group(1)) if af_match else 0.0

            ann_match = re.search(r'ANN=([^;]+)', info)
            if ann_match:
                first_ann = ann_match.group(1).split(',')[0].split('|')
                effect = first_ann[1] if len(first_ann) > 1 else 'Unknown'
                impact = first_ann[2] if len(first_ann) > 2 else 'Unknown'
                hgvs_p = first_ann[10] if len(first_ann) > 10 else ''
            else:
                effect, impact, hgvs_p = 'None', 'None', ''

            cftr2_val = 'Uncategorized_in_CFTR2'
            for mut, status in cftr2_known.items():
                if mut in hgvs_p:
                    cftr2_val = status
                    break

            cadd_match = re.search(r'CADD_PHRED=([0-9.]+)', info)
            sift_match = re.search(r'SIFT_pred=([ADBN])', info)
            polyphen_match = re.search(r'POLYPHEN_pred=([DPBN])', info)

            cadd_val = float(cadd_match.group(1)) if cadd_match else np.nan
            sift_val = sift_match.group(1) if sift_match else 'N/A'
            polyphen_val = polyphen_match.group(1) if polyphen_match else 'N/A'

            records.append({
                'CHROM': chrom,
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'CLNSIG': clnsig_val,
                'CLNDN': clndn_val,
                'AF': af_val,
                'EFFECT': effect,
                'IMPACT': impact,
                'HGVS_P': hgvs_p,
                'CFTR2_STATUS': cftr2_val,
                'CADD_PHRED': cadd_val,
                'SIFT': sift_val,
                'POLYPHEN': polyphen_val,
            })

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    return df


vcf_input_file = 'cftr_clinvar_annotated.vcf'
df_step3 = build_dataframe(vcf_input_file)

#STEP 4

#Scoring for clinical prioritization

import numpy as np
import pandas as pd


def filtering_and_scoring_vectorized(
    input_csv="scored_annotations.csv",
    output_csv="ranked_cftr_variants_final.csv",
):
    print("--- Starting Step 4 (Vectorized) ---")
    print(f"Loading data from '{input_csv}'...")

    df = pd.read_csv(input_csv)
    initial_count = len(df)

    # 1. RARE VARIANT FILTERING (AF < 0.01)
    # Ensure numeric type safely
    df["AF"] = pd.to_numeric(df.get("AF"), errors="coerce").fillna(0.0)
    df_filtered = df[df["AF"] < 0.01].copy()
    filtered_count = len(df_filtered)

    print(
        f"Filtered out {initial_count - filtered_count} common variants (AF >= 0.01)."
    )
    print(f"Remaining rare variants to score: {filtered_count}")

    # 2. SNPEFF IMPACT NUMERICAL MAPPING
    impact_map = {
        "HIGH": 4,
        "MODERATE": 3,
        "LOW": 2,
        "MODIFIER": 1,
        "UNKNOWN": 0,
        "NONE": 0,
    }
    impact_str = (
        df_filtered["IMPACT"].fillna("").astype(str).str.upper()
        if "IMPACT" in df_filtered.columns
        else pd.Series("", index=df_filtered.index)
    )
    df_filtered["IMPACT_SCORE"] = (
        impact_str.map(impact_map).fillna(0).astype(int)
    )

    # 3. VECTORIZED PATHOGENICITY RANK SCORE CALCULATION
    # Extract & normalize strings safely across entire Series
    clnsig_series = (
        df_filtered["CLNSIG"].fillna("").astype(str).str.lower()
        if "CLNSIG" in df_filtered.columns
        else pd.Series("", index=df_filtered.index)
    )

    cadd_series = (
        pd.to_numeric(df_filtered["CADD_PHRED"], errors="coerce").fillna(0.0)
        if "CADD_PHRED" in df_filtered.columns
        else pd.Series(0.0, index=df_filtered.index)
    )

    cftr2_series = (
        df_filtered["CFTR2_STATUS"].fillna("").astype(str)
        if "CFTR2_STATUS" in df_filtered.columns
        else pd.Series("", index=df_filtered.index)
    )

    sift_series = (
        df_filtered["SIFT"].fillna("").astype(str).str.upper()
        if "SIFT" in df_filtered.columns
        else pd.Series("", index=df_filtered.index)
    )

    polyphen_series = (
        df_filtered["POLYPHEN"].fillna("").astype(str).str.upper()
        if "POLYPHEN" in df_filtered.columns
        else pd.Series("", index=df_filtered.index)
    )

    # Evaluate boolean vector conditions
    is_pathogenic = clnsig_series.str.contains(
        "pathogenic", regex=False
    ) & ~clnsig_series.str.contains("benign", regex=False)
    is_cftr2_causing = cftr2_series == "CF-causing"
    is_high_cadd = cadd_series >= 20.0
    is_deleterious = (sift_series == "D") | (polyphen_series == "D")

    # Sum bonuses directly as integer arrays
    df_filtered["RANK_SCORE"] = (
        df_filtered["IMPACT_SCORE"]
        + np.where(is_pathogenic, 3, 0)
        + np.where(is_cftr2_causing, 3, 0)
        + np.where(is_high_cadd, 2, 0)
        + np.where(is_deleterious, 1, 0)
    )

    # 4. VECTORIZED PRECISION DRUG THERAPY MAPPING
    hgvs_series = (
        df_filtered["HGVS_P"].fillna("").astype(str)
        if "HGVS_P" in df_filtered.columns
        else pd.Series("", index=df_filtered.index)
    )

    effect_series = (
        df_filtered["EFFECT"].fillna("").astype(str).str.lower()
        if "EFFECT" in df_filtered.columns
        else pd.Series("", index=df_filtered.index)
    )

    # Condition 1: Trikafta candidates
    cond_trikafta = hgvs_series.str.contains(
        "F508del|Phe508del", regex=True
    )

    # Condition 2: Kalydeco candidates
    cond_kalydeco = hgvs_series.str.contains(
        "G551D|Gly551Asp|R117H|Arg117His", regex=True
    )

    # Condition 3: Read-through candidates
    cond_readthrough = (
        hgvs_series.str.contains("G542X|Gly542Ter", regex=True)
        | effect_series.str.contains("stop_gained", regex=False)
        | effect_series.str.contains("nonsense", regex=False)
    )

    # Condition 4: General CFTR modulators
    cond_general = df_filtered["IMPACT_SCORE"] >= 3

    # Define precedence conditions and matching output labels
    therapy_conditions = [
        cond_trikafta,
        cond_kalydeco,
        cond_readthrough,
        cond_general,
    ]

    therapy_choices = [
        "Elexacaftor/Tezacaftor/Ivacaftor (Trikafta)",
        "Ivacaftor (Kalydeco)",
        "Read-Through Candidate (Experimental / mRNA)",
        "General CFTR Modulator Candidate",
    ]

    # Select therapy based on prioritized array evaluation
    df_filtered["RECOMMENDED_THERAPY"] = np.select(
        therapy_conditions,
        therapy_choices,
        default="No Specific Targeted Modulator",
    )

    # 5. SORTING AND PRIORITIZATION
    df_sorted = df_filtered.sort_values(
        by=["RANK_SCORE", "CADD_PHRED", "IMPACT_SCORE"],
        ascending=[False, False, False],
    )

    # Export final table
    df_sorted.to_csv(output_csv, index=False)
    print(f"\nStep 4 Complete! Output saved to '{output_csv}'.")

    # Display Top Ranked Variants Preview
    print("\n--- TOP RANKED VARIANTS PREVIEW ---")
    preview_cols = [
        c
        for c in [
            "POS",
            "REF",
            "ALT",
            "HGVS_P",
            "CLNSIG",
            "RANK_SCORE",
            "RECOMMENDED_THERAPY",
        ]
        if c in df_sorted.columns
    ]
    print(df_sorted[preview_cols].head())

    return df_sorted


if __name__ == "__main__":
    filtering_and_scoring_vectorized("scored_annotations.csv")

#STEP 5

#Visualizing and final reporting

##pip install matplotlib
##pip install seaborn

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
