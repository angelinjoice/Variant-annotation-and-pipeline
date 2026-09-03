
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
