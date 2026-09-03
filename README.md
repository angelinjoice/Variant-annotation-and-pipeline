# Automated Computational Profiling of CFTR Gene Variants for Pathogenicity Assessment in Cystic Fibrosis

A bioinformatics pipeline designed to parse, annotate, filter, and prioritize *CFTR* (Cystic Fibrosis Transmembrane Conductance Regulator) genomic variants from VCF datasets for clinical pathogenicity assessment.

## Objective

The primary objective of this project is to build an automated, reproducible computational pipeline that extracts *CFTR* gene variants from high-throughput genomic data, enriches them with multi-database annotations (ClinVar, gnomAD, CFTR2), predicts functional impact using *in silico* scores (CADD, SIFT, PolyPhen-2), and prioritizes high-confidence pathogenic mutations for clinical interpretation and targeted therapeutic mapping.

## Genomic Domain & Coordinates
* **Target Gene:** *CFTR* (Cystic Fibrosis Transmembrane Conductance Regulator)
* **Organism:** *Homo sapiens*
* **Chromosome Region:** `chr7:117,480,025 - 117,668,665`
* **Genome Assembly:** GRCh38 / hg38 (`NC_000007.14`)

## Workflow & Methodology

The pipeline follows a structured 5-step computational workflow:

### Step 1: Data Acquisition
* Coordinates: Human *CFTR* gene on Chromosome 7 (GRCh38: `NC_000007.14:g.117480025-117668665`).
* Input: ClinVar VCF datasets and sample target files containing known *CFTR* variants.

### Step 2: Variant Annotation
* Tools: SnpEff integration.
* Added Context: Gene locations (exon, intron, UTRs), mutation types (missense, nonsense, frameshift), and protein-level changes.

### Step 3: Database Integration (Cross referencing)
* **ClinVar**: Extract reported clinical significance (Pathogenic, Benign, VUS).
* **CFTR2**: Cross-reference variants with disease-specific clinical data.
* **gnomAD**: Retrieve global Allele Frequencies ($AF$) to filter common polymorphisms ($AF < 0.01$).
* **In Silico Predictors**: Fetch SIFT, PolyPhen-2, and CADD impact scores.

### Step 4: Scoring Python Data Processing & Prioritization
Executed via custom Python modules (`scripts/parsing.py`, `scripts/pipeline.py`):
1. **Quality Filtering**: Retain high-confidence variant calls.
2. **Frequency Filtering**: Exclude common benign variants ($AF > 0.01$).
3. **Pathogenicity Prioritization**: Flag deleterious mutations co-predicted by SIFT/PolyPhen-2.
4. **Therapeutic Mapping**: Map actionable variants to known CFTR modulators (e.g., Trikafta/Ivacaftor).

### Step 5: Visualization & Final Reporting
Executed via `scripts/visualiz.py`:
* Bar charts of mutation distributions (Missense vs. Nonsense vs. Indels).
* Protein domain mapping (Lollipop plots mapping Nucleotide Binding and Transmembrane Domains).
* CADD pathogenicity score distribution plots.
* Exported final ranked CSV output (`final_ranked_cftr_summary.csv`).

---

## Repository Structure

```text
├── CfTr.py                   # Master wrapper script executing end-to-end pipeline
├── scripts/                  # Modular pipeline core scripts
│   ├── parsing.py            # VCF reading and variant extraction module
│   ├── pipeline.py           # Annotation parsing, scoring, and filtering script
│   └── visualiz.py           # Plotting and summary visualization module
├── cftr_test.vcf             # Sample VCF dataset for testing
├── requirements.txt          # Python dependencies
├── images/                   # Output visualization plots
│   ├── 1_mutation_types_barchart.png
│   ├── 2_cftr_domain_lollipop.png
│   └── 3_cadd_score_distribution.png
└── README.md                 # Project documentation
