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
