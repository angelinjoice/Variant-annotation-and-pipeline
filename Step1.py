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
