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
