#!/bin/bash
set -ev

# Create the various sets of VCF files with different population subsets.

# IN should be the directory containing the VCF files from
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
IN=1KG_high_coverage/

OUT=kgp_ooa3_i35
SAMPLE_FILE=samples.ooa3.i35.txt
mkdir -p ${OUT}
for i in $(seq 1 22); do
    BASE=1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel
    bcftools view --min-af 0.0000001:nref --max-af 0.9999999:nref -S ${SAMPLE_FILE} -v snps ${IN}/${BASE}.vcf.gz -O z -o ${OUT}/${BASE}.ooa3.vcf.gz
done

OUT=kgp_admix_pur_i35
SAMPLE_FILE=samples.admix_pur.i35.txt
mkdir -p ${OUT}
for i in $(seq 1 22); do
    BASE=1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel
    bcftools view --min-af 0.0000001:nref --max-af 0.9999999:nref -S ${SAMPLE_FILE} -v snps ${IN}/${BASE}.vcf.gz -O z -o ${OUT}/${BASE}.admix_pur.vcf.gz
done

OUT=kgp_admix_clm_i35
SAMPLE_FILE=samples.admix_clm.i35.txt
mkdir -p ${OUT}
for i in $(seq 1 22); do
    BASE=1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel
    bcftools view --min-af 0.0000001:nref --max-af 0.9999999:nref -S ${SAMPLE_FILE} -v snps ${IN}/${BASE}.vcf.gz -O z -o ${OUT}/${BASE}.admix_clm.vcf.gz
done

OUT=kgp_admix_mxl_i35
SAMPLE_FILE=samples.admix_mxl.i35.txt
mkdir -p ${OUT}
for i in $(seq 1 22); do
    BASE=1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel
    bcftools view --min-af 0.0000001:nref --max-af 0.9999999:nref -S ${SAMPLE_FILE} -v snps ${IN}/${BASE}.vcf.gz -O z -o ${OUT}/${BASE}.admix_mxl.vcf.gz
done

OUT=kgp_admix_clm_gbr_i35
SAMPLE_FILE=samples.admix_clm.gbr.i35.txt
mkdir -p ${OUT}
for i in $(seq 1 22); do
    BASE=1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel
    bcftools view --min-af 0.0000001:nref --max-af 0.9999999:nref -S ${SAMPLE_FILE} -v snps ${IN}/${BASE}.vcf.gz -O z -o ${OUT}/${BASE}.admix_clm_gbr.vcf.gz
done

OUT=kgp_admix_clm_ibs_i35
SAMPLE_FILE=samples.admix_clm.ibs.i35.txt
mkdir -p ${OUT}
for i in $(seq 1 22); do
    BASE=1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel
    bcftools view --min-af 0.0000001:nref --max-af 0.9999999:nref -S ${SAMPLE_FILE} -v snps ${IN}/${BASE}.vcf.gz -O z -o ${OUT}/${BASE}.admix_clm_ibs.vcf.gz
done

OUT=kgp_admix_pel_i35
SAMPLE_FILE=samples.admix_pel.i35.txt
mkdir -p ${OUT}
for i in $(seq 1 22); do
    BASE=1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel
    echo "Processing chromosome ${i}"
    bcftools view --min-af 0.0000001:nref --max-af 0.9999999:nref -S ${SAMPLE_FILE} -v snps ${IN}/${BASE}.vcf.gz | bcftools norm -m +any -O z -o ${OUT}/${BASE}.admix_pel.vcf.gz
    bcftools index ${OUT}/${BASE}.admix_pel.vcf.gz
done

OUT=kgp_ooa3fin_i35
SAMPLE_FILE=samples.ooa3fin.i35.txt
mkdir -p ${OUT}
for i in $(seq 1 22); do
    BASE=1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel
    echo "Processing chromosome ${i}"
    bcftools view --min-af 0.0000001:nref --max-af 0.9999999:nref -S ${SAMPLE_FILE} -v snps ${IN}/${BASE}.vcf.gz | bcftools norm -m +any -O z -o ${OUT}/${BASE}.ooa3fin.vcf.gz
    bcftools index ${OUT}/${BASE}.ooa3fin.vcf.gz
done
