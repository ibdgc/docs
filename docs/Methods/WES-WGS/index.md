# Methods
# For the vcf files from folders “SubsetHailJointCall”. These files didn’t have required FILTER tags in the header. Also, the GT field was named as LGT, which prevented accurate processing of these files with bcftools and plink. An example for how to correct the issues is shown below
 
ml bcftools
#Did the following steps for vcfs under “SubsetHailJointCall” folders to correct headers
bcftools view -h /sc/arion/projects/ibdgc/uchicago-transfer/data/wes/updates/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes_SubsetHailJointCall/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes.vcf.bgz > hdr.txt
##add necessary filter flags in the hdr.txt:
##FILTER=<ID=VQSRTrancheSNP99.00to99.90+,Description="VQSRTrancheSNP99.00to99.90+">
##FILTER=<ID=VQSRTrancheINDEL99.00to99.90+,Description="VQSRTrancheINDEL99.00to99.90+">
##FILTER=<ID=VQSRTrancheINDEL99.00to99.90,Description="VQSRTrancheINDEL99.00to99.90">
##FILTER=<ID=VQSRTrancheSNP99.00to99.90,Description="VQSRTrancheSNP99.00to99.90">
 
 
bcftools reheader -h hdr.txt -o header_fixed_anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes.vcf.bgz /sc/arion/projects/ibdgc/uchicago-transfer/data/wes/updates/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes_SubsetHailJointCall/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes.vcf.bgz

#Fix the GT tag in the FORMAT
bcftools view -Ov header_fixed_anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes.vcf.bgz | sed -e 's/LGT/GT/g' | bcftools view -Oz -o GT_fixed_anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes.vcf.bgz
