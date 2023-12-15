# Deduplication
#KING pipeline – December 2023
 
#1. European GSA - merged, imputed plink format

#Set unique names for each sample by adding the dataset name
awk '{print " gsa_allCHRmerged_imputed ""@"$2"\t"dataset"@"$2"\t"$3"\t"$4"\t"$5"\t"$6}' gsa_allCHRmerged_imputed.fam> tmp.fam
mv tmp.fam gsa_allCHRmerged_imputed.fam
 
#2. non-European GSA - chromosome-based imputed genotypes in vcf format
 
ml plink2/v2.00a3.3
ml plink
all_chroms=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
 
#Set dataset name (changed for each site)
dataset=niddk_cho_gsa
#Generated plink2 format pgen, pvar, psam files for each chromosome
for i in ${all_chroms[@]}
do
    	plink2 --vcf $dataset\_noneur_topmed_2022/$i.dose.vcf.gz --make-pfile --out $i
    	echo $i >> merge_list.txt
done
#Concatenated chromosome-based files.
plink2 --pmerge-list merge_list.txt --out merged
#Generated bim,bed,fam files.
plink2 --pfile merged --make-bed --out $dataset
rm chr* merge*
#Set unique names for each sample by adding the dataset name.
awk -v dataset="$dataset" '{print dataset"@"$2"\t"dataset"@"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $dataset.fam > tmp.fam
mv tmp.fam $dataset.fam
 
#Selected overlapping variant IDs and generate a file containing those variant IDs named NONEURGSAIDs.txt.
#Extracted overlapping variants from each site-specific plink file.
file_list=$(ls | grep bed$)
for i in $file_list
do
    	plink2 --bfile ${i%????} --make-bed --extract NONEURGSAIDs.txt --out selected_${i%????}
    	echo selected_${i%????} >> merge_list.txt
done
#Merged all site-specific plink files together (only available in plink v.1.9).
plink --merge-list merge_list.txt --out nonEUR_GSA_merged --keep-allele-order
 
#3. African American Array – chromosome-based imputed genotypes in vcf format (Omni and Axiom)
 
ml plink2/v2.00a3.3
ml bcftools
file_list=$(ls | grep vcf.gz$)
 
for file in $file_list
do
#VCF headers didn’t contain the required FILTER tag, “GENOTYPED”.
#Added the line: “##FILTER=<ID=GENOTYPED,Description="Site was genotyped">” in the header.
    	bcftools view -h $file > hdr.txt
    	sed -i '1,2d' hdr.txt
    	cat first.txt hdr.txt > newhdr.txt
    	bcftools reheader -h newhdr.txt -o head_corrected_$file $file
#Normalized, split multi-allelic variants and corrected the swapped reference-alternate alleles.
    	bcftools norm -c s -m -any head_corrected_$file -Oz -o norm_$file --fasta-ref ../Homo_sapiens.GRCh38.dna.toplevel.fa
    	name="${file%???????}"
#Generated plink2 format pgen, pvar, psam files for each chromosome.
    	plink2 --vcf norm_$file --make-pfile --out $name --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 1000
done
wait
#Concatenated chromosome-based files.
plink2 --pmerge-list omni_merge_list.txt --out mergedOmni
plink2 --pmerge-list axiom_merge_list.txt --out mergedAxiom
#Generated bim,bed,fam files.
plink2 --pfile mergedOmni --make-bed --out mergedOmni
 
#Set unique names for each sample by adding the dataset name.
awk '{print dataset"@"$2"\t"dataset"@"$2"\t"$3"\t"$4"\t"$5"\t"$6}' mergedOmni.fam > tmp.fam
mv tmp.fam mergedOmni.fam
plink2 --pfile mergedAxiom --make-bed --out mergedAxiom
awk '{print dataset"@"$2"\t"dataset"@"$2"\t"$3"\t"$4"\t"$5"\t"$6}' mergedAxiom.fam > tmp.fam
mv tmp.fam mergedAxiom.fam
 
#4. WGS – plink format

#Set unique names for each sample by adding the dataset name.
awk '{print "cutler-AA-WGS-big_daly""@"$2"\t"dataset"@"$2"\t"$3"\t"$4"\t"$5"\t"$6}' cutler-AA-WGS-big_daly.fam> tmp.fam
mv tmp.fam cutler-AA-WGS-big_daly.fam
 
#5. WES – vcf files
 
# I. For the 275 small vcf files from stampfer: /sc/arion/projects/ibdgc/uchicago-transfer/data/wes/updates/anvil_ccdg_broad_ai_ibd_daly_stampfer_wes/anvil_ccdg_broad_ai_ibd_daly_stampfer_wes/sharded_vcf/
#filelist.txt contains the names of all 275 files in each new line.
bcftools concat -f filelist.txt -Oz -o ibd_daly_stampfer_wes.filtered.vcf.gz
 
# II. For the vcf files that don’t have a problem in the header or GT field (and ibd_daly_stampfer_wes.filtered.vcf.gz).
ml plink2/v2.00a3.3
ml bcftools
file_list=$(ls | grep vcf.gz$)
 
for file in $file_list
do
#Normalized and split multi allelic variants.
    	bcftools norm -c s -m -any $file -Oz -o norm_$file --fasta-ref ../Homo_sapiens.GRCh38.dna.toplevel.fa
#Selected variant that passed the QC check that was already done.
    	bcftools view -f PASS norm_$file -Oz -o PASS_$file
    	name="${file%???????}"
#Converted into plink format and set variant ids.
    	plink2 --vcf PASS_$file --make-bed --out $name --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 1000
#Set unique names for each sample by adding the dataset name.
    	awk -v dataset="$name" '{print dataset"@"$2"\t"dataset"@"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $name.fam > tmp.fam
    	mv tmp.fam $name.fam
done
 
# III. For the vcf files from folders “SubsetHailJointCall”. These files didn’t have required FILTER tags in the header. Also, the GT field was named as LGT, which prevented accurate processing of these files with bcftools and plink.
 
ml bcftools
#Did the following steps for vcfs under “SubsetHailJointCall” folders to correct headers
bcftools view -h /sc/arion/projects/ibdgc/uchicago-transfer/data/wes/updates/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes_SubsetHailJointCall/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes.vcf.bgz > hdr.txt
##add necessary filter flags in the hdr.txt:
##FILTER=<ID=VQSRTrancheSNP99.00to99.90+,Description="VQSRTrancheSNP99.00to99.90+">
##FILTER=<ID=VQSRTrancheINDEL99.00to99.90+,Description="VQSRTrancheINDEL99.00to99.90+">
##FILTER=<ID=VQSRTrancheINDEL99.00to99.90,Description="VQSRTrancheINDEL99.00to99.90">
##FILTER=<ID=VQSRTrancheSNP99.00to99.90,Description="VQSRTrancheSNP99.00to99.90">
 
 
bcftools reheader -h hdr.txt -o anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes.vcf.bgz /sc/arion/projects/ibdgc/uchicago-transfer/data/wes/updates/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes_SubsetHailJointCall/anvil_ccdg_broad_ai_ibd_daly_newberry_share_wes.vcf.bgz
 
#Following steps for the header corrected files
ml plink2/v2.00a3.3
ml plink
ml bcftools
file_list=$(ls | grep ^anvil)
 
for file in $file_list
do
#Selected variant that passed the QC check that was already done. Correct GT field tag with sed and normalized and split multi-allelic variants.
    	bcftools view -O v -f PASS $file | sed -e 's/LGT/GT/g' | bcftools norm -c s -m -any -Oz -o PASS_norm_$file --fasta-ref ../ Homo_sapiens.GRCh38.dna.toplevel.fa
    	name="${file%???????}"
#Converted into plink format and set variant ids.
    	plink2 --vcf PASS_norm_$file --make-bed --out $name --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 1000
#Set unique names for each sample by adding the dataset name.
    	awk -v dataset="$name" '{print dataset"@"$2"\t"dataset"@"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $name.fam > tmp.fam
    	mv tmp.fam $name.fam
done
 
#Selected overlapping variant IDs and generate a file containing those variant IDs named WESIDs.txt.
#Extracted overlapping variants from each site-specific plink file.
file_list=$(ls | grep bed$)
for i in $file_list
do
    	plink2 --bfile ${i%????} --make-bed --extract WESIDs.txt --out selected_${i%????}
    	echo selected_${i%????} >> WES_merge_list.txt
done
#Merged all site-specific plink files together (only available in plink v.1.9).
plink --merge-list WES_merge_list.txt --out WES_merged --keep-allele-order
 
#6. Selected overlapping variants from all merged plink files and further merged them together
 
awk -F'\t' 'NR==FNR{c[$2]++;next};c[$2] > 0' gsa_allCHRmerged_imputed.bim nonEUR_GSA_merged.bim > tmp1
awk -F'\t' 'NR==FNR{c[$2]++;next};c[$2] > 0' tmp1 mergedOmni.bim > tmp2
awk -F'\t' 'NR==FNR{c[$2]++;next};c[$2] > 0' tmp2 mergedAxiom.bim > tmp3
awk -F'\t' 'NR==FNR{c[$2]++;next};c[$2] > 0' tmp3 cutler-AA-WGS-big_daly.bim > tmp4
awk -F'\t' 'NR==FNR{c[$2]++;next};c[$2] > 0' tmp4 WES_merged > allIDs.txt
rm tmp*
 
file_list=$(ls | grep bed$)
for i in $file_list
do
    	plink2 --bfile ${i%????} --make-bed --extract allIDs.txt --out selected_${i%????}
    	echo selected_${i%????} >> final_merge_list.txt
done
 
plink --merge-list final_merge_list.txt --out merged_all --keep-allele-order

#7. Ran king to check if the same sample IDs from the same sites were indeed the same samples (yes, they were!)

king --b WES_merged.bed --duplicate --prefix king_for_wes

#removed recurring samples from the same sites

plink2 --bfile merged_all --remove same_samples_from_the_samesites.txt --make-bed --out final_merged_for_king
 
#8. Ran KING for all data
 
ml king/2.2.5
 
king -b final_merged_for_king.bed --duplicate --prefix final_king
