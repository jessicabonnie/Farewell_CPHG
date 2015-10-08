#!/bin/bash
#!/usr/bin/bash

### This script runs initial relatedness checking on data after the initial QC step.
### It uses files in the 2_QC1 folder as well as the covariable lists created during pheno_inc.sh



##############################################################################
# Usage: sh  SHcode/relatedness_qc.sh ${nickname} ${covariatevalue} ${chip}
##############################################################################

####################################################################3
#Arguments
#nickname : string, the alias for the project, used in all filenames
#covariablevalue : integer, the number of WHICH covariable should be included in the table for checking (USUALLY this is the cohort, #1)
#chip : I/E, character value indicating whether data was generated from the Immunochip (I) or the HumanCoreExomeChip (E)
#
#####################################################################



echo -e "#########################################################################"
echo -e
echo -e "# GWAS-Pipedream - Standard GWAS QC Pipeline Package"
echo -e "# (c) 2014-2015 JBonnie, WMChen"
echo -e
echo -e "#########################################################################"

nickname=$1
covariablevalue=$2
chip=$3

#if [ ${covariablevalue} == '' ]; then covariablevalue=1 fi

covcol=$((${covariablevalue} + 2))
project_folder=$(pwd)
qc1_folder=${project_folder}/QC1
#qc3_folder=${project_folder}/kinship
nb_folder=${project_folder}/NB
#mkdir ${qc3_folder}
mkdir ${nb_folder}


raw_folder=${project_folder}/raw

cov=${qc1_folder}/${nickname}.cov


cd ${qc1_folder}
echo "moved"
#Create necessary tables for drawing and further relatedness QC if this is the first covariable
if [ ${covariablevalue} -eq 1 ]; then
king -b ${nickname}6.bed --related --degree 5 --errorrate 0.003 --prefix ${nickname}6
# > ${nickname}6kin.out
fi

#Add designated covariable col to both kin and kin0
field_count=$(awk 'NR==1{print NF}' ${nickname}6.kin0)
covariablecol=$((${field_count} + ${covcol}))
echo $(head -n1 ${nickname}6.kin0) "Covariable1 Covariable2"  > ${nickname}6kin0_cov${covariablevalue}.txt
awk -v cov=${cov} 'BEGIN{while((getline<cov)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${nickname}6.kin0 | awk -v covcol=${covariablecol} '{print $1,$2,$3,$4,$5,$6,$7,$8,$covcol}' > ${nickname}6kin0.tmp
awk -v cov=${cov} 'BEGIN{while((getline<cov)>0)l[$2]=$0}$4 in l{print $0"\t"l[$4]}' ${nickname}6kin0.tmp | awk -v covcol=${covariablecol} '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$(covcol+1)}' >> ${nickname}6kin0_cov${covariablevalue}.txt

#Note, we are changing the kin file slightly to make it match the kin0
field_countkin=$(awk 'NR==1{print NF}' ${nickname}6.kin)
covariablecolkin=$((${field_countkin} + ${covcol}))
echo $(head -n1 ${nickname}6.kin0) "Covariable1 Covariable2"  > ${nickname}6kin.txt
awk -v cov=${cov} 'BEGIN{while((getline<cov)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${nickname}6.kin | awk -v covcol=${covariablecolkin} '{print $1,$2,$1,$3,$4,$7,$8,$9,$covcol}' > ${nickname}6kin.tmp
awk -v cov=${cov} 'BEGIN{while((getline<cov)>0)l[$2]=$0}$4 in l{print $0"\t"l[$4]}' ${nickname}6kin.tmp | awk -v covcol=${covariablecol} '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$(covcol+1)}' >> ${nickname}6kin.txt



#Determine relatedness threshold
#Threshold is lower for the Immunochip because there is less coverage
usual_related_threshold=0.0884
imchip_related_threshold=0.1768
icaa_related_threshold=0.1
threshold=0.18
firstdegree=0.1


#Create relatedness files for checking
awk -v relat=${threshold} '$9!=$10 && $8>relat' ${nickname}6kin0_cov${covariablevalue}.txt > relat_btwn_covariables_cov${covariablevalue}.nb
awk -v relat=${threshold} '$8>relat' ${nickname}6kin0_cov${covariablevalue}.txt > related_individuals_cov${covariablevalue}.nb

#cp related_individuals_cov${covariablevalue}.nb relat_btwn_covariables_cov${covariablevalue}.nb ${nb_folder}


cd ${qc1_folder}
#Duplicate Identification
####MIGHT NEED TO MAKE A DECISION ABOUT MULTIDUPLICATES HERE####  
### ALSO, turns out that when we have related samples that are duplicates, they get left out if we don't also look in the kin file!!
awk 'NR>1 && $NF > 0.4' ${qc1_folder}/${nickname}6.kin0 | awk '{print $1, $2}' > dup1.txt
awk 'NR>1 && $9 > 0.4' ${qc1_folder}/${nickname}6.kin | awk '{print $1, $2}' >> dup1.txt
awk 'NR>1 && $8 > 0.4' ${qc1_folder}/${nickname}6kin0_cov${covariablevalue}.txt  > dupcovariables_cov${covariablevalue}.txt
awk 'NR>1 && $8 > 0.4' ${qc1_folder}/${nickname}6kin.txt  >> dupcovariables_cov${covariablevalue}.txt
awk 'NR>1 && $9!=$10 && $8 > 0.4' ${qc1_folder}/${nickname}6kin0_cov${covariablevalue}.txt  > dup_btwn_covariables_cov${covariablevalue}.txt
awk 'NR>1 && $9!=$10 && $8 > 0.4' ${qc1_folder}/${nickname}6kin.txt  >> dup_btwn_covariables_cov${covariablevalue}.txt
awk 'NR>1 && $9==$10 && $8 > 0.4' ${qc1_folder}/${nickname}6kin0_cov${covariablevalue}.txt  > dup_within_covariables_cov${covariablevalue}.txt
awk 'NR>1 && $9==$10 && $8 > 0.4' ${qc1_folder}/${nickname}6kin.txt  >> dup_within_covariables_cov${covariablevalue}.txt


#Here the relationships between related samples are defined
head -n1  ${qc1_folder}/${nickname}6kin0_cov${covariablevalue}.txt | awk '{print $1, $2, $3, $4, $7, $8, $9, $10,"Relationship"}'> relationship_cov${covariablevalue}.nb
awk 'NR > 1 && $8 > 0.4' ${qc1_folder}/${nickname}6kin0_cov${covariablevalue}.txt | awk '{print $1, $2, $3, $4, $7, $8, $9, $10, "Duplicate"}' >> relationship_cov${covariablevalue}.nb
awk 'NR > 1 && $8 > 0.4' ${qc1_folder}/${nickname}6kin.txt | awk '{print $1, $2, $3, $4, $7, $8, $9, $10, "Duplicate"}' >> relationship_cov${covariablevalue}.nb
awk -v FS=${firstdegree} 'NR > 1 && $8 < 0.4 && $8 > FS && $7 < 0.005' ${qc1_folder}/${nickname}6kin0_cov${covariablevalue}.txt | awk '{print $1, $2, $3, $4, $7, $8, $9, $10, "PO"}' >> relationship_cov${covariablevalue}.nb
awk -v FS=${firstdegree} 'NR > 1 && $8 < 0.4 && $8 > FS && $7 > 0.005' ${qc1_folder}/${nickname}6kin0_cov${covariablevalue}.txt | awk '{print $1, $2, $3, $4, $7, $8, $9, $10, "FS"}' >> relationship_cov${covariablevalue}.nb

awk -v FS=${firstdegree} -v SD=${threshold} 'NR > 1 && $8 <= FS && $8 > SD' ${qc1_folder}/${nickname}6kin0_cov${covariablevalue}.txt | awk '{print $1, $2, $3, $4, $7, $8, $9, $10, "2nd"}' >> relationship_cov${covariablevalue}.nb


yes | cp *.nb ${nb_folder}
