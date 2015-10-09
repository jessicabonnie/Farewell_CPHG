#!/bin/bash
#This script details how to create SNP position update file for the immunochip


#File downloaded on 4/22 from: ftp://ussd-ftp.illumina.com/downloads/ProductFiles/HumanImmuno/Immuno_BeadChip_11419691_B.csv
manifest=Immuno_BeadChip_11419691_B.csv

#This file is downloaded from T1Dbase http://www.t1dbase.org/downloads/data/immunochip/hg19_immunochip_all_diseases.gff3 on 7/21/2015
all_diseases=hg19_immunochip_all_diseases.gff3.txt


#The Quinlan Bible contains all the information we need for the regular SNPs
qb=quinlan-immunochip-snps-annotated-2013-Apr-24b.txt

sed 's/\r//g' ${qb} | awk 'NR>1{print $7,$3}' > mapupdate_Apr24b.tmp

sed 's/\r//g' ${qb} | awk 'NR>1{print $7,$1}' > mapupdate_Apr24b_chr.tmp

#For the INDELs, however, we need to use a combination of the manifest and the t1dbase reference file


grep "INDEL" ${all_diseases} > grepforINDEL.tmp

sed -i 's/;/\t/g; s/ID=//g; s/\r//g' grepforINDEL.tmp

#There is an error in the gff3 file for one of the long Illumina SNP IDs:
sed -i 's/seq-VH-94Z-1_M_R_177145924/seq-VH-94Z-1_M_R_1771459240/g' grepforINDEL.tmp

#We need to merge with the manifest in order to retrieve the usual SNP IDs to match with the positions
awk '{print $9, $5}'  grepforINDEL.tmp > t1dbase_indels.tmp

awk 'BEGIN {FS = ","}; NR>8 && NR<196534 {print $1,$2}' ${manifest} | sed 's/\r//g' > manifest_name_key.txt

join <(cat manifest_name_key.txt| sort) <( cat t1dbase_indels.tmp | sort) | awk '{print $2,$3}' > indels_pos_update.tmp

cat mapupdate_Apr24b.tmp > map_pos_update.txt
cat indels_pos_update.tmp >> map_pos_update.txt
