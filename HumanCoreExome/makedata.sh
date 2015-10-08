#!/bin/bash
#!/usr/bin/bash


echo -e "#########################################################################"
echo -e
echo -e "# GWAS-Pipedream - Standard GWAS QC Pipeline Package"
echo -e "# (c) 2014-2015 JBonnie, WMChen"
echo -e
echo -e "#########################################################################"

### This script transforms output from Genome Studio (Full Data Tables in given format -- e.g. Top Alleles, GType) into PLINK format.


##############################################################################
# Usage: sh  make_data.sh {DATAFILE} {nickname} {format}
##############################################################################

##################################################################################################################
#Arguments
# DATAFILE : the path to the full data table exported from Genome Studio.
#		It is expected that the table will contain more than one column per sample.
# nickname : string, the alias for the project, used in all filenames
# format : the column type of interest e.g. "Top Allele" or "Gtype"
#
########################################################################################################



# Definining/Reading user specified script parameters/variables/values
DATAFILE=$1
OUTNAME=$2
format=$3



echo -e "########################################################################################################################################"
echo -e
echo -e "# Title: MAKE DATA - Standard GWAS QC Pipeline"
echo -e "#"
 echo -e "#"
echo -e "# Disc: Transforms output from Genome Studio into plink files."
echo -e "#"
echo -e "# Note: PLINK output files will NOT contain sex, status, or family information, those are added during pheno_inc.sh"
echo -e "# Note: {format} must match the COLUMN NAME of the desired genotyping format."
echo -e "# Note: If {format} contains a space QUOTES must be used."
echo -e "# Note: DATAFILE expected to have more than one column per sample."
echo -e "# Note: The first two genotype columns are expected to occur within the first 20 columns of DATAFILE."
echo -e "#"
echo -e "# Usage: sh  ~/SHcode/make_data.sh {DATAFILE} {OUTNAME} {format}"
echo -e "#"
echo -e "# See Script for option/parameter details"
echo -e "#"
echo -e "# by:  jbonnie"
echo -e "# date: 10.25.13"
echo -e
echo -e "Full Data Table: ${DATAFILE}"
echo -e "Name/Path for plink output: ${OUTNAME}"
echo -e "Allele Format (column name): ${format}"
echo -e "Performed on:"
date
echo -e
echo -e "########################################################################################################################################"




# project_folder=$(pwd)


#It's easier to make variables to indicate these subfolders and files

raw_folder=$(dirname $OUTNAME)/makedata_tmp

echo -e
echo -e "\n------------------"
echo "Creating Raw Data Folder"
echo -e "------------------\n"
echo -e
mkdir ${raw_folder}
cd ${raw_folder}

echo -e
echo -e "\n------------------"
echo "Storing column numbers from the Full Data Table &"
echo "Determining Uniform Number of Columns for Each Sample"
echo -e "------------------\n"
echo -e


#Lets read some columns numbers!


#Locate the first two occurances of the pattern ${format}", and then find the difference.
#This tells us the uniform number of columns between each genotype column!

gencols=$(cut -f1-20 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk -v pattern="${format}" 'BEGIN {FS="\t"}; $0 ~ pattern {print NR}')
gen1=$( echo ${gencols} | awk '{print $1}')
gen2=$( echo ${gencols} | awk '{print $2}')
dif=$((${gen2} - ${gen1}))

#Locate the first occurance of the pattern ".Top Allele", that is the first data column of interest.
#top1=$(cut -f1-8 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/.Top Allele/ {print NR}')

#The second occurance of the pattern ".Top Allele" cannot be more than 10 columns in, right? Lets Find it!
#Once we know where it is we will know the uniform number of columns between each ".Top Allele"!
#dif=$(cut -f$((${top1}+1))-10 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/.Top Allele/ {print NR}')
#top2=$((${top1} + ${dif}))

#We assume that the preliminary columns also occur within the first 20 columns

chr=$(cut -f1-20 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/^Chr$/ {print NR}')
pos=$(cut -f1-20 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/^Position$/ {print NR}')
snp=$(cut -f1-20 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/^Name$/ {print NR}')


echo -e
echo -e "Column Numbers"
echo -e "Chromosome : ${chr}"
echo -e "SNP : ${snp}"
echo -e "SNP POSITION : ${pos}"
echo -e "First ${format} : ${gen1}"
echo -e "Second ${format} : ${gen2}"
echo -e "Difference (Number of Columns per Sample) : ${dif}"
echo -e




echo -e
echo -e "\n------------------"
echo "Creating PLINK TPED: Transferring Allele Data"
echo -e "------------------\n"
echo -e


## For each genotype column, print the allele information, replacing "--" with "0 0" and adding spaces between the alleles
## Also, for some reason, flip the alleles, this was how it was done in inherited skeleton, so it has been left that way.
formattitle=$(echo ${format} | sed 's/ //g')

awk 'NR>1' $DATAFILE | sed "s/\r//g" |\
  awk -v gen1=${gen1} -v gen2=${gen2} -v dif=${dif} '{printf("%s", $gen1); for(i=gen2;i<=NF;i+=dif) printf(" %s", $i);printf("\n");}' |\
  sed 's/--/0 0/g' |\
sed 's/\([A-Z]\)\([A-Z]\)/\2 \1/g' > ${raw_folder}/${formattitle}Allele_tped.tmp


echo -e
echo -e "\n------------------"
echo "Creating PLINK TPED: Transferring SNP Information"
echo -e "------------------\n"
echo -e


# The first four columns of the tped hold the SNP info from the table with "0" place holders in the "Genetic Distance" column (column 3)

awk 'NR>1' $DATAFILE | sed "s/\r//g" | awk -v chr=${chr} -v pos=${pos} -v snp=${snp} '{printf("%s %s 0 %d\n", $chr,$snp,$pos);}' > ${raw_folder}/snpInfo_tped.tmp

echo -e
echo -e "\n------------------"
echo "Creating PLINK TPED: Pasting TMP files"
echo -e "------------------\n"
echo -e

# create tped by pasting snpInfo_tped.tmp and ${formattitle}Allele_tped.tmp
paste -d" " ${raw_folder}/snpInfo_tped.tmp ${raw_folder}/${formattitle}Allele_tped.tmp > ${raw_folder}/data.tped



echo -e
echo -e "\n------------------"
echo "Creating PLINK TFAM"
echo -e "------------------\n"
echo -e

# -9 means unknown/missing affection status
# take header, replace tabs with newlines, and strip windows newlines; take only the lines containing 'GType' and then strip '.GType' to leave only sampleIDs;
# NOT SURE WHAT THIS IF STATEMENT WOULD CATCH.... maybe a colname that is "FID IID.GType"... 
# regardless, outside of the if statement the IID is doubled to also be FID and parents and affection status are defaulted to their missing values

head -1 $DATAFILE | sed "s/\r//g" | sed "s/\t/\n/g" |\
 grep "${format}" | sed "s/.$format//g" |\
 awk '{if(NF==2){print $1,$2,0,0,0,-9}else{print $1,$1,0,0,0,-9}}' > ${raw_folder}/data.tfam


#head -1 $DATAFILE | sed "s/\r//g" | sed "s/\t/\n/g" |\
#    grep "Top Alleles" | sed "s/.Top Alleles//g" |\
#    awk '{if(NF==2){print $1,$2,0,0,0,-9}else{print $1,$1,0,0,0,-9}}' > data.tfam




echo -e
echo -e "\n------------------"
echo "Creating PLINK BINARY FILES"
echo -e "------------------\n"
echo -e


plink --tfile ${raw_folder}/data --make-bed --out ${OUTNAME} --noweb --allow-no-sex


echo -e
echo -e "\n------------------"
echo "Delete Extraneous Files and Tar/Zip the Rest"
echo -e "------------------\n"
echo -e

#rmlist=$(ls data.* *.tmp *.hh *.nosex)
#rm -f $rmlist

cd ${raw_folder}/..
tar -czf ${raw_folder}.tar.gz ${raw_folder}
rm -rf ${raw_folder}
rm -f ${OUTNAME}.nof ${OUTNAME}.hh ${OUTNAME}.nosex
#tar --create -f makedata_tmp.tar ${raw_folder}

#gzip raw.tar


echo 'makedata.sh completed: '$(date)
