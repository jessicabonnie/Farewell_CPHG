#!/bin/bash
#!/usr/bin/bash

### This script runs basic sample and SNP QC on PLINK files.
### It produces tables containing samples for follow up with the genotyping lab.


##############################################################################
# Usage: sh  qc1.sh {nickname} {sample_remove}
##############################################################################

####################################################################3
#Arguments
#
#nickname : string, the alias for the project, used in all filenames
#sample_remove : T/F, Do we want to exclude a number of samples right at the beginning. This is normally used when qc is being run for a second time.
#		If so, it had better be contained in ${project_folder}/samplesremove_at_start.txt.
#
#####################################################################


# Definining/Reading user specified script parameters/variables/values
plinkraw=$1
nickname=$2
sample_remove=$3
project_folder=$(pwd)


echo -e "#########################################################################"
echo -e
echo -e "# GWAS-Pipedream - Standard GWAS QC Pipeline Package"
echo -e "# (c) 2014-2015 JBonnie, WMChen"
echo -e
echo -e "#########################################################################"


echo -e "########################################################################################################################################"
echo -e
echo -e "# Title: FIRST QC SCRIPT - Standard GWAS QC Pipeline"
echo -e "#"
 echo -e "#"
echo -e "# Disc: Run basic sample and SNP QC on PLINK files. Produces sample lists for follow up with genotyping lab."
echo -e "#"
echo -e "# Note: If duplicates samples are included in the samplesremove_at_start file, it will be removed TWICE."
echo -e "# Note: Data should be in binary plink format."
echo -e "#"
echo -e "# Usage: sh  ./qc1.sh {plinkraw} {nickname} {sample_remove}"
echo -e "#"
echo -e "# See Script for option/parameter details"
echo -e "#"
echo -e "# by:  jbonnie, WMChen"
echo -e "# date: 08.25.14"
echo -e
echo -e "Raw Plink File input: ${plinkraw}"
echo -e "Study Alias: ${nickname}"
echo -e "Remove Samples from Data at Start: ${sample_remove}"
echo -e "Performed on:"
date
echo -e
echo -e "########################################################################################################################################"




#It's easier to make variables to indicate these subfolders and files, also make them if they don't exist
qc1_folder=${project_folder}/QC1
nb=${project_folder}/NB


mkdir ${qc1_folder}
mkdir ${nb}

cd ${qc1_folder}

cp ${plinkraw}.bim ${nickname}1.bim
cp ${plinkraw}.bed ${nickname}1.bed
cp ${plinkraw}.fam ${nickname}1.fam
cp ${plinkraw}.log ${nickname}1.log

#add headers to lists of snps and samples to be removed, which will be appended to throughout the program
echo "FID IID REASON" > sampletoberemoved.txt
echo "SNP REASON" > snptoberemoved.txt


echo -e
echo -e "\n------------------"
echo "Adding Samples to the SampleRemove list before any analysis: ${sample_remove}"
echo -e "------------------\n"
echo -e

#if fed a sample remove list from the start -- note, if a duplicate is being removed, adding it to the sampletoberemoved at this point will lead to both being removed in the next two steps:
if [ "${sample_remove}" = "T" ]; then
  awk '{print $1, $2, $3}' ${project_folder}/samplesremove_at_start.txt >> sampletoberemoved.txt
  #cp ${nickname}1.log ${nickname}1orig.log
  #cp ${nickname}1.bim ${nickname}1orig.bim
  #cp ${nickname}1.bed ${nickname}1orig.bed
  #cp ${nickname}1.fam ${nickname}1orig.fam
  #plink --noweb --make-bed --bfile ${nickname}1orig --out ${nickname}1 --remove ${project_folder}/samplesremove_at_start.txt


fi

echo -e
echo -e "\n------------------"
echo "STEP 1: KING by SNP"
echo -e "------------------\n"
echo -e

#use king to check snps
king -b ${nickname}1.bed --bySNP --prefix ${nickname}1

#identify and remove non-Ychr missing SNPs and first set of monomorphic snps
awk '$11 < 0.8 && $2 != "Y"' ${nickname}1bySNP.txt | awk '{print $1}' > allmissingSNP.txt
awk '$8+$9==0 || $9+$10==0' ${nickname}1bySNP.txt | awk '{print $1}' > monomorphicSNP.txt

awk '{print $1, "CallRateLessThan80"}' allmissingSNP.txt >> snptoberemoved.txt
awk 'NR>1{print $1, "Monomorphic"}' monomorphicSNP.txt >> snptoberemoved.txt

plink --bfile ${nickname}1 --exclude snptoberemoved.txt --make-bed --out ${nickname}2 --noweb --remove sampletoberemoved.txt

#identify indels
awk '$5=="I" || $5=="D"' ${nickname}2.bim | awk '{print $2}' > indels.txt

king -b ${nickname}2.bed --bysample --prefix ${nickname}2

#samples which failed more than %5 of the time added to lessthan95.txt and then removed
awk 'NR > 1 && $7 > 0.05' ${nickname}2bySample.txt | awk '{print $1,$2}' > lessthan95.txt

plink --bfile ${nickname}2 --remove lessthan95.txt --make-bed --out ${nickname}3 --noweb

awk '{print $1,$2, "MissingMoreThan5"}' lessthan95.txt >> sampletoberemoved.txt


king -b ${nickname}3.bed --bySNP --prefix ${nickname}3

# check for low call rate snps and poor Y snps

#Determine the relative ratio of males to males+females in order to check for poor Y-SNPs
males=$(awk '$5==1' ${nickname}3.fam | wc -l)
females=$(awk '$5==2' ${nickname}3.fam | wc -l)
mratio=$(echo ${males}/$((${males}+${females}))|bc -l)
ratio=$(echo ${mratio}*1.1|bc -l)
echo "The ratio of males:(males+females) is " ${mratio} ", which means we will cut Y SNPs with call rates greater than " ${ratio}
awk -v r=${ratio} '$2=="Y" && $11 > r' ${nickname}3bySNP.txt | awk '{print $1}' > poorYSNP.txt

#awk '$2=="Y" && $11 > 0.6' ${nickname}3bySNP.txt | awk '{print $1}' > poorYSNP.txt



# check for low call rate snps and poor Y snps
awk '$11 < 0.95 && $2!="Y"' ${nickname}3bySNP.txt | awk '{print $1}' > lowcallrateSNP.txt

cat lowcallrateSNP.txt > lowcall_poorY.txt
cat poorYSNP.txt >> lowcall_poorY.txt

plink --bfile ${nickname}3 --exclude lowcall_poorY.txt --make-bed --out ${nickname}4a --noweb


#Additional Y-chromosome SNP QC
#Note: This is being done before removing "misreported" sex samples, so some "Good" Y SNPs might be removed.
#Since we don't give a fweep about the Y chromosome for our analyses in complex diseases and are only doing these Y SNP QC steps
#in order to improve our gender QC, we are leaving it here. If someone else is using this, you could move the second
#(filtered for females only) step to a point after the mislabeled gender and also fiddle with the fudge factor (now 1.1)
#to make it less strict.

#filter females and then use call rates on that set to determine second set of poor Y SNPs
plink --bfile ${nickname}4a --filter-females --out ${nickname}4afemale --noweb --make-bed
king -b ${nickname}4afemale.bed --bySNP --prefix ${nickname}4afemale
awk '$2=="Y" && $11 > 0.04' ${nickname}4afemalebySNP.txt | awk '{print $1}' > poorYSNP2.txt

plink --bfile ${nickname}4a --exclude poorYSNP2.txt --make-bed --out ${nickname}4 --noweb

#Add all low call rate snps and poor y snps to master list
awk '{print $1, "CallRateLessThan95"}' lowcallrateSNP.txt >> snptoberemoved.txt
awk '{print $1, "YInFemales"}' poorYSNP.txt >> snptoberemoved.txt
awk '{print $1, "YInFemales"}' poorYSNP2.txt >> snptoberemoved.txt

king -b ${nickname}4.bed --bysample --bySNP --prefix ${nickname}4 > ${nickname}4bySample.log


#Determine thresholds for use in identifying erroneaously male and female samples, originally this was hard coded to 700 ysnps. Also thresholds determinded for midsex and X0 errors
ysnps=$(grep "Y-chromosome SNPs" ${nickname}4bySample.log | awk -F, '{print $3}'| awk -F' ' '{print $1}')

half_ysnps=$(echo ${ysnps} / 2 | bc )
third_ysnps=$(echo ${ysnps} / 3 | bc )
twothird_ysnps=$(echo "2 * ${third_ysnps}" | bc)
sixth_ysnps=$(echo ${ysnps} / 6 | bc)
fivesixth_ysnps=$(echo ${ysnps} - ${sixth_ysnps} | bc)

max_xhetero=$(awk 'NR == 2 {max=$10 ; min=$10} $10 >= max {max = $10} $10 <= min {min = $10} END { print max }' ${nickname}4bySample.txt)
half_max_xhetero=$(echo ${max_xhetero} / 2 | bc -l )

awk -v hy=${half_ysnps} 'NR>1 && $5==2 && $11 > hy' ${nickname}4bySample.txt | awk '{print $1,$2}' > femaleerror.txt

awk -v hy=${half_ysnps} 'NR>1 && $5==1 && $11 < hy' ${nickname}4bySample.txt | awk '{print $1,$2}' > maleerror.txt
cat femaleerror.txt > sexerror.txt
cat maleerror.txt >> sexerror.txt

##lists specifically to be sent to GSL to look at based on gender graph
## experimentally edited to use the ysnp count for boundaries as well as half the max x-heterozygosity
#awk 'NR>1 && $11 > 500 && $11 < 1100' ${nickname}4bySample.txt | awk '{print $1,$2}' > midsexerror.nb
#awk 'NR>1 && $10 < 0.05 && $11 < 200' ${nickname}4bySample.txt | awk '{print $1,$2}' > X0error.nb
awk -v twothirdsy=${twothird_ysnps} -v onethirdy=${third_ysnps} 'NR>1 && $11 > onethirdy && $11 < twothirdsy' ${nickname}4bySample.txt | awk '{print $1,$2}' > midsexerror.nb
head -n1 ${nickname}4bySample.txt > midsexerror_suna.nb
cat midsexerror.nb >> sexerror.txt
awk -v twothirdsy=${twothird_ysnps} -v onethirdy=${third_ysnps} 'NR>1 && $11 > onethirdy && $11 < twothirdsy' ${nickname}4bySample.txt >> midsexerror_suna.nb
#awk -v sixthy=${sixth_ysnps} -v halfhetero=${half_max_xhetero} 'NR>1 && $10 < halfhetero && $11 < sixthy' ${nickname}4bySample.txt | awk '{print $1,$2}' > X0error.nb
awk -v sixthy=${sixth_ysnps}  'NR>1 && $10 < 0.05 && $11 < sixthy' ${nickname}4bySample.txt | awk '{print $1,$2}' > X0error.nb
cat X0error.nb >> sexerror.txt
head -n1  ${nickname}4bySample.txt > X0error_suna.nb
awk -v sixthy=${sixth_ysnps} 'NR>1 && $10 < 0.05 && $11 < sixthy' ${nickname}4bySample.txt >> X0error_suna.nb

awk -v fivesixthy=${fivesixth_ysnps} 'NR>1 && $10 > 0.05 && $11 > fivesixthy' ${nickname}4bySample.txt | awk '{print $1,$2}' > XXYerror.nb
cat XXYerror.nb >> sexerror.txt
head -n1 ${nickname}4bySample.txt > XXYerror_suna.nb
awk -v fivesixthy=${fivesixth_ysnps} 'NR>1 && $10 > 0.05 && $11 > fivesixthy' ${nickname}4bySample.txt >> XXYerror_suna.nb



##The the half the number of Y snps in the dataset is used as a threshold to define missing gender. Originally 700 was used.
awk '$5==0' ${nickname}4bySample.txt | awk -v hy=${half_ysnps} '{if($11<hy){print $1, $2, 2}else{print $1,$2,1}}' > updatesex2.txt

plink --bfile ${nickname}4 --remove sexerror.txt --update-sex updatesex2.txt --make-bed --out ${nickname}5 --noweb

#Create lists for checking
awk '{print $1, $2, "MislabeledAsFemale"}' femaleerror.txt >> sampletoberemoved.txt
awk '{print $1, $2, "MislabeledAsMale"}' maleerror.txt >> sampletoberemoved.txt
awk '{print $1, $2, "X0error"}' X0error.nb >> sampletoberemoved.txt
awk '{print $1, $2, "XXYerror"}' XXYerror.nb >> sampletoberemoved.txt
awk '{print $1, $2, "MidSexError"}' midsexerror.nb >> sampletoberemoved.txt

king -b ${nickname}5.bed --bySNP --prefix ${nickname}5


#Any Y SNPs with a number of heterozygous samples greater than one tenth the number of males, should be removed
males=$(awk '$5==1' ${nickname}5.fam | wc -l)
males10=$(echo "$males *.1" | bc -l)
echo "Any Y SNPs with a number of heterozygous samples greater than one tenth the number of males, should be removed. We have" ${males} "males"
awk -v males10=${males10} '$2=="Y" && $9>males10' ${nickname}5bySNP.txt | awk '{print $1}' > YwHeterozygosity10.txt

plink --bfile ${nickname}5 --exclude YwHeterozygosity10.txt --noweb --make-bed --out ${nickname}6a
awk '{print $1, "YwHeterozygosity10"}' YwHeterozygosity10.txt >> snptoberemoved.txt

#Check again for monomorphic SNPs now that samples have been removed and filter them out
king -b ${nickname}6a.bed --bySNP --prefix ${nickname}6a
awk '$8+$9==0 || $9+$10==0' ${nickname}6abySNP.txt | awk '{print $1}' > monomorphicSNP2.txt
plink --bfile ${nickname}6a --exclude monomorphicSNP2.txt --noweb --make-bed --out ${nickname}6
awk 'NR>1{print $1, "Monomorphic2"}' monomorphicSNP2.txt >> snptoberemoved.txt



echo "Completed qc steps in qc.sh, now moving select files to NB folder, " ${nb} "for checking: " $(date)
yes | cp YwHeterozygosity10.txt *.nb sexerror.txt sampletoberemoved.* snptoberemoved.* ${nb}

#echo "Now moving select files to storage.tar: " $(date)
#rmlist=$(ls ${nickname}2* ${nickname}3* *.tmp *.hh *.nosex *TMP* ${nickname}4.*)
#tar --create -f storage.tar ${rmlist}

#echo "Now zipping storage.tar: " $(date)
#gzip storage.tar
#rm -f ${rmlist}


