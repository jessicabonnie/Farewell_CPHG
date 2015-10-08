#!/bin/bash


script_folder=/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Farewell/IMCHIP

echo -e "#########################################################################"
echo -e
echo -e "# GWAS-Pipedream - Standard GWAS QC Pipeline Package"
echo -e "# (c) 2014-2015 JBonnie, WMChen"
echo -e
echo -e "#########################################################################"

#cd data
project_folder=$(pwd)
overall_title=$1
nickname=$2
chip=$3
covariablevalue=$4
covcol=$((${covariablevalue} + 2))

covariablefile=${project_folder}/covariables${covariablevalue}.list

colorfile=${project_folder}/colors.list

qc1_folder=${project_folder}/QC1
cov=${qc1_folder}/${nickname}.cov
graph_folder=${project_folder}/graphs_cov${covariablevalue}
nb_folder=${project_folder}/NB


mkdir ${graph_folder}

#cd ${qc1_folder}
cd ${graph_folder}
echo "Now beginning qc_pdf.sh: " $(date)
gender_table=${qc1_folder}/${nickname}4bySample.txt
gender_covariable_table=${qc1_folder}/${nickname}4bySample_covariables${covariablevalue}.txt
echo ${covariablefile}
field_count=$(awk 'NR==1{print NF}' ${gender_table})
covariablecol=$(( ${field_count} + ${covcol} ))
#echo $(head -n1 ${gender_table} ) "Covariable"  > ${gender_covariable_table}
head -n1 ${gender_table} | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,"Covariable"}'> ${gender_covariable_table}
awk -v cov=${cov} 'BEGIN{while((getline<cov)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${gender_table} | awk -v covcol=${covariablecol} 'NR>1{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$covcol}' >> ${gender_covariable_table}

echo $cov

#awk -v cov=${cov} 'BEGIN{while((getline<cov)>0)l[$1]=$0}$1 in l{print $0"\t"l[$1]}' ${gender_table}| awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$NF}' > ${gender_covariable_table}

new_gender_table=${qc1_folder}/${nickname}4bySample_R.txt

cp ${gender_covariable_table} ${new_gender_table}

sed -i 's/#/_/g' ${new_gender_table}

gender_table=${new_gender_table}

gender_graph=${graph_folder}/${nickname}gender.ps
#draw gender pdf
echo ${covariablefile}
R CMD BATCH "--args ${gender_table} ${gender_graph} ${overall_title} ${nickname} ${covariablefile} ${colorfile}" ${scripts}/gender.R

rawrel_table=${qc1_folder}/${nickname}6kin0_cov${covariablevalue}.txt
new_rawrel_table=${qc1_folder}/${nickname}6kin0_cov${covariablevalue}_R.txt
cp ${rawrel_table} ${new_rawrel_table}
sed -i 's/#/_/g' ${new_rawrel_table}
rawrel_table=${new_rawrel_table}
echo $rawrel_table

btwn_covariables_table=${qc1_folder}/relat_btwn_covariables_cov${covariablevalue}.nb
Rbtwn_covariables_table=${qc1_folder}/relat_btwn_covariables_cov${covariablevalue}_R.nb
cp ${btwn_covariables_table} ${Rbtwn_covariables_table}
sed -i 's/#/_/g' ${Rbtwn_covariables_table}
btwn_covariables_table=${Rbtwn_covariables_table}
echo ${btwn_covariables_table}

#This is the basic relatedness graph for the whole data set colored by relationship using Exomechip thresholds
R CMD BATCH "--args ${rawrel_table}  ${nickname} ${overall_title}" ${scripts}/rawrel_relat.R



R CMD BATCH "--args ${rawrel_table} ${btwn_covariables_table} ${nickname} ${overall_title} ${covariablefile} ${colorfile}" ${scripts}/rawrel_imchip.R




#draw maf pdf

maf_table=${qc1_folder}/${nickname}5bySNP.txt

maf_graph=${graph_folder}/${nickname}maf.ps
R CMD BATCH "--args ${maf_table} ${maf_graph} ${overall_title} ${chip}" ${scripts}/maf.R 
#ps2pdf ${maf_graph}


#mv *.pdf ${graph_folder}
for psfile in $(ls *.ps); do ps2pdf ${psfile}; done
rm -f *.ps

echo "Completed qc_pdf.sh: " $(date)
