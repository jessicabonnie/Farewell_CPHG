#!/bin/bash
script_folder=/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Farewell/IMCHIP


echo -e "#########################################################################"
echo -e
echo -e "# GWAS_Pipedream - Standard GWAS QC Pipeline Package"
echo -e "# (c) 2014-2015 JBonnie, WMChen"
echo -e
echo -e "#########################################################################"

echo -e "# Disc: Run checks for Hardy Weinberg Disequilibrium on Northern European Population. Produces SNP list to be added to SNP removal list as well as graphs for visual representation of the population used."


nickname=$1
overall_title=$2
chip=$3
covariablevalue=$4
#covcol=$((${covariablevalue} + 2))
project_folder=$(pwd)
#log=${project_folder}/data_qc.log
qc2_folder=${project_folder}/QC2_HWE
structure_folder=${project_folder}/structure
graph_folder=${project_folder}/graphs_cov${covariablevalue}

qc1_folder=${project_folder}/QC1
nb_folder=${project_folder}/NB

imchip_hapmap3=/h4/t1/projects/IMCHIP/Hapmap3
hapmap_table=${structure_folder}/hapmappc.txt
#projpc_table=${structure_folder}/projpc.txt

cov=${qc1_folder}/${nickname}.cov
covariablefile=${project_folder}/covariables${covariablevalue}.list
colorfile=${project_folder}/colors.list
pca_odd=F
hwe_remove=F
hapmap3_folder=${imchip_hapmap3}

mkdir ${qc2_folder}
cd ${qc2_folder}

AAslope=1.3
AAint=0.03
NElist_loc=${qc2_folder}/NElist_alt.txt
NElist_covariables=${qc2_folder}/NElist_alt_covariables${covariablevalue}.nb
NOT_NElist_loc=${qc2_folder}/NOT_NElist_alt.txt
NOT_NElist_covariables=${qc2_folder}/NOT_NElist_alt_covariables${covariablevalue}.nb


AAlist_loc=${qc2_folder}/AAlist_alt.txt
AAlist_covariables=${qc2_folder}/AAlist_alt_covariables${covariablevalue}.nb
Europeanlist_loc=${qc2_folder}/Europeanlist_alt.txt
Europeanlist_covariables=${qc2_folder}/Europeanlist_alt_covariables${covariablevalue}.nb

awk '$8 > .043' ${structure_folder}/projpc${covariablevalue}.txt | awk '{print $1,$2}' > ${NElist_loc}
awk '$8 > .043' ${structure_folder}/projpc${covariablevalue}.txt | awk '{print $1,$2,$NF}' > ${NElist_covariables}
awk '$8 < .043' ${structure_folder}/projpc${covariablevalue}.txt | awk '{print $1,$2}' > ${NOT_NElist_loc}
awk '$8 < .043' ${structure_folder}/projpc${covariablevalue}.txt | awk '{print $1,$2,$NF}' > ${NOT_NElist_covariables}
    
awk '$8 > .043' ${structure_folder}/projpc${covariablevalue}.txt | awk '{print $1,$2}' > ${Europeanlist_covariables}


    
awk -v b=${AAint} -v m=${AAslope} '$7 < 0.0 && $8 < .043 && $8 > (m*$7 + b)' ${structure_folder}/projpc${covariablevalue}.txt | awk '{print $1,$2}' >${AAlist_loc}
awk -v b=${AAint} -v m=${AAslope} '$7 < 0.0 && $8 < .043 && $8 > (m*$7 + b)' ${structure_folder}/projpc${covariablevalue}.txt | awk '{print $1,$2,$NF}' >${AAlist_covariables}
awk -v b=${AAint} -v m=${AAslope} '$7 < 0.0 && $8 < .043 && $8 < (m*$7 + b)' ${structure_folder}/projpc${covariablevalue}.txt | awk '{print $1,$2,$NF}' >> ${nb_folder}/AA_exclude.nb

    #awk '$8 > .045' ${structure_folder}/projpc.txt > NEdraw.txt
    
HW_list=${NElist_loc}
HW_covariables=${NElist_covariables}

if [ "${covariablevalue}" -eq "1" ]; then
  plink --bfile ${qc1_folder}/${nickname}6 --keep ${HW_list} --make-bed --out HW1 --noweb
  king -b HW1.bed --unrelated --prefix HW1 > HW1unrelated.log


#Doing HWE on controls only, so now we don't need to worry about HLA
  hw_threshold=1e-6


  if [ "${hwe_remove}" = "T" ]; then
    plink --bfile HW1 --keep HW1unrelated.txt --hwe ${hw_threshold} --filter-controls --hardy --noweb --out HW2 --remove ${project_folder}/hwe_remove.list 
  else
    plink --bfile HW1 --keep HW1unrelated.txt --hwe ${hw_threshold} --filter-controls --hardy --noweb --out HW2
  fi

  echo $(head -n1 HW2.hwe) "POS"  > HW2pos.hwe
  awk 'BEGIN{while((getline<"HW1.bim")>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' HW2.hwe | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13}' >> HW2pos.hwe
  grep "ALL" HW2.hwe | awk -v hwt=${hw_threshold} '$NF < hwt' | awk '{print $2}' > hweSNP.txt
  grep "ALL" HW2pos.hwe | awk -v hwt=${hw_threshold} '$9 < hwt' | awk '{print $1,$2,$9,$NF}' > hweSNP.nb
  grep "ALL" HW2pos.hwe > HW2_plinkALL_hwe.nb
  grep "ALL" HW2.hwe > HW2_plinkALL.hwe 

fi
echo ${NElist_loc}
echo ${NOT_NElist_loc}

RNElist_loc=${NElist_loc}_R
RNOT_NElist_loc=${NOT_NElist_loc}_R

RAAlist_loc=${AAlist_loc}_R
RHW_list=${HW_list}_R
cp ${NElist_loc} ${RNElist_loc}
sed -i 's/#/_/g' ${RNElist_loc}
cp ${NOT_NElist_loc} ${RNOT_NElist_loc}
sed -i 's/#/_/g' ${RNOT_NElist_loc}

cp ${AAlist_covariables} ${RAAlist_loc}
sed -i 's/#/_/g' ${RAAlist_loc}


cp ${HW_covariables} ${RHW_list}
sed -i 's/#/_/g' ${RHW_list}

Rprojpc=${structure_folder}/projpc${covariablevalue}R.txt


R CMD BATCH "--args ${Rprojpc} ${hapmap_table} ${RNElist_loc} ${nickname} ${overall_title} ${covariablefile} ${colorfile} Europeans" ${script_folder}/NE1EURmergepca.R
R CMD BATCH "--args ${Rprojpc} ${hapmap_table} ${RNOT_NElist_loc} ${nickname} ${overall_title} ${covariablefile} ${colorfile} Non-Europeans" ${script_folder}/NE1EURmergepca.R

R CMD BATCH "--args ${Rprojpc} ${hapmap_table} ${RHW_list} ${nickname} ${overall_title} ${covariablefile} ${colorfile} HardyWeinbergPopulation" ${script_folder}/NE1EURmergepca.R



for psfile in $(ls *.ps); do ps2pdf ${psfile}; done

rm -f *.ps

mv *.pdf ${graph_folder}
