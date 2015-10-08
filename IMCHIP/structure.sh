#!/bin/bash



echo -e "#########################################################################"
echo -e
echo -e "# GWAS_Pipedream - Standard GWAS QC Pipeline Package"
echo -e "# (c) 2014-2015 JBonnie, WMChen"
echo -e
echo -e "#########################################################################"

nickname=$1
overall_title=$2
chip=$3
covariablevalue=$4

covcol=$((${covariablevalue} + 2))

project_folder=$(pwd)
script_folder=/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Farewell/IMCHIP

structure_folder=${project_folder}/structure
qc1_folder=${project_folder}/QC1
cov=${qc1_folder}/${nickname}.cov

covariablefile=${project_folder}/covariables${covariablevalue}.list
colorfile=${project_folder}/colors.list
hapmap_table=${structure_folder}/hapmappc.txt
Rhapmap_table=${structure_folder}/hapmappcR.txt


imchip_hapmap3=/h4/t1/projects/IMCHIP/Hapmap3
hapmap3_folder=${imchip_hapmap3}
echo "Now beginning structure.sh: "

graphfolder=graphs_cov${covariablevalue}
mkdir ${graphfolder}
mkdir ${structure_folder}
cd ${structure_folder}

#In some cases we need to update immunochip snps to hg19 in order to map properly. Is this one of those cases?
update_hapmap=T

#If we've already done this stuff... why do it again?
if [ "${covariablevalue}" -eq "1" ]; then


echo "Hapmap files are here: " ${hapmap3_folder}

#copy dn6.fam, dn6.bed, dn6.bim, created in qc.bsh, to dn6b in this folder
	awk '{print $1,$2,$3,$4,$5,2}' ${qc1_folder}/${nickname}6.fam > ${nickname}6b.fam
	cp ${qc1_folder}/${nickname}6.bed ${nickname}6b.bed
	cp ${qc1_folder}/${nickname}6.bim ${nickname}6b.bim

#Prepare the files	
if [ "${update_hapmap}" == "F" ]; then
	king -b ${nickname}6b.bed --merlin --prefix ${nickname} > ${nickname}6bmerlin.log
else
	plink --bfile ${nickname}6b --update-map ${hapmap3_folder}/updatemap.txt --make-bed --out ${nickname}6bmap0 --noweb
	plink --bfile ${nickname}6bmap0 --update-map ${hapmap3_folder}/updatechr.txt --update-chr --make-bed --out ${nickname}6bmap1 --noweb
	plink --bfile ${nickname}6bmap1 --update-map ${hapmap3_folder}/updatename.txt --update-name --make-bed --out ${nickname}6bmap --noweb
	king -b ${nickname}6bmap.bed --merlin --prefix ${nickname} > ${nickname}6bmapmerlin.log

fi

gzip ${nickname}.*

echo "Now merging merlin files with hapmap3 files from "${hapmap3_folder}": "
#merge hapmap3 and dn
king -d ${hapmap3_folder}/hapmap3.dat.gz,${nickname}.dat.gz \
     -p ${hapmap3_folder}/hapmap3.ped.gz,${nickname}.ped.gz \
     -m ${hapmap3_folder}/hapmap3.map.gz,${nickname}.map.gz  \
     --callrate 0.95 --prefix merge --autoflip --plink > merge.log


#pca projection?
king -b merge.bed --pca --project --prefix proj

fi
projpc=${structure_folder}/projpc${covariablevalue}.txt
Rprojpc=${structure_folder}/projpc${covariablevalue}R.txt


#DRAW THINGS
field_count=$(awk 'NR==1{print NF}' projpc.ped)
covariablecol=$((${field_count} + ${covcol}))
awk -v cov=${cov} 'BEGIN{while((getline<cov)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' projpc.ped | awk -v covcol=${covariablecol} '{print $1,$2,$3,$4,$5,$6,$7,$8,$covcol}' > ${projpc}

grep "NA" projpc.ped > ${hapmap_table}
#Because there is something funky going on with the pca-analysis, we're doing something different here
#awk '$7 > .002 && $8 > .044 && !($9==1)' projpc.txt > NEdraw.txt


cp ${projpc} ${Rprojpc}
sed -i 's/#/_/g' ${Rprojpc}


R CMD BATCH "--args ${Rprojpc} ${hapmap_table} ${nickname} ${overall_title} ${covariablefile} ${colorfile}" ${scripts}/mergepca.R



### EUROPEAN POP STUCTURE




echo "Tarring files to structure.tar: " $(date)
rmlist=$(ls ${nickname}6b* *TMP* *.hh *.tmp ${nickname}.*)

#tar --create -f structure.tar *
#gzip structure.tar

echo "Removing select files: " $(date)
rm -f $rmlist
for psfile in $(ls *.ps); do ps2pdf ${psfile}; done
rm -f *.ps

mv *.pdf ${project_folder}/${graphfolder}/
echo "Completed structure.sh: " $(date)


