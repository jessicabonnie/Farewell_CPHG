#!/bin/bash

#release=${outer_project_folder}/release_24Sept2015

release=$(pwd)
kinship=${release}/kinship

outer_project_folder=/t121/hellacious/projects/HELLACIOUS/workfiles/Jessica/
scripts=${outer_project_folder}/scripts




#### NOTE ##### THESE FILES MAY HAVE BEEN MOVED #### THEY ARE THE POST ZCALL RELEASE FILES FROM DATAFREEZE 1 and 2

 wm_12v1=/t121/hellacious/projects/HELLACIOUS/workfiles/Wei-Min/release_23Sept2015/12-1.0/QC_Z
 wm_12v11=/t121/hellacious/projects/HELLACIOUS/workfiles/Wei-Min/release_23Sept2015/12-1.1/QC_Z
 wm_24v1=/t121/hellacious/projects/HELLACIOUS/workfiles/Wei-Min/release_23Sept2015/24-1.0/QC_Z
 
 ### NOTE ### THESE FILES MAY HAVE BEEN MOVED ### THEY ARE THE ZCALLED PLINK FILES FROM DATAFREEZE 1 and 2
 jb_12v1=/t121/hellacious/projects/HELLACIOUS/workfiles/Jessica/dataFreeze1/release_23Sept2015/V1.0/zCalled/raw_zcall
 jb_12v11=/t121/hellacious/projects/HELLACIOUS/workfiles/Jessica/dataFreeze1/release_23Sept2015/V1.1/zCalled/raw_zcall
 jb_24v1=/t121/hellacious/projects/HELLACIOUS/workfiles/Jessica/dataFreeze2/release_24Sept2015/zCalled/raw_zcall
 
 ### PHENOTYPE FILES ###
 pheno12v1=/t121/hellacious/projects/HELLACIOUS/HumanCoreExome_12v1_HELLAC_Upload/phenofiles/HELLAC_DN_array_v1_phenotype_20150319_v2.txt
pheno12v11=/t121/hellacious/projects/HELLACIOUS/HumanCoreExome_12v1.1_HELLAC_Upload/phenofiles/HELLAC_DN_array_v1.1_phenotype_20150319.txt
pheno24v1=/t121/hellacious/projects/HELLACIOUS/HumanCoreExome_24v1_HELLACIOUS_Freeze2/HELLAC_Freeze2_Phenotype_20150717_query_v2.txt
 
mkdir ${kinship} 

###########################
#
#   Kinship
#
###############################

#Navigate to folder
cd ${kinship}

 ### Update all of the plink files with given lists

plink --noweb --make-bed --bfile ${jb_24v1} --update-sex ${wm_24v1}/updatesex.txt --remove ${wm_24v1}/sampletoberemoved.txt --exclude ${wm_24v1}/snptoberemoved.txt --out ${kinship}/24v1clean --allow-no-sex


plink --noweb --make-bed --bfile ${jb_12v1} --update-sex ${wm_v1}/updatesex.txt --remove ${wm_v1}/sampletoberemoved.txt --exclude ${wm_v1}/snptoberemoved.txt --out ${kinship}/12v1.0clean --allow-no-sex

# remove duplicate samples from chip 12v1.1 (which are also in chip 12v1.0)
cat  ${jb_12v1}.fam ${jb_12v11}.fam | awk 'BEGIN { FS = "\t" }{print $1, $2}' | sort -k2 | uniq -f1 -d  > ${kinship}/purposeful_duplicates_remove.tmp
cat ${wm_v11}/sampletoberemoved.txt ${kinship}/purposeful_duplicates_remove.tmp > ${kinship}/allsamplestoremove_fromv11.txt

plink --noweb --make-bed --bfile ${jb_12v11} --update-sex ${wm_v11}/updatesex.txt --remove ${kinship}/allsamplestoremove_fromv11.txt --exclude ${wm_v11}/snptoberemoved.txt --out ${kinship}/12v1.1clean --allow-no-sex


 ### CREATE A SAMPLECOHORT FILE with ALL SAMPLES with their COHORT/FREEZE/CHIP information for use in DIVIDING THE RELEASE

 ### There is something wrong with XCohortX's sample IDs, so we need to fix them
 sed 's/\r//g' ${pheno12v1}|awk 'BEGIN { FS = "\t" } NR>1{print $2,$7,$6,$9,$10,$12}'  > ${kinship}/samplescohorts12v1_plus.tmp
 
 grep -v SubStudyX ${kinship}/samplescohorts12v1_plus.tmp | awk '{print $1,$2,$3}' > ${kinship}/samplescohorts12v1.tmp
 
 grep SubStudyX ${kinship}/samplescohorts12v1_plus.tmp  | awk ' {sub(/SN_/, "", $3) }1;' |awk '{print $1,$2,$3}' >> ${kinship}/samplescohorts12v1.tmp


 sed 's/\r//g' ${pheno12v11} |awk 'BEGIN { FS = "\t" } NR>1{print $2,$7,$5}' > ${kinship}/samplescohorts12v11.tmp
 sed 's/\r//g' ${pheno24v1} |awk 'BEGIN { FS = "\t" } NR>1{print $2,$7,$5}' > ${kinship}/samplescohorts24v1.0.tmp

echo "FID IID Cohort SampleID Freeze Chip" > ${kinship}/samplescohorts.txt
 awk -v fam=${jb_12v1}.fam 'BEGIN{while((getline<fam)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${kinship}/samplescohorts12v1.tmp | awk '{print $4,$5,$1,$3,"1","12v1.0"}' >> ${kinship}/samplescohorts.txt

awk -v fam=${jb_12v11}.fam 'BEGIN{while((getline<fam)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${kinship}/samplescohorts12v11.tmp | awk '{print $4,$5,$1,$3,"1","12v1.1"}' >> ${kinship}/samplescohorts.txt

awk -v fam=${jb_24v1}.fam 'BEGIN{while((getline<fam)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${kinship}/samplescohorts24v1.0.tmp | awk '{print $4,$5,$1,$3,"2","24v1.0"}' >> ${kinship}/samplescohorts.txt


 ### CREATE A SAMPLECOHORT FILE with ALL CLEAN SAMPLES and their cohort/freeze/chip information for use in dividing the RELATIONSHIP FILES
echo "FID IID Cohort SampleID Freeze Chip" > ${kinship}/clean_samplescohorts.txt
awk -v fam=${kinship}/12v1.0clean.fam 'BEGIN{while((getline<fam)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${kinship}/samplescohorts12v1.tmp | awk '{print $4,$5,$1,$3,"1","12v1.0"}' >> ${kinship}/clean_samplescohorts.txt

awk -v fam=${kinship}/12v1.1clean.fam 'BEGIN{while((getline<fam)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${kinship}/samplescohorts12v11.tmp | awk '{print $4,$5,$1,$3,"1","12v1.1"}' >> ${kinship}/clean_samplescohorts.txt

awk -v fam=${kinship}/24v1.0clean.fam 'BEGIN{while((getline<fam)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${kinship}/samplescohorts24v1.0.tmp | awk '{print $4,$5,$1,$3,"2","24v1.0"}' >> ${kinship}/clean_samplescohorts.txt




#This function takes a chip and a list of chips to loop through in order to create all necessary kinship files involving that chip. It also draws relatedness graphs between chips.
relationship (){
    chipversion=$1
    chiplist=$2
    index=0
    for chip in ${chiplist[@]}
    do
        if [ $chip == $chipversion ]; then
            if [[ ! -s ${kinship}/${chipversion}.kin0 ]]; then
             king -b ${kinship}/${chipversion}clean.bed --related --degree 3 --prefix ${kinship}/${chipversion}
            fi
            echo $(head -n1 ${kinship}/${chipversion}.kin0) "Cohort1 SampleID1 Freeze_ID1 Chip_ID1 Cohort2 SampleID2 Freeze_ID2 Chip_ID2"  > ${kinship}/${chipversion}kin0_plus.txt
            
            LANG=en_EN join -1 2 -2 2 <(LANG=en_EN sort -k2 ${kinship}/${chipversion}.kin0) <(LANG=en_EN sort -k2 ${kinship}/clean_samplescohorts.txt) | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13}' > ${kinship}/${chipversion}kin0_plus.tmp
            echo "1"

            LANG=en_EN join -1 4 -2 2 <(LANG=en_EN sort -k4 ${kinship}/${chipversion}kin0_plus.tmp) <(LANG=en_EN sort -k2 ${kinship}/clean_samplescohorts.txt) | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$14,$15,$16,$17}' >> ${kinship}/${chipversion}kin0_plus.txt
            echo "2"
            head -n1  ${kinship}/${chipversion}kin0_plus.txt | awk '{print $0,"Relationship"}'>  ${kinship}/${chipversion}relationship.txt
            
            awk 'NR > 1 && $8 > 0.4' ${kinship}/${chipversion}kin0_plus.txt| awk '{print $0, "Duplicate"}' >>  ${kinship}/${chipversion}relationship.txt
            
            awk 'NR > 1 && $8 < 0.4 && $8 > 0.177 && $7<0.005' ${kinship}/${chipversion}kin0_plus.txt | awk '{print $0, "PO"}' >>  ${kinship}/${chipversion}relationship.txt
            
            awk 'NR > 1 && $8 < 0.4 && $8 > 0.177 && $7>0.005' ${kinship}/${chipversion}kin0_plus.txt | awk '{print $0, "FS"}' >>  ${kinship}/${chipversion}relationship.txt
            
            awk 'NR > 1 && $8 <=0.177 && $8 > 0.0884' ${kinship}/${chipversion}kin0_plus.txt| awk '{print $0, "2nd"}' >>  ${kinship}/${chipversion}relationship.txt
            
            awk '$9!=$13' ${kinship}/${chipversion}relationship.txt > ${kinship}/${chipversion}_relat_btwn_cohorts_plus.txt
            
            


        else
            if [[ ! -s ${kinship}/${chipversion}with${chip}_degree3.log ]]; then
               king -b ${kinship}/${chipversion}clean,${kinship}/${chip}clean --related --degree 3 --prefix ${kinship}/${chipversion}with${chip}_d3 > ${kinship}/${chipversion}with${chip}_degree3.log
            fi
            echo $(head -n1 ${kinship}/${chipversion}with${chip}_d3.kin0) "Cohort1 SampleID1 Freeze_ID1 Chip_ID1 Cohort2 SampleID2 Freeze_ID2 Chip_ID2"  > ${kinship}/${chipversion}with${chip}kin0_plus.txt

            LANG=en_EN join -1 2 -2 2 <(LANG=en_EN sort -k2 ${kinship}/${chipversion}with${chip}_d3.kin0) <(LANG=en_EN sort -k2 ${kinship}/clean_samplescohorts.txt) | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13}' > ${kinship}/${chipversion}with${chip}kin0_plus.tmp
            echo "3"

            LANG=en_EN join -1 4 -2 2 <(LANG=en_EN sort -k4 ${kinship}/${chipversion}with${chip}kin0_plus.tmp) <(LANG=en_EN sort -k2 ${kinship}/clean_samplescohorts.txt) | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$14,$15,$16,$17}' >> ${kinship}/${chipversion}with${chip}kin0_plus.txt
            echo "4"

            #Define Relationships in kin0
            head -n1  ${kinship}/${chipversion}with${chip}kin0_plus.txt | awk '{print $0,"Relationship"}'>  ${kinship}/relationship_${chipversion}with${chip}.txt
            
            awk 'NR > 1 && $8 > 0.4' ${kinship}/${chipversion}with${chip}kin0_plus.txt| awk '{print $0, "Duplicate"}' >>  ${kinship}/relationship_${chipversion}with${chip}.txt
            
            awk 'NR > 1 && $8 < 0.4 && $8 > 0.177 && $7<0.005' ${kinship}/${chipversion}with${chip}kin0_plus.txt | awk '{print $0, "PO"}' >>  ${kinship}/relationship_${chipversion}with${chip}.txt
            
            awk 'NR > 1 && $8 < 0.4 && $8 > 0.177 && $7>0.005' ${kinship}/${chipversion}with${chip}kin0_plus.txt | awk '{print $0, "FS"}' >>  ${kinship}/relationship_${chipversion}with${chip}.txt
            
            awk 'NR > 1 && $8 <=0.177 && $8 > 0.0884' ${kinship}/${chipversion}with${chip}kin0_plus.txt| awk '{print $0, "2nd"}' >>  ${kinship}/relationship_${chipversion}with${chip}.txt


            awk '$9!=$13' ${kinship}/relationship_${chipversion}with${chip}.txt > ${kinship}/${chipversion}with${chip}_relat_btwn_cohorts_plus.txt
            
            
            ### Draw Graphs
            cd ${kinship}
            R CMD BATCH "--args ${kinship}/${chipversion}with${chip}_d3.kin0 ${chipversion}and${chip}_ALL HELLAC-DN_-_ALL_in_${chipversion}_and_${chip}" ${scripts}/rawrel_relat.R
            
            awk '$11!=$15' ${kinship}/${chipversion}with${chip}kin0_plus.txt > ${kinship}/${chipversion}with${chip}_graph.tmp
            
            R CMD BATCH "--args ${kinship}/${chipversion}with${chip}_graph.tmp ${chipversion}and${chip} HELLAC-DN_between_Chips_${chipversion}_and_${chip}" ${scripts}/rawrel_relat.R
            
        fi
    
    done

}

 
#################################
#
#   Create Relationship Files
#
##################################


chiplist=( 24v1.0 12v1.0 12v1.1 )
relationship ${chiplist[0]} ${chiplist}
relationship ${chiplist[1]} ${chiplist[@]:1}
relationship ${chiplist[2]} ${chiplist[@]:2}


for psfile in $(ls ${kinship}/*.ps); do ps2pdf ${psfile} ${psfile/%ps/pdf}; done

rm -f ${kinship}/*.tmp ${kinship}/*.hh ${kinship}/*.nosex ${kinship}/*.ps 



 ### Create Master Lists ###
 
 
awk '$11!=$15' ${kinship}/relationship_24v1.0with12v1.1.txt > ${kinship}/24v1.0with12v1.1_relat_btwn_freezes_plus.txt
awk '$11!=$15' ${kinship}/relationship_24v1.0with12v1.0.txt > ${kinship}/24v1.0with12v1.0_relat_btwn_freezes_plus.txt

awk '$12!=$16' ${kinship}/relationship_12v1.0with12v1.1.txt > ${kinship}/12v1.0with12v1.1_relat_btwn_chips_plus.txt
awk '$12!=$16' ${kinship}/relationship_24v1.0with12v1.0.txt > ${kinship}/24v1.0with12v1.0_relat_btwn_chips_plus.txt
awk '$12!=$16' ${kinship}/relationship_24v1.0with12v1.1.txt > ${kinship}/24v1.0with12v1.1_relat_btwn_chips_plus.txt

 
 
 cat ${kinship}/24v1.0relationship.txt <(awk 'NR>1{print $0}' ${kinship}/24v1.0with12v1.0_relat_btwn_freezes_plus.txt)  <(awk 'NR>1{print $0}' ${kinship}/24v1.0with12v1.1_relat_btwn_freezes_plus.txt)  | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$12,$13,$15,$16,$17}' > ${kinship}/all_new_unreported_relationships.txt
 
  cat ${kinship}/24v1.0relationship.txt \
  <(awk 'NR>1{print $0}' ${kinship}/12v1.0relationship.txt) \
  <(awk 'NR>1{print $0}' ${kinship}/12v1.1relationship.txt) \
  <(awk 'NR>1{print $0}' ${kinship}/24v1.0with12v1.0_relat_btwn_chips_plus.txt) \
  <(awk 'NR>1{print $0}' ${kinship}/24v1.0with12v1.1_relat_btwn_chips_plus.txt) \
  <(awk 'NR>1{print $0}' ${kinship}/12v1.0with12v1.1_relat_btwn_chips_plus.txt) | \
  awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$12,$13,$15,$16,$17}' > ${kinship}/all_unreported_relationships.txt


 head -n1 ${kinship}/all_unreported_relationships.txt > ${kinship}/relationship_btwn_cohortskin0.txt
 awk '$9!=$12' ${kinship}/all_unreported_relationships.txt >> ${kinship}/relationship_btwn_cohortskin0.txt
 
 head -n1 ${kinship}/all_unreported_relationships.txt > ${kinship}/relationship_within_cohortskin0.txt
 awk '$9==$12' ${kinship}/all_unreported_relationships.txt >> ${kinship}/relationship_within_cohortskin0.txt
 
 
 ### Related Samples in 12v1.0
 
 echo $(head -n1 ${kinship}/12v1.1.kin) "Relationship"  > ${kinship}/12v1.1_reported.txt
 awk 'NR > 1 && $9 > 0.4' ${kinship}/12v1.1.kin | awk '{print $0, "Duplicate"}' >> ${kinship}/12v1.1_reported.txt      
 awk 'NR > 1 && $9 < 0.4 && $9 > 0.177 && $8<0.005' ${kinship}/12v1.1.kin | awk '{print $0, "PO"}' >> ${kinship}/12v1.1_reported.txt
 awk 'NR > 1 && $9 < 0.4 && $9 > 0.177 && $8>0.005' ${kinship}/12v1.1.kin | awk '{print $0, "FS"}' >>  ${kinship}/12v1.1_reported.txt
 awk 'NR > 1 && $9 <=0.177 && $9 > 0.0884' ${kinship}/12v1.1.kin| awk '{print $0, "2nd"}' >>  ${kinship}/12v1.1_reported.txt
 
 
 
#################################
#
#   Create Release: Split into Cohorts
#
##################################



chiplist=(24v1.0 12v1.0 12v1.1)
for chip in ${chiplist[@]}
  do
    release_chip=${release}/${chip}
    split=${release_chip}
    mkdir ${split}
    cd ${split}
    
    if [ ${chip} == "24v1.0" ]; then
    
    wmfolder=${wm_24v1}
    pheno=${pheno24v1}
    jbfolder=${jb_24v1}
    
    elif [ ${chip} == "12v1.0" ]; then
    
    wmfolder=${wm_12v1}
    pheno=${pheno12v1}
    jbfolder=${jb_12v1}
    
    elif [ ${chip} == "12v1.1" ]; then
    
    wmfolder=${wm_12v11}
    pheno=${pheno12v11}
    jbfolder=${jb_12v11}
    
    else
    echo "I HAVE NO CHIP VERSION!"
    exit
    
    fi

samplesremove=${wmfolder}/sampletoberemoved.txt
snpsremove=${wmfolder}/snptoberemoved.txt
sexupdate=${wmfolder}/updatesex.txt
samplescohorts=${kinship}/samplescohorts.txt
#plinkbase=${jb_24v1}



cd ${release_chip}
yes | cp ${samplesremove} ${release_chip}/sampletoberemoved.txt
samplesremove=${release_chip}/sampletoberemoved.txt
yes | cp ${snpsremove} ${release_chip}/snpstoberemoved.txt
snpsremove=${release_chip}/snpstoberemoved.txt

yes | cp ${sexupdate} ${release_chip}/updatesex.txt
sexupdate=${release_chip}/updatesex.txt


# make sample id key
grep ${chip} ${samplescohorts} | awk '{print $4,$2,$3}'  > ${release_chip}/SampleIDKey.txt

echo "FID IID Cohort" > ${split}/samplescohorts.txt
grep ${chip} ${samplescohorts} | awk '{print $1,$2,$3}' >> ${split}/samplescohorts.txt

cd ${split}

#sed 's/\r//g' ${pheno} |awk 'BEGIN { FS = "\t" } {print $2,$7}' > ${split}/samplescohorts.tmp

#echo "FID IID Cohort" > ${split}/samplescohorts.txt

#awk -v fam=${jbfolder}.fam 'BEGIN{while((getline<fam)>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${split}/samplescohorts.tmp | awk '{print $3,$4,$1}' >> ${split}/samplescohorts.txt



###################


 ### Determine List of Cohorts
awk 'NR>1{print $3}' ${split}/samplescohorts.txt| sed '/^\s*$/d' |sort | uniq > ${split}/cohorts.list
cohorts=($( <${split}/cohorts.list))
count=$(echo ${#cohorts[@]})

#include cohort information in removal and update lists
head -n1 ${samplesremove} | awk '{print $1,$2,$3,"Cohort"}'> ${split}/samplesremove_cohorts.txt
awk -v scohort=${split}/samplescohorts.txt 'BEGIN{while((getline<scohort)>0)l[$2]=$0}$2 in l{print $0" "l[$2]}' ${samplesremove} | awk '{print $1,$2,$3,$6}' >> ${split}/samplesremove_cohorts.txt

echo "FID IID Sex Cohort" > ${split}/updatesex_cohorts.txt
awk -v scohort=${split}/samplescohorts.txt 'BEGIN{while((getline<scohort)>0)l[$2]=$0}$2 in l{print $0" "l[$2]}' ${sexupdate} | awk '{print $1,$2,$3,$6}' >> ${split}/updatesex_cohorts.txt


# split the lists and the plink files and the relationship files
for i in $(seq 0 $((${count}-1))); do
#     if [ "${cohorts[$i]}" = "Joslin_Trios" ]; then
    mkdir ${split}/${cohorts[$i]}
    grep -w ${cohorts[$i]} ${split}/samplescohorts.txt | awk '{print $1,$2}' > ${split}/${cohorts[$i]}/${cohorts[$i]}_samples.txt
    grep -w ${cohorts[$i]} ${split}/samplesremove_cohorts.txt | awk '{print $1,$2,$3}' > ${split}/${cohorts[$i]}/${cohorts[$i]}_samplesremove.txt
    plink --noweb --keep ${split}/${cohorts[$i]}/${cohorts[$i]}_samples.txt --bfile ${jbfolder} --make-bed --out ${split}/${cohorts[$i]}/${cohorts[$i]}_raw

    head -n1 ${kinship}/relationship_btwn_cohortskin0.txt > ${split}/${cohorts[$i]}/${cohorts[$i]}_relationships_btwn_cohorts.txt
    grep -w ${cohorts[$i]} ${kinship}/relationship_btwn_cohortskin0.txt >> ${split}/${cohorts[$i]}/${cohorts[$i]}_relationships_btwn_cohorts.txt

    head -n1 ${kinship}/relationship_within_cohortskin0.txt > ${split}/${cohorts[$i]}/${cohorts[$i]}_unreported_relationships_within_cohort.txt
    grep -w ${cohorts[$i]} ${kinship}/relationship_within_cohortskin0.txt >> ${split}/${cohorts[$i]}/${cohorts[$i]}_unreported_relationships_within_cohort.txt

    
    if [ "${cohorts[$i]}" = "EDIC" ]; then
        
        yes | cp ${kinship}/12v1.1_reported.txt ${cohorts[$i]}/${cohorts[$i]}_reported_relationships_within_cohort.txt

    fi

    
    

    #Make Study SampleID Key
    

    echo "StudySampleID UVASampleID" > ${split}/${cohorts[$i]}/${cohorts[$i]}_SampleIDKey.txt
    grep -w ${cohorts[$i]} ${release_chip}/SampleIDKey.txt| awk '{print $1,$2}' >> ${split}/${cohorts[$i]}/${cohorts[$i]}_SampleIDKey.txt
    
    #Make study update sex file
     grep -w ${cohorts[$i]} ${split}/updatesex_cohorts.txt | awk '{print $1,$2,$3}' >> ${split}/${cohorts[$i]}/${cohorts[$i]}_updatesex.txt
    
    #Copy SNP removal list into folder
    yes | cp ${snpsremove} ${split}/${cohorts[$i]}/.
    
    rm -f ${split}/${cohorts[$i]}/*.nosex ${split}/${cohorts[$i]}/*.hh ${split}/${cohorts[$i]}/*.nof
#     fi
done



rm -f ../*.tmp rm *.tmp
done

