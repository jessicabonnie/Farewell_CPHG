#!/bin/bash

chip=$1

#Results will be written in folders on parent and aunt levels (above the current directory and lateral to that parent directory)because that is where wei-min's scripts write them

## Thus it is recommended that you create a "QC_temp folder and then create a "tempdata" folder inside of that one. Then go there.


 ### Human Boundaries will need to be manually changed (obviously)
 ### The Step12 and Step34 scripts have HARD CODED PATHS IN THEM
 ### HARD CODED JBSCRIPT FOLDER HERE:
jbscript=/t121/hellacious/projects/HELLACIOUS/workfiles/Jessica/scripts 
 ### HARD CODED WMC SCRIPT FOLDER HERE:
 wmscript=/t121/hellacious/projects/HELLACIOUS/workfiles/Wei-Min/RanyQC/QC_protocol/SHcode
 
 ### HARD CODED PATH TO HAPMAP POPULATION FILE
 
 popfile=/t121/hellacious/projects/HELLACIOUS/workfiles/Jessica/Reference/relationships_w_pops_041510.txt
 
 tmpfolder=$(pwd)
 # establish locations of output folders for wei-min's scripts
final_folder=${tmpfolder}/../../QC_final
qctemp=${tmpfolder}/..


## automatic boundary lists:

 autolist_24v1=( Cohort1 Cohort2 Cohort3)
 autolist_12v1=( Cohort4 Cohort5 Cohort6 Cohort7 Cohort8 )
 autolist_12v11=( Cohort9 Cohort10 Cohort11 Cohort12 )

draw_boundaries(){

    study=$1
    chip=$2
    flag=$3
    xmin=$4
    xmax=$5
    ymin=$6
    ymax=$7
    studytmp=${qctemp}/${study}
    studyfinal=${final_folder}/${study}
    
    
    projfile=${studytmp}/${study}_hapmap3pc.ped
    rmlist=${studytmp}/postPCA_removallist.txt
    
    echo "FID IID" > $rmlist 
    awk -v xmin=${xmin} -v xmax=${xmax} '$6==2 && ($7<xmin || $7 > xmax){print$1,$2}' $projfile >> $rmlist
    awk -v xmin=${xmin} -v xmax=${xmax} -v ymin=${ymin} -v ymax=${ymax} '$6==2 && ($7>xmin && $7 < xmax) && ($8<ymin || $8>ymax){print$1,$2}' $projfile >> $rmlist
    
    LANG=en_EN join -1 2 -2 1  <(awk '$6 == 1' $projfile | LANG=en_EN sort -k2) <(awk '{print $2,$7}' ${popfile} | LANG=en_EN sort -k1 )  > ${studytmp}/${study}_hapmap3pc.txt
    awk -v study=${study} '$6 == 2 {print $0, study}' $projfile >> ${studytmp}/${study}_hapmap3pc.txt

    R CMD BATCH "--args ${studytmp}/${study}_hapmap3pc.txt $study $chip 27 ${flag} ${xmin} ${xmax} ${ymin} ${ymax}" ${jbscript}/PCA_PLOTTER_JB.R
    
    for psfile in $(ls ${studytmp}/*.ps)
    do 
        ps2pdf ${psfile} ${psfile/%ps/pdf}
        yes | cp ${psfile/%ps/pdf} .
        rm -f ${psfile}
    done

    }


draw_scatters(){

    study=$1
    chip=$2
    studytmp=${qctemp}/${study}
    studyfinal=${final_folder}/${study}
    
    
    unprojfile=${studytmp}/${study}_finalpca_tmp_pc.ped
    rmlist=${studytmp}/postPCA_removallist.txt

    R CMD BATCH "--args ${studytmp}/${study}_hapmap3pc.txt $study $chip TRUE ${studytmp}/${study}_hapmap_rooted_PCA_PLOTS" ${jbscript}/PCA_SCATTER_JB_grid.R   
    
    grep -v -f ${rmlist} ${studytmp}/${study}_hapmap3pc.txt > ${studytmp}/${study}_hapmap3pc_clean.txt
    R CMD BATCH "--args ${studytmp}/${study}_hapmap3pc_clean.txt $study $chip TRUE ${studytmp}/${study}_hapmap_rooted_PCA_PLOTS_CLEANED" ${jbscript}/PCA_SCATTER_JB_grid.R
    
    
    ceutsi=${studytmp}/${study}_hapmap3_europeanpc.ped
    LANG=en_EN join -1 2 -2 1  <(awk '$6 == 1' $ceutsi | LANG=en_EN sort -k2) <(awk '{print $2,$7}' ${popfile} | LANG=en_EN sort -k1 )  > ${studytmp}/${study}_hapmap3_europeanpc.txt
    awk -v study=${study} '$6 == 2 {print $0, study}' $ceutsi >> ${studytmp}/${study}_hapmap3_europeanpc.txt
    
    #R CMD BATCH "--args ${studytmp}/${study}_hapmap3_europeanpc.txt $study $chip TRUE ${studytmp}/${study}_CEU_TSI_rooted_PCA_PLOTS_JB_wrap" ${jbscript}/PCA_SCATTER_JB_wrap.R
    #R CMD BATCH "--args ${studytmp}/${study}_hapmap3_europeanpc.txt $study $chip TRUE ${studytmp}/${study}_CEU_TSI_rooted_PCA_PLOTS" ${jbscript}/PCA_SCATTER_JB_grid.R
    

    #what are these xun samples? I don't want to draw them, right??
    awk '$6==2' $unprojfile > ${studytmp}/${study}_finalpca_tmp_pc.txt.tmp
    
    LANG=en_EN join -1 2 -2 2  <(LANG=en_EN sort -k2,2 ${studytmp}/${study}_finalpca_tmp_pc.txt.tmp) <( LANG=en_EN sort -k2,2 ${rmlist} ) | cut -d" " -f1-26 | awk '{print $0, "Removed"}' > ${studytmp}/${study}_finalpca_tmp_pc.txt
    
    grep -v -f ${rmlist} ${studytmp}/${study}_finalpca_tmp_pc.txt.tmp  | awk '{print $0, "Kept"}' >> ${studytmp}/${study}_finalpca_tmp_pc.txt
    
    R CMD BATCH "--args ${studytmp}/${study}_finalpca_tmp_pc.txt ${study} ${chip} ${studytmp}/${study}_Unrooted_PCA_PLOTS_CLEANED" ${jbscript}/PCA_SCATTER_PLAIN.R
    
    grep -v -f ${rmlist} ${studytmp}/${study}_hapmap3_europeanpc.txt > ${studytmp}/${study}_hapmap3_europeanpc_clean.txt
    R CMD BATCH "--args ${studytmp}/${study}_hapmap3_europeanpc_clean.txt $study $chip TRUE ${studytmp}/${study}_CEU_TSI_rooted_PCA_PLOTS_CLEANED" ${jbscript}/PCA_SCATTER_JB_grid.R
    
    for psfile in $(ls ${studytmp}/*.ps)
    do 
        ps2pdf ${psfile} ${psfile/%ps/pdf}
        yes | cp ${psfile/%ps/pdf} .
        rm -f ${psfile}
    done

    }


completely_automatic(){
    chipv=$1
    study=$2
    studytmp=${qctemp}/${study}
    studyfinal=${final_folder}/${study}

    echo -e
    echo -e "Starting Step12 for $study on chip ${chipv}. Writing intermediate files here: ${studytmp}"
     bash ${wmscript}/Step12_HELLAC_QC_shared_v5.sh ${study} ${chipv} 

    
    echo -e
    echo "Now starting step34.  Writing final result files here: ${studyfinal}"
    bash ${wmscript}/Step34_HELLAC_QC_shared_v5.sh ${study} ${chipv}
    
    echo -e
    echo "Now drawing separate PCA graphs."
    draw_boundaries ${study}  ${chipv} "auto_bound" $(cat ${studytmp}/${study}_PC_boundary.txt)
    
    echo -e
    echo "Now drawing scatter plots."
    draw_scatters ${study}  ${chipv}
    
    
    
    yes | cp ${study}_projpca_auto_bound.pdf ${studyfinal}/${study}_hapmap_rooted_PCA_PLOTS.pdf
    
    #rm -f ${studyfinal}/*JB*.pdf
    yes | cp ${studytmp}/${study}_hapmap_rooted_PCA_PLOTS_CLEANED.pdf ${studyfinal}/.
    yes | cp ${studytmp}/${study}_CEU_TSI_rooted_PCA_PLOTS_CLEANED.pdf ${studyfinal}/.
    yes | cp ${studytmp}/${study}_Unrooted_PCA_PLOTS_CLEANED.pdf ${studyfinal}/.
    
    #yes | cp ${studytmp}/*JB*.pdf ${studyfinal}/.
    echo -e
    echo -e "${study} on ${chipv} now complete."
    }
    
not_automatic() {
    chipv=$1
    study=$2
    args=($@)
    boundaries=$(echo ${args[@]:2})
    studytmp=${qctemp}/${study}
    studyfinal=${final_folder}/${study}
    echo $chipv
    echo $study
    echo $boundaries
    
    echo -e
    echo -e "Starting Step12 for $study on chip $chipv. Writing intermediate files here: ${studytmp}"
     bash ${wmscript}/Step12_HELLAC_QC_shared_v5.sh ${study} ${chipv} 
    
    
    echo -e
    echo "Now starting step34.  Writing final result files here: ${studyfinal}"
    bash ${wmscript}/Step34_HELLAC_QC_shared_v5.sh ${study} ${chipv} $(echo ${boundaries})
     
    echo -e
    echo "Now drawing separate PCA graph with these automated boundaries: $(cat ${studytmp}/${study}_PC_boundary.txt)"
    
    draw_boundaries ${study}  ${chipv} "auto_bound" $(cat ${studytmp}/${study}_PC_boundary.txt)
    
    echo -e
    echo "Now drawing separate PCA graph with these new, non-automated human-determined boundaries: $(echo ${boundaries})"
    draw_boundaries ${study}  ${chipv} "manual_bound" $(echo ${boundaries})
    
     
    echo -e
    echo "Now drawing scatter plots."
    draw_scatters ${study}  ${chipv}
    
    yes | cp ${study}_projpca_manual_bound.pdf ${studyfinal}/${study}_hapmap_rooted_PCA_PLOTS.pdf
    
    #rm -f ${studyfinal}/*JB*.pdf
    yes | cp ${studytmp}/${study}_hapmap_rooted_PCA_PLOTS_CLEANED.pdf ${studyfinal}/.
    yes | cp ${studytmp}/${study}_CEU_TSI_rooted_PCA_PLOTS_CLEANED.pdf ${studyfinal}/.
    yes | cp ${studytmp}/${study}_Unrooted_PCA_PLOTS_CLEANED.pdf ${studyfinal}/.
    
    
    echo -e
    echo -e "${study} on ${chipv} now complete."
    }

    
    
    
 ### First do the automatics
    
if [[ ${chip} == 24-1.0 ]]; then

    autolist=(${autolist_24v1[@]})
    
elif [[ ${chip} == 12-1.0 ]]; then
    
    autolist=(${autolist_12v1[@]})
#     echo skipping automatics

elif [[ ${chip} == 12-1.1 ]]; then
    
    autolist=(${autolist_12v11[@]})
#         echo skipping automatics
else 
    echo -e "Chip Version $chip does not exist"
    exit 1
    
fi
    
   for coh in ${autolist[@]};
 do
    completely_automatic ${chip} ${coh}
 done


 ### Now do the non-automatics
 
     
if [[ ${chip} == 24-1.0 ]]; then

    not_automatic 24-1.0 Cohort15 0.0015 0.0072 0.0455 0.054   
    
    not_automatic 24-1.0 Cohort16 0.0025 0.0115 0.0369 0.0571
    
    not_automatic 24-1.0 Cohort17 0.0013 0.0084 0.0432 0.0599
    
    not_automatic 24-1.0  Cohort18 0 0.0087 0.0433 0.0549
    
    
    
elif [[ ${chip} == 12-1.0 ]]; then
    
    not_automatic 12-1.0 Cohort19 0.0028 0.0084 0.0432 0.0566
    
    not_automatic 12-1.0 Cohort20 0.0015 0.0077 0.044 0.0554
    
    not_automatic 12-1.0 Cohort21 0.0025 0.0115 0.0369 0.054
#     echo skipping non automatics
    

elif [[ ${chip} == 12-1.1 ]]; then
    
    not_automatic 12-1.1 Cohort22 0.0018 0.0075 0.0434 0.0555
    
    not_automatic 12-1.1 Cohort23 0.0025 0.0115 0.0369 0.054
    
    not_automatic 12-1.1 Cohort24 6e-04 0.008 0.0432 0.0589
#         echo skipping non automatics

else
    
    echo -e "Chip Version Does Not Exist"
    exit 1
fi
    
