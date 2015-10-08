#!/bin/bash -x
#!/usr/bin/bash
##############################################################################
#
# Title: Script commands to run zCall on JDRF Samples
#
# Desc:  Script to Process GenomeStudio SNP report and run 
#        zCall for JDRF samples genotyped at University of Virginia
#
# Note: See "Using zCall on raw Genotype Files " SOP
#
# Date: Feb, 2015
# By:   jbonnie

##############################################################################

zcall_fulldata=$1
chip_version=$2
sampletoberemoved=$3
rawplink=$4


#I assume that we are in the tempfolder already, final results will be put in folder above
workdir=$(pwd)
xPATH="${workdir}/../zCalled"



# if [ $# == 5 ]; then
#   cPATH=$5
# else
# cPATH="/t121/jdrfdn/projects/JDRFDN/workfiles/Jessica/dataFreeze2/scripts"
# fi

echo -e
echo "Here are the arguments passed to the script: "
echo "____________________________________________"
echo -e
echo "Raw GS Report Location: " $zcall_fulldata
echo "HumanCoreExome Chip Version (12v1, 12v1.1, or 24v1): " $chip_version
echo "List of Samples to be removed because of prior analysis: " $sampletoberemoved
echo "Raw Plink Files used in prior analysis: " $rawplink
# echo "Script folder being used: " $cPATH

#zCall scripts are here:
zPATH=/t121/jdrfdn/projects/JDRFDN/apps/zCall/Version3_GenomeStudio
Rscript='/h4/t1/apps/stat/R-3.0.3/bin/Rscript'

echo -e
echo "Here are hard coded paths to programs: "
echo "____________________________________________"
echo -e
echo "Path to Rscript: " $zPATH
echo "Path to zCall: " $Rscript


#I assume that we are in the tempfolder already

#basicdaughter=${workdir}/basicThresholds
basicdaughter=${workdir}
echo -e
echo "Top Level Directory is : ${workdir}"
echo -e


# Find the Name of the Report File

# reportname=$(basename "$zcall_fulldata")
# extension="${reportname##*.}"
# reportname="${reportname%.*}"


reportname=zcallreport
report_qced=${workdir}/${reportname}_QCed.txt


echo -e
echo -e "\n------------------"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo -e
echo "I. Create Files Needed"
echo -e
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo -e "------------------\n"
echo -e


echo -e
echo -e "\n------------------"
echo "Step I.1. QC Raw Report for Feeding Into Zcall"
echo -e "------------------\n"
echo -e

  echo -e
  echo -e "\n------------------"
  echo "Step I.1.A. Removing List of Known Bad Samples from Report"
  echo -e "------------------\n"
  echo -e
  # Step 1a. Removing List of Samples from Report

    #Edit the sample list to be used in during zcall
    sampletoberemovedZ=${workdir}/sampletoberemovedZ.txt
    awk '{print $2}' ${sampletoberemoved} > ${sampletoberemovedZ}

    #Here is where we will write the output when we remove the samples
    echo ${workdir}/${reportname}_QCed1a.txt

    #Use zcall's additional scripts version of dropSamplesFromReport (THE OTHER ONE DOESN'T WORK) to remove the given samples
    python  ${zPATH}/../additionalScripts/dropSamplesFromReport_FasterVersion.py ${zcall_fulldata} ${sampletoberemovedZ} > ${workdir}/${reportname}_QCed1a.txt
    echo "Complete"
    #DID NOT REALIZE HOW GINORMOUS A FILE THAT IS TO KEEP AROUND ... perhaps zip it up or remove it after usage?

  echo -e
  echo -e "\n------------------"
  echo "Step I.1.B. Applying Subject Call Rate Threshold Filter (98%)"
  echo -e "------------------\n"
  echo -e
    # Step 1b. Applying Subject Call Rate Threshold Filter (98%)

    python  ${zPATH}/qcReport.py -R  ${workdir}/${reportname}_QCed1a.txt -C 0.98 > ${report_qced}
    echo "Complete"

  echo -e
  echo -e "\n------------------"
  echo "Step I.1.C. Counting number of subjects - Original and after filtering"
  echo -e "------------------\n"
  echo -e
  # 1c. Counting number of subjects - Original and after filtering
  # Original File (divide number by 3 for subject count)
    echo "Subject Counts:"
    echo "Original Report: "$(head -1 ${zcall_fulldata} | tr '\t ' '\n' | awk 'NR > 3' | wc)

    # File with subjects removed from list (divide number by 3 for subject count)
    echo "Minus Predetermined Bad Samples: " $(head -1 ${workdir}/${reportname}_QCed1a.txt | tr '\t ' '\n' | awk 'NR > 3' | wc)

  # Subject Call rate filtered File (divide number by 3 for subject count)
    echo "After Subject Callrate Filter (98%):" $(head -1 ${report_qced} | tr '\t ' '\n' | awk 'NR > 3' | wc)
    
  echo -e
  echo -e "\n------------------"
  echo "Step I.1.D. Zip Up Intermediate QC report"
  echo -e "------------------\n"
  echo -e
    gzip -1qf ${workdir}/${reportname}_QCed1a.txt
    echo "Complete"
    

  echo -e
  echo -e "\n------------------"
  echo "Step I.2. Calculate Mean & Standard Deviations for each Homozygote Cluster"
  echo -e "------------------\n"
  echo -e
    meanfile=${workdir}/${reportname}_MEAN.SD.txt
    python  ${zPATH}/findMeanSD.py -R ${report_qced} > ${meanfile}
    echo "Complete"
    echo "Line Count of ${meanfile}: " $(wc -l ${meanfile})

  echo -e
  echo -e "\n------------------"
  echo "I.3.: Create Rare SNP List using un-zCalled raw plink files (those used by Wei-Min in Initial QC)"
  echo -e "------------------\n"
  echo -e

    echo -e
    echo -e "\n------------------"
    echo "Step I.3.A: Calculate Minor Allele Frequency and Missingness in raw samples not in the samplestoberemoved list"
    echo -e "------------------\n"
    echo -e
      
      # I am not currently using this king approach, but here is the beginning of implementing it
#       plink --noweb --bfile ${rawplink} --remove ${sampletoberemoved} --out ${workdir}/raw_cleaned --make-bed
    
#       king -b raw_cleaned.bed --bySNP --prefix ${workdir}/raw_genotype_cleaned
      
#       awk '($6 >= 0.99 || $6 <=0.01) && $11 >= 0.99' ${workdir}/raw_genotype_cleanedbySNP.txt > ${workdir}/rareSNP_details_cleaned.txt
#       awk 'NR>1{print $1}' ${workdir}/rareSNP_details_cleaned.txt > ${workdir}/rareSNPs_cleaned.txt
      
      
      plink --noweb --bfile ${rawplink} --freq --out ${workdir}/rawgeno --remove ${sampletoberemoved}
      plink --noweb --bfile ${rawplink} --missing --out ${workdir}/rawgeno_miss --remove ${sampletoberemoved}
      
      LANG=en_EN join -1 2 -2 2 <(LANG=en_EN sort -k2,2 rawgeno.frq) <(LANG=en_EN sort -k2,2 rawgeno_miss.lmiss) > ${workdir}/raw_geno_snpdetails.txt

      
    echo -e
    echo -e "\n------------------"
    echo "I.3.B: Create List of Rare (MAF<1%) SNPs with Callrates of 99% or Greater"
    echo -e "------------------\n"
    echo -e
      
      awk '$5 <=.01 && $10 <=.01' ${workdir}/raw_geno_snpdetails.txt > ${workdir}/rare99.txt
      echo $(cat ${workdir}/rare99.txt | wc -l) "Rare SNPs"
      


  echo -e
  echo -e "\n------------------"
  echo "I.4.: Create Allele and Flipping Update Files from Manifest"
  echo -e "------------------\n"
  echo -e
  
  
    echo -e
    echo -e "\n------------------"
    echo "Step I.4.A. Identify correct manifest files based on chip version"
    echo -e "------------------\n"
    echo -e

      manifest_loc=XXXXXXXXXXXX
      snp_table=XXXXXXXXXXXX
      if [ "${chip_version}" = "12v1" ] || [ "${chip_version}" = "12v1.0" ]; then 
        manifest_loc=/t121/jdrfdn/projects/JDRFDN/HumanCoreExome_12v1_JDRF_Upload/HumanCoreExome-12-v1-0-D.csv
        snp_table=XXXXXXXXXXXXXXX
      
      elif [ "${chip_version}" = "12v1.1" ]; then 
        manifest_loc=/t121/jdrfdn/projects/JDRFDN/HumanCoreExome_12v1.1_JDRF_Upload/HumanCoreExome-12-v1-1-C.csv
        snp_table=/t121/jdrfdn/projects/JDRFDN/HumanCoreExome_12v1.1_JDRF_Upload/12082014-JDRFCoreExome1.1_SNP_Table.txt
 
      elif [ "${chip_version}" = "24v1" ] || [ "${chip_version}" = "24v1.0" ]; then 
        manifest_loc=/t121/jdrfdn/projects/JDRFDN/HumanCoreExome_24v1_JDRFDN_Freeze2/HumanCoreExome-24v1-0_A.csv
        snp_table=XXXXXXXXXXXXXXX
      fi
      echo "Manifest Loc set to ${manifest_loc}"

    echo -e
    echo -e "\n------------------"
    echo "Step I.4.B. Prepare File to Update Illumina Alleles (A/B) to ATCG Later"
    echo -e "------------------\n"
    echo -e
      illumina_allele_update=${workdir}/${chip_version}_illumina_allele_update.txt
      cat ${manifest_loc} | awk 'BEGIN {FS=","};NR > 8 {print $2, $4}' | sed -e 's/\[//g' -e 's/\]//g' -e 's/\//\t/g' | awk '{print $1, "\t", "A", "\t", "B", "\t", $2, "\t", $3}' > ${illumina_allele_update}

    echo -e
    echo -e "\n------------------"
    echo "Step I.4.C. Prepare Master File to Aid in Flipping to Match Top Allele Format Later"
    echo -e "------------------\n"
    echo -e
      illumina_topflip_update=${workdir}/${chip_version}_illumina_topflip_update.txt

       cat ${manifest_loc} | awk 'BEGIN {FS=","};NR>7 && NF > 4{print $2, $3, $21}' | sed 's/^Name/SNP/' >  ${illumina_topflip_update}

    echo -e
    echo -e "\n------------------"
    echo "Step I.5. EXCLUDE RARE SNPs from Raw PLINK Files"
    echo -e "------------------\n"
    echo -e

      plink --noweb --exclude ${workdir}/rare99.txt --bfile ${rawplink} --make-bed --out ${workdir}/raw_excluderare99 --allow-no-sex


    echo -e
    echo -e "\n------------------"
    echo "I.6.: Create Update Files from Raw PLINK File"
    echo -e "------------------\n"
    echo -e
      awk '{print $2,$2,$1,$2}' ${rawplink}.fam > updateids.txt
      awk '{print $1,$2,$5}' ${rawplink}.fam > updatesex.txt
      awk '{print $1,$2,$3,$4}' ${rawplink}.fam > updateparents.txt
      awk '{print $1,$2,$6}' ${rawplink}.fam > updatestatus.txt



echo -e
echo -e "\n------------------"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo -e
echo "II. Use Thresholds That are Calculated Using zCall's built in Thresholds for SNPs"
echo -e
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo -e "------------------\n"
echo -e


#basicdaughter=${workdir}/basicThresholds
basicdaughter=${workdir}

 ### basicdaughter is a remnant of the totally ridiculous folder structure... we are now making all the intermediate files write to the same work/temp folder






    echo -e
    echo -e "\n------------------"
    echo "Step II.1. Recalculate BETAs using ${meanfile} file"
    echo -e "------------------\n"
    echo -e
      $Rscript ${zPATH}/findBetas.r ${meanfile}  ${workdir}/BETAS.txt 1
      echo "Complete"

    echo -e
    echo -e "\n------------------"
    echo "Step II.2. Find the thresholds used for calling No Calls"
    echo -e "------------------\n"
    echo -e

      python ${zPATH}/findThresholds.py -R ${report_qced} -B ${workdir}/BETAS.txt -Z 7 -I 0.2 > ${workdir}/thresholds.7X
      echo "Complete"
      
    echo -e
    echo -e "\n------------------"
    echo "Step II.3. Recall all NoCalls using these Thresholds"
    echo -e "------------------\n"
    echo -e

        python ${zPATH}/zCall.py -R ${zcall_fulldata} -T ${workdir}/thresholds.7X -O ${workdir}/zcall_preupdate
        
        echo "Complete"

    echo -e
    echo -e "\n------------------"
    echo "Step II.4. Convert zCalled file to PLINK Bim/Bed Format, Keeping only RARE SNPs"
    echo -e "------------------\n"
    echo -e
 ### NOTE - one could really extract the raw snps at the step and speed up the whole process... but in the interest of not fixing what ain't broke, I'm leaving it as is.
        plink --noweb --tfile ${workdir}/zcall_preupdate  --make-bed --out  ${workdir}/zcall_preupdate_nopheno --allow-no-sex --extract ${workdir}/rare99.txt
#         plink --noweb --tfile ${workdir}/zcall_preupdate  --make-bed --out  ${workdir}/zcall_preupdate_nopheno --allow-no-sex 
     
    echo -e
    echo -e "\n------------------"
    echo "Step II.5. Update Phenotype Info: Sex/IDs/Parents/Status"
    echo -e "------------------\n"
    echo -e
        
        plink --noweb --bfile ${workdir}/zcall_preupdate_nopheno --update-ids ${workdir}/updateids.txt --make-bed --out ${workdir}/zcall_preupdateA 
        plink --noweb --bfile ${workdir}/zcall_preupdateA --update-sex ${workdir}/updatesex.txt --make-bed --out  ${workdir}/zcall_preupdateB 
        plink --noweb --bfile ${workdir}/zcall_preupdateB --update-parents ${workdir}/updateparents.txt --make-bed --out  ${workdir}/zcall_preupdateC 
        plink --noweb --bfile ${workdir}/zcall_preupdateC --pheno ${workdir}/updatestatus.txt --make-bed --out  ${workdir}/zcall_preupdate 
         
    echo -e
    echo -e "\n------------------"
    echo "Step II.6. Update Alleles"
    echo -e "------------------\n"
    echo -e
    
        plink --noweb --bfile ${workdir}/zcall_preupdate  --make-bed --out  ${workdir}/zcall_ALLELES --update-alleles ${illumina_allele_update} --allow-no-sex

    echo -e
    echo -e "\n------------------"
    echo "Step II.7. Flip Certain Rare SNPs to return to TOP ALLELE format"
    echo -e "------------------\n"
    echo -e
      
      
        echo -e
        echo -e "\n------------------"
        echo "Step II.7.A. Run frequency test on zCalled data -- Remember these are RARE SNPs only."
        echo "Step II.7.B. Check Whether A2 Values of Rare SNPs Match between un-zCalled and zCalled frequency files"
        echo -e "------------------\n"
        echo -e
        
        plink --noweb --bfile ${workdir}/zcall_ALLELES --out ${workdir}/zcall_ALLELES_rare --freq --remove ${sampletoberemoved}
        #plink --noweb --bfile ${workdir}/zcall_ALLELES --out ${workdir}/zcall_ALLELES_rare --freq --remove ${sampletoberemoved} --extract ${workdir}/rare99.txt
        
        
        LANG=en_EN join -1 2 -2 2 <( awk '{print $1,$2,$3,$4,$5}' ${workdir}/zcall_ALLELES_rare.frq | LANG=en_EN sort -k2 ) <(awk '{print $1,$2,$3,$4,$5}' ${workdir}/rawgeno.frq | LANG=en_EN sort -k2)  > ${workdir}/zcall_rare_FreqCompare.tmp
        grep -w SNP ${workdir}/zcall_rare_FreqCompare.tmp | sed 's/SNP CHR A1 A2 MAF CHR A1 A2 MAF/SNP CHR_Z A1_Z A2_Z MAF_Z CHR_GS A1_GS A2_GS MAF_GS/' > ${workdir}/zcall_rare_FreqCompare.txt
        grep -v -w SNP ${workdir}/zcall_rare_FreqCompare.tmp >> ${workdir}/zcall_rare_FreqCompare.txt
        
        
        LANG=en_EN join -1 1 -2 1 <( LANG=en_EN sort -k1 ${workdir}/zcall_rare_FreqCompare.txt) <(LANG=en_EN sort -k1 ${illumina_topflip_update})  > ${workdir}/zcall_rare_FreqCompare_withStrand.tmp
        
        grep -w SNP zcall_rare_FreqCompare_withStrand.tmp > ${workdir}/zcall_rare_FreqCompare_withStrand.txt
        grep -w -v SNP zcall_rare_FreqCompare_withStrand.tmp >> ${workdir}/zcall_rare_FreqCompare_withStrand.txt

        
        awk '$4 != $8' ${workdir}/zcall_rare_FreqCompare_withStrand.txt > ${workdir}/zcall_SNPs_A2_NOTMATCHED.txt 
        
        
        echo -e
        echo -e "\n------------------"
        echo "Step II.7.C. Flip ONLY those Rare SNPS whose A2 values DO NOT MATCH value in un-zCalled TOP ALLELE Frequency File"
        echo -e "------------------\n"
        echo -e
        
#         LANG=en_EN join -1 2 -2 2 <( awk '{print $1,$2,$3,$4,$5}' ${workdir}/zcall_ALLELES_rare.frq | LANG=en_EN sort -k2 ) <(awk '{print $1,$2,$3,$4,$5}' ${workdir}/rawgeno.frq | LANG=en_EN sort -k2) | awk '$4 != $8' > ${workdir}/zcall_SNPstoFlip.txt
        
#         grep -w -f <( awk '{print $1}' ${workdir}/zcall_SNPstoFlip.txt) ${illumina_topflip_update} > ${workdir}/rare_topflip_update.txt
        

#         plink --bfile ${workdir}/zcall_ALLELES --noweb --flip ${workdir}/zcall_SNPs_A2_NOTMATCHED.txt --make-bed --out ${workdir}/zcall_TOP --allow-no-sex
        
        plink --bfile ${workdir}/zcall_ALLELES --noweb --flip ${workdir}/zcall_SNPs_A2_NOTMATCHED.txt --make-bed --out ${workdir}/zcall_rare99

#     echo -e
#     echo -e "\n------------------"
#     echo "Step II.8. EXTRACT RARE lists from zCalled PLINK FILES"
#     echo -e "------------------\n"
#     echo -e

#       plink --bfile ${workdir}/zcall_TOP --noweb --extract ${workdir}/rare99.txt --make-bed --out ${workdir}/zcall_rare99 --allow-no-sex


    echo -e
    echo -e "\n------------------"
    echo "Step II.8. Merge zCalled Rare SNPs with Subset of unZCalled raw (from which rare SNPs have been EXCLUDED)"
    echo -e "------------------\n"
    echo -e

      plink --noweb --bfile ${workdir}/zcall_rare99 --bmerge ${workdir}/raw_excluderare99.bed ${workdir}/raw_excluderare99.bim ${workdir}/raw_excluderare99.fam --out ${workdir}/merge99 --make-bed --allow-no-sex

      
    echo -e
    echo -e "\n------------------"
    echo "Step II.7. Create Results from EXTRACTING RARES from zCall and EXCLUDING RARES from raw"
    echo -e "------------------\n"
    echo -e
      mkdir ${xPATH}

            for plinky in $(ls ${workdir}/merge99.*)
            do 
                cp ${plinky} ${xPATH}/$(basename ${plinky/merge99/raw_zcall})
            done


    echo -e
    echo -e "\n------------------"
    echo "Step III. Run Threshold and Calibration Checking"
    echo -e "------------------\n"
    echo -e
#       bash $cPATH/zCalibrationLoop.sh ${report_qced} ${workdir}/BETAS.txt $zPATH

