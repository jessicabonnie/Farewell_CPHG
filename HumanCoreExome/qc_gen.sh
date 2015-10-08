#!/bin/bash

#################################################################
# Usage: qc_gen.sh {nickname} {rawfile} {chip} {overall_title} {pre-specified remove list}
#################################################################

nickname=$1
inrawfile=$2
chip=$3
overall_title=$4

hPATH="/t121/jdrfdn/projects/JDRFDN/hapmap"
cPATH=/h2/t74/cphgdesk/share/cphg_OnengutLab/Jessica/Farewell/HumanCoreExome


#cPATH=$(pwd)/../scripts
#"/t121/jdrfdn/projects/JDRFDN/workfiles/Jessica/dataFreeze2/scripts"
#overall_title="T1DGC"


temp_folder=$(pwd)
xPATH="${temp_folder}/../QC"
logfile="${xPATH}/${nickname}.log"
rawfile="${nickname}raw"

KING="king"
PLINK="plink"
R="R"

outsamplefile="$xPATH/sampletoberemoved.txt"
outsnpfile="$xPATH/snptoberemoved.txt"
outsexfile="$xPATH/updatesex.txt"

if [ ! -d $xPATH ]; then
  mkdir $xPATH
fi
echo "FID IID REASON" > $outsamplefile
echo "SNP REASON" > $outsnpfile

nout=$(cat ${inrawfile}.fam | wc -l)
nout2=$(cat ${inrawfile}.bim | wc -l)
echo "There are $nout samples with $nout2 SNPs in the raw dataset $inrawfile"

if [ $# == 5 ]; then
  insamplefile=$5
  step=0
  nextstep=$((${step}+1))
  echo "Step $step: remove pre-specified samples"
  date
  errorfile="errorprespecified.txt"
  awk '{if(NF>2){print $1, $2, $3}else if(NF==2){print $1, $2, "ColumnError"}}' $insamplefile > $errorfile
  nout=$(cat $errorfile | wc -l)
  $PLINK --bfile $inrawfile --remove $errorfile --make-bed --out $rawfile --noweb > $logfile
  echo "# Samples to be removed due to prespecified error is $nout"
  cat $errorfile >> $outsamplefile
  echo "qc_gen.sh $1 $2 $3 $4 $5" > $logfile
elif [ $# == 4 ]; then 
  rawfile=$inrawfile
  echo "qc_gen.sh $1 $2 $3 $4" > $logfile
else
  echo "Usage: qc_gen.sh {nickname} {rawfile} {chip} {overall_title} {pre-specified remove list}"
  exit
fi

echo -e
step=1
nextstep=$((${step}+1))
echo "Step $step: remove very poor quality SNPs"
date
$KING -b ${rawfile}.bed --bySNP --prefix ${nickname}raw >> $logfile
awk '$11 < 0.8 && $2 != "Y"' ${nickname}rawbySNP.txt | awk '{print $1}' > allmissingSNP.txt
#awk '$8+$9==0 || $9+$10==0' ${nickname}rawbySNP.txt | awk '{print $1}' > monomorphicSNP.txt
#awk '$4=="I" || $4=="D"' ${nickname}rawbySNP.txt | awk '{print $1}' > indels.txt
nout=$(cat allmissingSNP.txt | wc -l)
echo "# SNPs to be removed due to severe missingness is $nout"
#nout=$(cat monomorphicSNP.txt | wc -l)
#echo "# SNPs to be removed due to monomorphism is $nout"

cat allmissingSNP.txt > snptoberemovedTMP.txt
#cat monomorphicSNP.txt >> snptoberemovedTMP.txt
#cat indels.txt >> snptoberemovedTMP.txt

$PLINK --bfile $rawfile --exclude snptoberemovedTMP.txt --make-bed --out ${nickname}${nextstep} --noweb >> $logfile
awk '{print $1, "CallRateLessThan80"}' allmissingSNP.txt >> $outsnpfile
#awk '{print $1, "DiscordantInDuplicates"}' poorSNPinDuplicate.txt >> $outsnpfile
#awk '{print $1, "Monomorphic"}' monomorphicSNP.txt >> $outsnpfile
#awk '{print $1, "INDELS"}' indels.txt >> $outsnpfile

echo -e
step=2
nextstep=$((${step}+1))
echo "Step $step: remove poor samples"
date
$KING -b ${nickname}${step}.bed --bysample --prefix ${nickname}${step} >> $logfile
awk 'NR > 1 && $7 > 0.05' ${nickname}${step}bySample.txt | awk '{print $1,$2}' > lessthan95.txt
awk 'NR > 1 && $7 > 0.01' ${nickname}${step}bySample.txt | awk '{print $1,$2}' > lessthan99.txt
$PLINK --bfile ${nickname}${step} --remove lessthan99.txt --make-bed --out ${nickname}${nextstep} --noweb >> $logfile 
awk '{print $1, $2, "MissingMoreThan5"}' lessthan95.txt >> $outsamplefile
nout=$(cat lessthan95.txt | wc -l)
echo "# samples to be removed due to call rate < 95% is $nout"
nout=$(cat lessthan99.txt | wc -l)
echo "# samples to be dropped from the following SNP-level QC step due to call rate < 99% is $nout"

#$R CMD BATCH -"${nickname}${step}bySample.txt" ${cPATH}/samplecallrate.R
$R CMD BATCH "--args ${nickname}${step}bySample.txt ${overall_title} ${chip}" ${cPATH}/samplecallrate_gen.R
pdffile="${xPATH}/samplecallrate.pdf"
ps2pdf samplecallrate.ps $pdffile
echo "Histogram of sample-level call rate is saved in $pdffile"

echo -e
step=3
nextstep=$((${step}+1))
echo "Step $step: remove poor quality SNPs in best quality samples"
date
snpqcfile="${nickname}${step}bySNP.txt"
$KING -b ${nickname}${step}.bed --cluster --bySNP --prefix ${nickname}${step} >> $logfile      

callrate=$(head -n1 $snpqcfile | sed 's/ /\n/g'| awk '$0 ~/^CallRate$/ {print NR}')
chr=$(head -n1 $snpqcfile| sed 's/ /\n/g'| awk '$0 ~/^Chr$/ {print NR}')
af=$(head -n1 $snpqcfile| sed 's/ /\n/g'| awk '$0 ~/^Freq_A$/ {print NR}')   

awk -v "cr=$callrate" -v "chr=$chr" '$cr < 0.95 && $chr!="Y"' $snpqcfile | awk '{print $1}' > lowcallrateSNP.txt
nout=$(cat lowcallrateSNP.txt | wc -l)
echo "# SNPs to be removed due to call rate < 95% is $nout"
awk 'NR>1 && $13>1 && ($15>0.01||$16>0.1)' $snpqcfile | awk '{print $1}' > poorSNPinDuplicate.txt
nout=$(cat poorSNPinDuplicate.txt | wc -l)
echo "# SNPs to be removed due to not concordant in duplicates (error rate>1% or adjusted error rate>10%) is $nout"
awk 'NR>1 && $19>1 && ($20>0.01||$21>0.1)' $snpqcfile | awk '{print $1}' > poorSNPinPO.txt
nout=$(cat poorSNPinPO.txt | wc -l)
echo "# SNPs to be removed due to MI in POs is $nout"
awk 'NR>1 && $24>1 && ($25>0.01||$26>0.1)' $snpqcfile | awk '{print $1}' > poorSNPinTrio.txt
nout=$(cat poorSNPinTrio.txt | wc -l)
echo "# SNPs to be removed due to MI in trios is $nout"
awk -v "cr=$callrate" -v "af=$af" '($af<=0.01||$af>=0.99)&&($cr<=0.99&&$cr>=0.95)' $snpqcfile | awk '{print $1}' > lowcallrateRareSNP.txt
awk -v "cr=$callrate" -v "af=$af" '($af>0.01&&$af<0.05) && ($cr<1-$6)' $snpqcfile | awk '{print $1}' >> lowcallrateRareSNP.txt
awk -v "cr=$callrate" -v "af=$af" '($af<0.99&&$af>0.95) && ($cr<$6)' $snpqcfile | awk '{print $1}' >> lowcallrateRareSNP.txt
nout=$(cat lowcallrateRareSNP.txt | wc -l)
echo "# SNPs to be removed due to low call rate in rare variants $nout"
#$R CMD BATCH -${snpqcfile} ${cPATH}/snpqc.R
$R CMD BATCH "--args ${snpqcfile} ${overall_title} ${chip}" ${cPATH}/snpqc_gen.R
pdffile="${xPATH}/snpqc.pdf"
ps2pdf snpqc.ps $pdffile
echo "Plots of SNP-level QC are saved in $pdffile"

awk '{print $1, "CallRateLessThan95"}' lowcallrateSNP.txt >> $outsnpfile
awk '{print $1, "ConcordantRateLessThan99"}' poorSNPinDuplicate.txt >> $outsnpfile
awk '{print $1, "MendelInconsistency"}' poorSNPinPO.txt >> $outsnpfile
awk '{print $1, "MendelInconsistency"}' poorSNPinTrio.txt >> $outsnpfile
awk '{print $1, "CallRateLessThanAlleleFrq"}' lowcallrateRareSNP.txt >> $outsnpfile

$PLINK --bfile $inrawfile --remove $outsamplefile --exclude $outsnpfile --make-bed --out ${nickname}${nextstep} --noweb >> $logfile

echo -e
step=4
nextstep=$((${step}+1))
echo "Step $step: remove gender errors"
date
$PLINK --bfile ${nickname}${step} --filter-males --out ${nickname}${step}male --noweb --make-bed >> $logfile
$KING -b ${nickname}${step}male.bed --bySNP --prefix ${nickname}${step}male >> $logfile
awk '$2=="Y" && $11 < 0.8' ${nickname}${step}malebySNP.txt | awk '{print $1}' > poorYSNPM1.txt

$PLINK --bfile ${nickname}${step} --filter-females --out ${nickname}${step}female --noweb --make-bed >> $logfile
$KING -b ${nickname}${step}female.bed --bySNP --prefix ${nickname}${step}female >> $logfile
awk '$2=="Y" && $11 > 0.1' ${nickname}${step}femalebySNP.txt | awk '{print $1}' > poorYSNPF1.txt

cat poorYSNPM1.txt >> snptoberemovedTMP.txt
cat poorYSNPF1.txt >> snptoberemovedTMP.txt
$PLINK --bfile ${nickname}${step} --exclude snptoberemovedTMP.txt --make-bed --out ${nickname}${nextstep} --noweb >> $logfile

$KING -b ${nickname}${nextstep}.bed --bysample --bySNP --prefix ${nickname}${nextstep} >> $logfile
sampleqcfile="${nickname}${nextstep}bySample.txt"
#Determine thresholds for use in identifying erroneaously male and female samples, originally this was hard coded to 700 ysnps. Also thresholds determinded for midsex and X0 errors
ysnps=$(awk '$2=="Y"' ${nickname}${nextstep}bySNP.txt | wc -l)
half_ysnps=$((${ysnps}/2))
onethird_ysnps=$((${ysnps}/3))
twothird_ysnps=$((${onethird_ysnps}*2))
xheterozygosity=0.1
echo "# Y-SNPs to be used as the cutoff value between males vs females is ${half_ysnps}"
awk -v hy=${half_ysnps} 'NR>1 && $5==2 && $11 > hy' $sampleqcfile | awk '{print $1,$2}' > femaleerror.txt
awk -v hy=${half_ysnps} 'NR>1 && $5==1 && $11 < hy' $sampleqcfile | awk '{print $1,$2}' > maleerror.txt
awk -v hy=${half_ysnps} '$5==0 && $11 > hy' $sampleqcfile | awk '{print $1, $2, 1}' > $outsexfile
awk -v hy=${half_ysnps} '$5==0 && $11 <= hy' $sampleqcfile | awk '{print $1, $2, 2}' >> $outsexfile

awk -v "y13=$onethird_ysnps" -v hy=${half_ysnps} 'NR>1 && $5==2 && $11>y13 && $11<hy' $sampleqcfile | awk '{print $1,$2}' > gendererror.txt
awk -v "y23=$twothird_ysnps" -v hy=${half_ysnps} 'NR>1 && $5==1 && $11>hy && $11<y23' $sampleqcfile | awk '{print $1,$2}' > gendererror.txt 
awk -v xh=$xheterozygosity -v hy=${half_ysnps} 'NR>1 && $5==2 && $10 < xh && $11 < hy' $sampleqcfile | awk '{print $1, $2}' >> gendererror.txt
awk -v xh=$xheterozygosity -v hy=${half_ysnps} 'NR>1 && $5==1 && $10 > xh && $11 > hy' $sampleqcfile | awk '{print $1, $2}' >> gendererror.txt
nout=$(cat femaleerror.txt | wc -l)
echo "# samples to be removed for being mislabeled as a female is $nout"
nout=$(cat maleerror.txt | wc -l)
echo "# samples to be removed for being mislabeled as a male is $nout"
nout=$(cat gendererror.txt | wc -l)
echo "# samples to be removed for other gender QC errors is $nout"

cat femaleerror.txt > sexerror.txt
cat maleerror.txt >> sexerror.txt
cat gendererror.txt >> sexerror.txt

$PLINK --bfile ${nickname}${step}male --remove sexerror.txt --update-sex $outsexfile --make-bed --out ${nickname}${step}maleB --noweb >> $logfile
$PLINK --bfile ${nickname}${step}female --remove sexerror.txt --make-bed --out ${nickname}${step}femaleB --noweb >> $logfile
$KING -b ${nickname}${step}maleB.bed --bySNP --prefix ${nickname}${step}maleB  >> $logfile
awk '$2=="Y" && $11 < 0.95' ${nickname}${step}maleBbySNP.txt | awk '{print $1}' > poorYSNPM2.txt
$KING -b ${nickname}${step}femaleB.bed --bySNP --prefix ${nickname}${step}femaleB >> $logfile
awk '$2=="Y" && $11 > 0.02' ${nickname}${step}femaleBbySNP.txt | awk '{print $1}' > poorYSNPF2.txt
nout=$(cat poorYSNPM2.txt | wc -l)
echo "# Y-chr SNPs to be removed due to call rate < 95% among males is $nout"
nout=$(cat poorYSNPF2.txt | wc -l)
echo "# Y-chr SNPs to be removed due to being present in females (>2%) is $nout" 

cat poorYSNPM2.txt > snptoberemovedTMP.txt
cat poorYSNPF2.txt >> snptoberemovedTMP.txt
#$R CMD BATCH -${sampleqcfile} -${half_ysnps} ${cPATH}/gender.R
$R CMD BATCH "--args ${sampleqcfile} ${half_ysnps} ${overall_title} ${chip}" ${cPATH}/gender_gen.R
pdffile="${xPATH}/gender.pdf"
ps2pdf gender.ps $pdffile
echo "Gender checking plot is saved in $pdffile"

$PLINK --bfile ${nickname}${step} --remove sexerror.txt --exclude snptoberemovedTMP.txt --make-bed --out ${nickname}${nextstep} --noweb >> $logfile
awk '{print $1, $2, "MislabeledAsFemale"}' femaleerror.txt >> $outsamplefile
awk '{print $1, $2, "MislabeledAsMale"}' maleerror.txt >> $outsamplefile
awk '{print $1, $2, "GenderQC"}' gendererror.txt >> $outsamplefile
awk '{print $1, "YCallRateLessThan95"}' poorYSNPM2.txt >> $outsnpfile
awk '{print $1, "YInFemales"}' poorYSNPF2.txt >> $outsnpfile

echo -e
step=5
nextstep=$((${step}+1))
echo "Step $step: population structure analysis"
date
awk '{print $1,$2,$3,$4,$5,2}' ${nickname}${step}.fam > ${nickname}${step}a.fam
cp ${nickname}${step}.bed ${nickname}${step}a.bed
cp ${nickname}${step}.bim ${nickname}${step}a.bim

echo "Step ${step}a: projecting to all HapMap samples"

#$KING -b ${nickname}${step}a,$cPATH/../hapmap/hapmap --pca --projection --prefix ${nickname}${step}a >> $logfile
$KING -b ${nickname}${step}a,$hPATH/hapmap --pca --projection --prefix ${nickname}${step}a >> $logfile
awk '$6==2' ${nickname}${step}apc.ped | awk '$8 > 0.045' | awk '{print $1, $2}' > keeplist.txt
awk '$6==2' ${nickname}${step}apc.ped | awk '$8 <= 0.045' | awk '{print $1, $2}' > nonEUR.txt    
nout=$(cat nonEUR.txt | wc -l)
echo "# non-European samples is $nout"
$PLINK --bfile ${nickname}${step}a --keep keeplist.txt --make-bed --out ${nickname}${step}b --noweb >> $logfile
#$R CMD BATCH -"${nickname}${step}apc.ped" $cPATH/projpca.R
$R CMD BATCH "--args ${nickname}${step}apc.ped  ${overall_title} ${chip}" $cPATH/projpca_gen.R

pdffile="${xPATH}/projpca.pdf"
ps2pdf projpca.ps $pdffile
echo "Population structure projected to HapMap is saved in $pdffile"

echo -e
echo "Step ${step}b: projecting to European HapMap (CEU+TSI) samples"
$KING -b ${nickname}${step}b,$hPATH/EURhapmap --pca --projection --prefix ${nickname}${step}b >> $logfile
posCount=$(awk '$6==2 && $7 > 0' ${nickname}${step}bpc.ped | wc -l)
negCount=$(awk '$6==2 && $7 <= 0' ${nickname}${step}bpc.ped | wc -l)
diffCount=$(($posCount-$negCount))
awk '$6==2' ${nickname}${step}bpc.ped | awk -v diff=$diffCount '{if((diff>0 && $7>0)||(diff<=0 && $7<=0)){print $1, $2}}' > keeplist.txt
nout=$(cat keeplist.txt | wc -l)
echo "# samples to be used for HWE checking is $nout"
$PLINK --bfile ${nickname}${step} --keep keeplist.txt --make-bed --out ${nickname}${step}c --noweb >> $logfile
$KING -b ${nickname}${step}c.bed --unrelated --prefix ${nickname}${step}c >> $logfile
$PLINK --bfile ${nickname}${step} --keep ${nickname}${step}cunrelated.txt --hwe 0.000001 --hardy --noweb --out ${nickname}${step}c >> $logfile

cat ${nickname}${step}.bim > NE1.bim
echo $(head -n1 ${nickname}${step}c.hwe) "POS"  > NE2pos.hwe
awk 'BEGIN{while((getline<"NE1.bim")>0)l[$2]=$0}$2 in l{print $0"\t"l[$2]}' ${nickname}${step}c.hwe | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13}' >> NE2pos.hwe
grep "ALL" NE2pos.hwe | awk '$9 < 1E-6'| awk '{if($1==6 && $NF>24482000 && $NF<33859000){if($9<1E-20){print $2}}else{print $2}}'>hweSNP_20.txt
nout=$(cat hweSNP_20.txt | wc -l)
echo "# SNPs to be removed due to HWE P<1E-6 (or <1E-20 for MHC) is $nout"
awk '{print $1, "HWDInEUR"}' hweSNP_20.txt >> $outsnpfile
$PLINK --bfile ${nickname}${step} --exclude hweSNP_20.txt --make-bed --out ${nickname}${nextstep} --noweb >> $logfile
#$R CMD BATCH -"${nickname}${step}bpc.ped" $cPATH/EURprojpca.R
$R CMD BATCH "--args ${nickname}${step}bpc.ped ${overall_title} ${chip}" $cPATH/EURprojpca_gen.R

pdffile="${xPATH}/EURprojpca.pdf"
ps2pdf EURprojpca.ps $pdffile
echo "Population structure projected to CEU+TSI is saved in $pdffile"

echo -e
time=$(date)
echo "QC is done at $time"

finalstep=6
nout=$(cat $outsamplefile | wc -l)
echo "List of samples to be removed (N=${nout}) is saved in $outsamplefile"
nout=$(cat $outsnpfile | wc -l)
echo "List of SNPs to be removed (N=${nout}) is saved in $outsnpfile"
nout=$(cat $outsexfile | wc -l)
if [ $nout > 0 ]; then
  echo "List of samples to have sex updated (N=${nout}) is saved in $outsexfile"
finalfile=${nickname}${finalstep}
nout=$(cat ${finalfile}.fam | wc -l)
nout2=$(cat ${finalfile}.bim | wc -l)
echo "There are $nout samples with $nout2 SNPs in the final dataset $finalfile"  
fi


