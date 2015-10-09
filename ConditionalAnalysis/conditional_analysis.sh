#!/bin/bash

locuszoom=/m/cphg-expr1/cphg-expr1/locuszoom/bin/locuszoom
#alias python=/home/apps/python2.7/bin/python2.7
export PATH=/home/apps/bin/:$PATH

usage()
{
cat << EOF
usage: $0 options

This script runs conditional analysis on a region, drawing a graph at each stage. It relies on proper installation of LocusZoom and PLINK.

OPTIONS:
   -h      Show this message
   -d      DATA -- path to PLINK file for study
   -r      REGION -- the REGION ID of a single region of interest.
   -l      REGION List -- list of regions of interest.
   -f      REGION FILE - tab-delimited table of region information in megabases (colnames: Chr, Start, End, Title, ID) (N.B. only regions in the region list will be drawn, unless no region list is provided.)
   -x      Extra arguments file to be passed to locuszoom. If file is not given, there will be no title to the graph. The word "REGION" should be given in the title argument, to be substituted by the Region Title later.
   -i      INTERRUPTED - indicates that the region should be analyzed using files already provided, because the previous analysis was interrupted. (MHC ANYONE?)
   -o      output parent directory. (DEFAULT= current directory)
   -m      MAXIMUM P-VALUE - the maximum significance theshold. (DEFAULT = 0.05)
   -p      PLOT ONLY - If you are SURE that all the necessary files are present, only draw graphs. **********THIS IS NOT CURRENTLY IMPLEMENTED*********
   -n      NO PLOT - Skip the plotting and just find the significant hits.
   -t      PREFIX - String to add to front of graph name -- in case graphs with different arguments are desired in the same folder.
EOF
}

clean=1
plinkfile=
extraargs=''
regionFile=
maxp=0.05
xDIR=$(pwd)
plot=1
prefix=''
onlyplot=0
while getopts “hd:r:l:f:x:io:m:npt:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         d)
             plinkfile=$OPTARG
             ;;
         r)
             regionSingle=$OPTARG
             ;;
         l)
             regionList=$OPTARG
             ;;
         f)
             regionFile=$OPTARG
             ;;             
         x)
             extraargs=$OPTARG
             ;;
         i)
             clean=0
             ;;
         o)
             xDIR=$OPTARG
             ;;
         m)
             maxp=$OPTARG
             ;;
         n)
             plot=0
             ;;
         p)
             onlyplot=1
             ;;
         t)
             prefix=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z "$plinkfile" ]] || [[ -z "$regionFile" ]];
then
     usage
     exit 1
fi


if [[ -z $regionList ]] && [[ -z $regionSingle ]];
then
    idcol=$(head -n1 ${regionFile}| sed 's/\r//g' | sed 's/\t/\n/g'| awk 'BEGIN {FS="\t"}; $0 ~/^ID$/ {print NR}')
    awk -v idcol=${idcol} 'BEGIN{FS="\t"} NR>1{print $idcol}' ${regionFile} >  regions.tmp
    IFS=$'\n' read -d '' -r -a regions < regions.tmp
    #rm -f regions.tmp
elif [[ ! -z ${regionSingle} ]];
then
    regions[0]=${regionSingle}
else
    IFS=$'\n' read -d '' -r -a regions < ${regionList}
fi


regionalAnalysis(){

    regionID=$1
    
    echo -e
    echo -e "Now Begininning Conditional Analysis on Region ${regionID}: $(date)"
    echo -e
    

    ld_caches=$(pwd)/ld_caches
    mkdir ${ld_caches}
    
    keepgoing=1
    ##############Retrieve Region Information######################


    #Use column names to determine columns of interest
    titlecol=$(head -n1 ${regionFile}| sed 's/\r//g' | sed 's/\t/\n/g'| awk 'BEGIN {FS="\t"}; $0 ~/^Title$/ {print NR}')
    chrcol=$(head -n1 ${regionFile}| sed 's/\r//g' | sed 's/\t/\n/g'| awk 'BEGIN {FS="\t"}; $0 ~/^Chr$/ {print NR}')
    startcol=$(head -n1 ${regionFile}| sed 's/\r//g' | sed 's/\t/\n/g'| awk 'BEGIN {FS="\t"}; $0 ~/^Start$/ {print NR}')
    endcol=$(head -n1 ${regionFile}| sed 's/\r//g' | sed 's/\t/\n/g'| awk 'BEGIN {FS="\t"}; $0 ~/^End$/ {print NR}')
    idcol=$(head -n1 ${regionFile}| sed 's/\r//g' | sed 's/\t/\n/g'| awk 'BEGIN {FS="\t"}; $0 ~/^ID$/ {print NR}')

    

    chr=$(awk -v col=${chrcol} -v idcol=${idcol} -v id=${regionID} '$idcol==id {print $col}' ${regionFile})
    start=$(awk -v col=${startcol} -v idcol=${idcol} -v id=${regionID} '$idcol==id {print $col}' ${regionFile})
    end=$(awk -v col=${endcol} -v idcol=${idcol} -v id=${regionID} '$idcol==id {print $col}' ${regionFile})
    regtitle=$(awk -v col=${titlecol} -v idcol=${idcol} -v id=${regionID} '$idcol==id {print $col}' ${regionFile})

    startfloat=$(echo "${start}*1000000"| bc -l)
    endfloat=$(echo "${end}*1000000"| bc -l)

    echo -e "Region Chr: ${chr}"
    echo -e "Region Start: ${startfloat%.*} "
    echo -e "Region End: ${endfloat%.*} "
    echo -e

#############Extract region SNPs#######################

    plink --noweb --bfile $plinkfile --make-bed --chr ${chr} --from-mb ${start} --to-mb ${end} --out ${regionID}/${regionID}

    # Prepare labeling file
    echo -e "snp\tstring\tcolor" > ${regionID}/denotemarkers.txt
    awk 'BEGIN { FS="\t"; OFS="\t"}; {print "chr"$1":"$4,$2,"purple"}' ${regionID}/${regionID}.bim >> ${regionID}/denotemarkers.txt


##############Make Decisions######################


    if [ $clean -eq 1 ]; then
        > ${regionID}/${regionID}.txt 
        > ${regionID}/${regionID}_conditionlist.txt
    elif [ ! -f ${regionID}/${regionID}.txt ]; then
        >&2 echo -e "Although clean is FALSE, there is no file ${regionID}/${regionID}.txt. \nclean will be set to TRUE."
        clean=1
        > ${regionID}/${regionID}.txt 
        > ${regionID}/${regionID}_conditionlist.txt
    else
        echo -e "Clean is set to FALSE! The process will begin on this loop: "$(cat ${regionID}/${regionID}.txt | wc -l)
    fi


#############Region Specific Pre Loop#######################

    while [[ $keepgoing -eq 1 ]]; do

        counter=$(cat ${regionID}/${regionID}.txt | wc -l)
        awk '{print $2}' ${regionID}/${regionID}.txt > ${regionID}/${regionID}_conditionlist.txt
        covararg=""
        if [[ -z ${plinkfile}.cov ]]; then
            covararg="--covar ${plinkfile}.cov"
        fi
        
        plink --noweb --bfile ${regionID}/${regionID} --logistic --no-sex --out ${regionID}/${regionID}_n${counter} --hide-covar --condition-list ${regionID}/${regionID}_conditionlist.txt --maf 0.005 --ci 0.95 --sex ${covararg}
        
        # Because there are certain cases where this file will need to be altered, but we still will want to preserve the integrity of the original output, we are going to set this filename to a variable.
        currentassoc=${regionID}/${regionID}_n${counter}.assoc.logistic

    # --covar ${plinkfile}.cov
 
#Find necessary column numbers
        pcol=$(head -n1 ${currentassoc} | sed 's/\r//g;s/ /\n/g' | sed '/^$/d'| awk '$0 ~/^P$/ {print NR}')
        snpcol=$(head -n1 ${currentassoc} | sed 's/\r//g;s/ /\n/g' | sed '/^$/d'| awk '$0 ~/^SNP$/ {print NR}')
        bpcol=$(head -n1 ${currentassoc} | sed 's/\r//g;s/ /\n/g' | sed '/^$/d'| awk '$0 ~/^BP$/ {print NR}')
        ccol=$(head -n1 ${currentassoc} | sed 's/\r//g;s/ /\n/g' | sed '/^$/d'| awk '$0 ~/^CHR$/ {print NR}')
        orcol=$(head -n1 ${currentassoc} | sed 's/\r//g;s/ /\n/g' | sed '/^$/d'| awk '$0 ~/^OR$/ {print NR}')
        statcol=$(head -n1 ${currentassoc} | sed 's/\r//g;s/ /\n/g' | sed '/^$/d'| awk '$0 ~/^STAT$/ {print NR}')
    
        
        lowp=$(awk -v pcol=${pcol} 'min=="" || $pcol < min {min=$pcol}; END{ print min}' ${currentassoc})
        
    
        
        if [[ ${lowp} == "NA" ]];
        then keepgoing=0
            echo -e
            echo -e "There are no more SNPs in the region with any p-value at all. Are you sure you have the right region boundaries and/or the correct maximum p-value threshold? You should rethink your decisions." 
            echo -e 
            continue
        fi
    
        lowpcount=$(awk -v pcol=${pcol} '{print $pcol}' ${currentassoc} | grep "^${lowp}$" | wc -l)
        
        if [[ ${lowpcount} -gt 1 ]]; then
        
            # if there are multiple SNPs with p-values equal to the lowest, we will use |STAT| (alternatively sqrt(STAT^2) ) to determine which one is the "top"
            awk -v pcol=${pcol} -v lowp=${lowp} '$pcol==lowp' ${currentassoc} \
            | awk -v statcol=${statcol} -v counter=${counter} 'max=="" || sqrt($statcol^2) > max {max=$statcol; maxline=$0}; END{ print maxline,counter}' | sed 's/ /\t/g'  >> ${regionID}/${regionID}.txt
        
        else
            awk -v pcol=${pcol} -v counter=${counter} -v lowp=${lowp} '$pcol==lowp {print $0,counter}' ${currentassoc} | sed 's/ /\t/g' >> ${regionID}/${regionID}.txt
        fi
        
        
        if [[ ${lowp} -eq 0 ]];
        then 
            nextlowp=$(awk -v pcol=${pcol} 'min=="" || $pcol < min && $pcol != 0 {min=$pcol}; END{ print min}' ${currentassoc})
            newlowp=$(echo $nextlowp | awk '{print $1 / 10}')
            echo -e
            echo -e "NOTE: PLINK has determined the lowest p-value to be 0. In order to draw, locuszoom will be told that the pvalue for any SNP with p=0 is: ${newlowp}" 
            echo -e 
            awk -v pcol=$pcol -v new=$newlowp '{sub(/^0$/, new, $pcol) }1;' ${currentassoc} > ${currentassoc}.adapted
            
            echo -e "IMPORTANT: the value has been changed in ${currentassoc}.adapted. The original file remains here: ${currentassoc}" 
            
            ### NOTE: The value of ${currentassoc} is being changed here
            
            currentassoc=${currentassoc}.adapted
            
            #yes | cp ${regionID}/${regionID}_n${counter}.assoc.logistic ${regionID}/${regionID}_n${counter}.assoc.logistic.original
            #yes | mv ${regionID}/${regionID}_n${counter}.assoc.logistic.adapted ${regionID}/${regionID}_n${counter}.assoc.logistic
            
        fi
        
        echo -e "SNP\tEXTREME\tP" > ${regionID}/${regionID}_n${counter}lz.txt
        
        awk -v snp=${snpcol} -v pval=${pcol} -v chr=${ccol} -v bp=${bpcol} \
        'BEGIN { OFS="\t"}; NR>1{print $snp,"chr"$chr":"$bp,$pval}' ${currentassoc} >> ${regionID}/${regionID}_n${counter}lz.txt
    
        moreargs=$(sed "s/REGION/${regtitle} Loop ${counter}/g" ${extraargs})
        echo ${moreargs}
      
        if [[ $plot -eq 1 ]];
        then
            curdir=$(pwd)
            cd ${regionID}
            echo ${endfloat%.*}
            $locuszoom \
                --metal ${regionID}_n${counter}lz.txt \
                --prefix ${prefix}${regionID}_n${counter} \
                --markercol "EXTREME" --pvalcol "P" \
                --plotonly  --delim tab --no-date \
                --chr ${chr} --start ${startfloat%.*} --end ${endfloat%.*}  \
                --cache ${ld_caches}/${regionID}.ldcache.db \
                ${moreargs} \
                #--plotonly 
                #--cache ${regionID}.ld.cache \
                #--pop EUR --build hg19 --source 1000G_Nov2010 \
                #--denote-markers ${regionID}/denotemarkers.txt
                #--gwas-cat whole-cat_significant-only \
                #--pop EUR --build hg19 --source #1000G_March2012 \
                #--add-refsnps ${extramarkers}
                #
            cd ${curdir}
        fi
        echo "HEREHEREHERE"
        echo $lowp
    #check if pvalue is lower than maximum pvalue
        keepgoing=$(echo "${lowp}<${maxp}"| sed -e 's/[eE]+*/*10^/'| bc -l)
        echo $keepgoing
    # cat lz_arguments.txt | sed 's/\"X$/\: $counter SNPs conditioned/g'

    done
    
    conlength=$(cat ${regionID}/${regionID}_conditionlist.txt | wc -l)
    echo $conlength
    if [[ "${conlength}" -gt 0 ]] && [[ "${plot}" -eq 1 ]]; then
        echo -e
        echo -e "There was at least one significant hit in this region! A plot will be drawn using all independent hits for reference!"
        cat ${regionID}/${regionID}_conditionlist.txt
        echo -e
        
        echo -e "snp\tstring\tcolor" > ${regionID}/denotemarkers_top.txt
        grep -f ${regionID}/${regionID}_conditionlist.txt ${regionID}/denotemarkers.txt >> ${regionID}/denotemarkers_top.txt
        addref=$(awk 'NR>2{print $1}' ${regionID}/denotemarkers_top.txt | tr '\n' ','| sed 's/,$//g')
        moreargs=$(sed "s/REGION/${regtitle}/g" ${extraargs})
        echo "HERE IS WHAT I THINK addref is :"$addref
        echo '"'"$addref"'"'
        
        cd ${regionID}
        
        $locuszoom \
            --metal ${regionID}_n0lz.txt\
            --prefix ${prefix}${regionID} \
            --markercol "EXTREME" --pvalcol "P" \
            --plotonly --delim tab --no-date \
            --cache ${ld_caches}/${regionID}.ldcache.db \
            --chr ${chr} --start ${startfloat%.*} --end ${endfloat%.*}\
            --denote-markers denotemarkers_top.txt \
            --add-refsnps $addref \
            ${moreargs}\
            #--cache ${ld_caches}/${regionID}.ld.cache \
            
            #--pop EUR --build hg19 --source 1000G_March2012 \
        cd ..
    fi

}

plotAll(){
    regionID=$1
    
    echo -e
    echo -e "Now PLOTTING Region ${regionID}: $(date)"
    echo -e

    ld_caches=ld_caches
    looptotal=$(cat ${regionID}/${regionID}.txt | wc -l)
    for j in ${looptotal}; do
    continue
    
    done
}




 for i in $(seq $(echo ${#regions[@]})); do
    index=$(($i-1))
    echo $index
    region=${regions[${index}]}
    echo $region
    mkdir ${region}

    if [[ $clean -eq 0 ]]; then
        regionalAnalysis ${region} | tee -a ${region}/${region}_condition.log
    else
        regionalAnalysis ${region} | tee ${region}/${region}_condition.log
    fi
  done

  
#Create a summary file of the top hits
echo -e $(head -n1 ${regions[0]}/${regions[0]}_n0.assoc.logistic | sed 's/\r//g;s/ /\n/g' | sed '/^$/d'| sed 's/\n/\t/g')"\tLoop\tREGION"  > ${xDIR}/TOPSNPS_PER_REGION.txt
  
  for regiondir in $(ls -d ${xDIR}/*.*/); do
    reggie=$(basename ${regiondir})
    cat ${reggie}/${reggie}.txt | sed 's/\r//g;s/ /\n/g' | sed '/^$/d'| sed 's/\n/\t/' | awk -v region=${reggie} 'FS="\t"{print $0,region}' >> ${xDIR}/TOPSNPS_PER_REGION.txt
    
   done
   