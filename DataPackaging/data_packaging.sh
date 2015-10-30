#!/bin/bash

#tracking_sheet=$1
#chip=$2
#study=$3
#superbase=$4
#geno_loc=$5
#plink_base=$6
#release_date=$7


function usage () { 
  cat << EOF

    $0 -t TRACKING_SHEET -C CHIP -s STUDY -b BASE -g GS_EXPORT -p PLINK -d RELEASE_DATE -m MANIFEST -c CLUSTER[-h | --help]

    OPTIONS:
      -h    This help text
      -t    TRACKING_SHEET, path to .csv tracking sheet
      -C    CHIP (I or E1.0 or E1.1), referring to immunochip or exomechip (versions 1.0 or 1.1)
      -s    STUDY
      -b    BASE, location where packaged folder should be written
      -g    GS_EXPORT, path to raw Genome Studio export file
      -d    RELEASE_DATE, given in YearMonthDay format (e.g. 20130431)
      -p    PLINK, path to PLINK root
      -m    MANIFEST, path to Illumnia SNP manifest
      -c    CLUSTER, path to Illumnia cluster file



EOF
}


while getopts t:C:s:b:g:p:d:m:c:h option
  do
    case "${option}" in
       t) tracking_sheet=${OPTARG};;
       C) chip=${OPTARG};;
       s) study=${OPTARG};;
       b) superbase=$OPTARG;;
       g) geno_loc=$OPTARG;;
       p) plink_base=$OPTARG;;
       d) release_date=$OPTARG;;
       m) manifest_loc=$OPTARG;;
       c) cluster_loc=$OPTARG;;
       \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
       h) usage; exit;;
    esac
done




line_count=$(cat ${tracking_sheet} | wc -l | cat)
track_count=$((line_count-9))
#today=$(date +"%Y%m%d")

## Hard Code the idat folder and manifest locations based on the chip argument
if [ "${chip}" == "I" ]; then
	chipname='IC'
	idat_loc='/h3/t5/cphgcore/ImmunoChip_Project/Immunochip_Image_Data'
	chiplong='Immunochip'
	snp_count=196,524

elif [ "${chip}" == "E1.0" ]; then
echo I know to use Human Exome Chip
	chipname='HCEchipV1.0'
	idat_loc='/h3/t5/cphgcore/JDRF_DN_Project/CoreExome_Image_Data'
	manifest_loc='/h3/t5/cphgcore/JDRF_DN_Project/Jessica_Manifest_files/HumanCoreExome-12-v1-0-D.csv'
	chiplong='HumanCoreExomeV1.0'
	snp_count=538,448
	
elif [ "${chip}" == "E1.1" ]; then
echo I know to use Human Exome Chip
	chipname='HCEchipV1.1'
	idat_loc='/h3/t5/cphgcore/JDRF_DN_Project/CoreExome_Image_Data'
	manifest_loc='/h3/t5/cphgcore/JDRF_DN_Project/Jessica_Manifest_files/HumanCoreExome-12-v1-1-C.csv'
	chiplong='HumanCoreExomeV1.1'
	snp_count=542,585
fi


## Name and make the directories
basename=${chipname}_${study}_${track_count}_UVA_${release_date}
basefolder=${superbase}/${basename}
idat_folder=${basefolder}/iDat
plink_folder=${basefolder}/PLINK_GenotypeFiles
illumina_folder=${basefolder}/IlluminaReferenceFiles
genotype_folder=${basefolder}/RawGenomeStudioExport


mkdir -p ${idat_folder} ${plink_folder} ${illumina_folder} ${genotype_folder}

#assign all the other filenames to be used in this script

echo $basefolder
path_list_loc=${idat_folder}/${basename}_files.list
geno_id_loc=${genotype_folder}/${basename}_genoids.list
md5_loc=${idat_folder}/${basename}.md5
tar_loc=${idat_folder}/${basename}.tar
tar_md5_loc=${idat_folder}/${basename}.tar.md5
zip_loc=${idat_folder}/${basename}.tar.gz
test_list=${idat_folder}/test_list.txt
test_track=${idat_folder}/test_track.csv
#plink_md5=${plink_folder}/plink.md5
plink_exts=(.bed .bim .fam .cov)
plink_fam=${plink_base}.fam
plink_bim=${plink_base}.bim
readme=${basefolder}/README
echo ${readme}

function move_and_check(){
	original_file=$1
	new_folder=$2
	filename=$(basename ${original_file})
	first_md5=$( cat ${original_file}| md5sum)
	cp ${original_file} ${new_folder}
	new_loc=${new_folder}/${filename}
	second_md5=$( cat ${new_loc} | md5sum)
	[ "$first_md5" = "$second_md5" ] && echo  "${filename} Successfully Copied" || echo "${filename}'s md5sum after copy has changed! Please Check"

	}


#Transfer a copy of the tracking sheet into the basefolder and check its md5sum
move_and_check ${tracking_sheet} ${basefolder}
#track_md5=$(cat ${tracking_sheet} | md5sum)
#cp ${tracking_sheet} ${basefolder}
#check_md5=$(cat ${basefolder}/*.csv | md5sum)
#[ "$track_md5" = "$check_md5" ] && echo  "Tracking Sheet Successfully Copied" || echo "Tracking Sheet's md5sum after copy has changed! Please Check"


#### Transfer Illumina Files and check their md5sums

for illumina in ${manifest_loc} ${cluster_loc}
  do
    move_and_check ${illumina} ${illumina_folder}
  done





# Time to deal with the raw genome studio export
#head -n1 $geno_loc | awk 
if [ "${geno_loc}" == "jkdlg" ]; then
IFS=$'\t' read -ra harray <<< "$(head -n1 ${geno_loc})" #Split title columns in the header

IFS=$'\n'; shorten=($(printf '%s\n' "${harray[@]}" |sed '/.Top Alleles/!d')) #Keep only titles containing ".Top Alleles"
ids=(${shorten[@]%.Top Alleles}) # Remove the ".Top Alleles" tag leaving only the sample IDs

#xargs echo < ${ids[@]} > ${geno_id_loc}
for i in "${ids[@]}"; do echo $i; done  > ${geno_id_loc}

# Transfer, check and zip the Raw Genome Studio Export
move_and_check ${geno_loc} ${genotype_folder}

#cp ${geno_loc} ${genotype_folder}
#geno_md5=$(cat ${geno_loc} | md5sum)
#geno_check=$(cat ${genotype_folder}/*.txt | md5sum)
#[ "${geno_md5}" = "${geno_check}" ] && echo  "Raw Genome Studio Export Successfully Copied" || echo "Raw Genome Studio Export's md5sum after copy has changed! Please Check"
gzip ${genotype_folder}/$(basename ${geno_loc}) 
fi
#FOR TESTING PURPOSES RESTRICT THE SIZE OF THE TRACKING SHEET

#head -n39 ${tracking_sheet} > ${test_track}
#tracking_sheet=${test_track}
#path_list_loc=${test_list}
#echo $tracking_sheet
#echo $path_list_loc

#Transfer and check the md5sums of the plink files
#snp_count=$(cat ${plink_bim} | wc -l | cat)
if [ "${plink_base}" == "nfjksl" ]; then
for ext in "${plink_exts[@]}"; do
	plink_file=${plink_base}"${ext}"
	move_and_check ${plink_file} ${plink_folder}
	#cp "${plink_file}" ${plink_folder}
	#plink_check=$(cat "${plink_file}" | md5sum)
	#plink_md5=$(cat ${plink_folder}/*"${ext}" | md5sum)
	#[ "${plink_md5}" = "${plink_check}" ] && echo  "Plink ${ext} File Successfully Copied" || echo "Plink ${ext} File's md5sum after copy has changed! Please Check"
done

fi






## COPY/CHECK/TAR IDATS
cd ${idat_loc}
python /home/jkb4y/h4t1/programs/packaging_data.py ${tracking_sheet} ${path_list_loc} ${geno_id_loc} ${plink_fam} # 1. make list of idat files using tracking sheet info, compare ids with the genotype and plink files, and check for duplicates
#exit
tar --create -Wf ${tar_loc} --files-from ${path_list_loc} # 2a. tar up iDats from list
md5sum ${tar_loc} > ${tar_md5_loc} # 2b. create md5sum of tar for end user

python /home/jkb4y/h4t1/programs/tarsum.py ${tar_loc} -o ${md5_loc} # 3. create tarsum of transferred idats

md5sum --quiet --check ${md5_loc} # 4. check the tarsum against md5sums of original files

gzip ${tar_loc} -q # 5. Zip the tar
cd ${basefolder}
##OLD STUFF
# make a checksum file of those files: xargs md5sum < ${list_loc} > ${md5_loc}

#ALTERNATIVELY, we could both tar and md5sum at the same time:
#tar -cWvpf ${tar_loc} --files-from ${list_loc} | xargs  -I '{}' sh -c "test -f '{}' && md5sum '{}'" | tee ${md5_loc}


##### WRITE THE README FILE  #####
echo ${readme}
cat > ${readme} << EOF
                                  
README  

Center for Public Health Genomics, UVA

${chiplong} Bead Array ${study} Release

Release Date: ${release_date:0:4}/${release_date:4:2}/${release_date:6:2}


CONTENTS

This drive contains raw data from the ${study} ${chiplong} Project.
Number of Samples: $(echo " ${track_count}" | sed -r ':L;s=\b([0-9]+)([0-9]{3})\b=\1,\2=g;t L')
Number of SNPs: ${snp_count}


DIRECTORIES
The drive holds the following directories and files

1. $(basename ${illumina_folder})/
	A. Illumina ${chiplong} cluster file (.egt)
	B. Illumina ${chiplong} manifest files (.bpm and .csv)

2. $(basename ${plink_folder})/
	A. Genotype files in PLINK binary format

3. $(basename ${idat_folder})/
	A. Zipped Tar of intensity data files (*.idat)
	B. File containing the MD5 checksum of the zipped tar
	C. File containing the list of idat files
	D. File containing the MD5 checksums of the idat files

4. $(basename ${genotype_folder})/
	A. Raw genotype data exported from GenomeStudio using GeneTrain2

5. $(basename ${tracking_sheet})



NOTES



EOF






