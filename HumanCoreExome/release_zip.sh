#!/bin/bash
freezenum=$1


writefolder=/t121/hellacious/projects/HELLACIOUS/workfiles/JessicaQC_Final_Zips
qc_final_folder=/t121/hellacious/projects/HELLACIOUS/workfiles/Jessica/QC_final
reference=/t121/hellacious/projects/HELLACIOUS/workfiles/Jessica/Reference

if [[ ${freezenum} -eq 1 ]]; then
    password_key=${writefolder}/folders_passwords_df1.txt
    passfile=${reference}/passwords_df1.txt
    folderfile=${reference}/studyfolders_df1.txt
    outfolder=${writefolder}/dataFreeze1
elif [[ ${freezenum} -eq 2 ]]; then
    password_key=${writefolder}/folders_passwords_df2.txt
    passfile=${reference}/passwords_df2.txt
    folderfile=${reference}/studyfolders_df2.txt
    outfolder=${writefolder}/dataFreeze2
else
echo -e "${freezenum} is not a Date Freeze Number I recognize."
fi

mkdir ${outfolder}

#ls ${qc_final_folder} | sed 's/\t/\n/g' | grep -v .txt | grep -v .sh > studyfolders.txt


IFS=$'\n' read -d '' -r -a passwords < ${passfile}
IFS=$'\n' read -d '' -r -a studies < ${folderfile}
count=$(cat ${folderfile} | wc -l)

echo -e Folder'\t'Password > ${password_key}
for i in $(seq ${count}); do
    index=$(($i-1))
    echo $index
    study=${studies[${index}]}
    echo $study
    password=${passwords[${index}]}
    echo $password
    tar  --directory ${qc_final_folder}  -cvf ${outfolder}/${study}.tar ${study}
    echo -e "${study}\t${password}" >> ${password_key}
    md5sum ${qc_final_folder}/${study}/* > ${outfolder}/${study}.md5
    tar  --directory ${outfolder} -vuf ${outfolder}/${study}.tar ${study}.md5 
    #gzip gzip -cv $writefolder/$study.tar | openssl enc -aes-256-cbc -e -k $password > $writefolder/$study.tar.gz.enc

    curdir=$(pwd)
    cd ${outfolder}
    zip -P "${password}" ${study}.zip ${study}.tar
    cd ${curdir}
    rm -f ${outfolder}/${study}.tar
  done



