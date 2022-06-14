#!/bin/bash

CWD=$(pwd)
echo $CWD

endings=(*.g *.g.2.* *.g.3.* *.g.4.* *.post *.blot *.e *.e.4.*)

for ending in "${endings[@]}"
do 
    # echo $ending
    matches=$(grep -r -I -H --include="${ending}" "../${ending}" .)

    for lines in ${matches}
    do
        # echo ${lines}
        IFS=':'
        read -a strarr <<< "$lines"
        file=${strarr[0]}
        link=${strarr[1]}

        dir="$(dirname "${file}")"
        file="$(basename "${file}")"

        echo Relink ${link} ${dir} ${file}
        
        cd $dir
        rm ${file}
        ln -s ${link} ${file}
        cd ${CWD}

    done
done

# for d in test/*/*/
# do      
#     (cd "$d" && dos2unix *.comp)
# done
