#!/bin/bash

CWD=$(pwd)
echo $CWD

for d in test/verification/*/*/ ; do

    echo "$d"
    cd $d

    endings=("*.g" "*.g.2.*" "*.g.3.*" "*.g.4.*" "*.post" "*.blot" "*.e" "*.e.4.*" "*.txt")

    for ending in "${endings[@]}"
    do 
        # echo $ending
        matches=$(grep -r -I -H --include="${ending}" "../${ending}" .)

        # echo  matches: $matches end
        IFS='
'
        for lines in ${matches}
        do
            # echo lines: ${lines} end
            IFS=':'
            read -a strarr <<< "$lines"
            file=${strarr[0]}
            link=${strarr[1]}

            dir="$(dirname "${file}")"
            file="$(basename "${file}")"

            echo Relink ${link} ${dir} ${file}
            
            # cd $dir
            rm ${file}
            ln -s ${link} ${file}

        done
    done

    cd $CWD
done

for d in test/regression/*/*/ ; do

    echo "$d"
    cd $d

    endings=("*.g" "*.g.2.*" "*.g.3.*" "*.g.4.*" "*.post" "*.blot" "*.e" "*.e.4.*" "*.txt")

    for ending in "${endings[@]}"
    do 
        # echo $ending
        matches=$(grep -r -I -H --include="${ending}" "../${ending}" .)

        # echo  matches: $matches end
        IFS='
'
        for lines in ${matches}
        do
            # echo lines: ${lines} end
            IFS=':'
            read -a strarr <<< "$lines"
            file=${strarr[0]}
            link=${strarr[1]}

            dir="$(dirname "${file}")"
            file="$(basename "${file}")"

            echo Relink ${link} ${dir} ${file}
            
            # cd $dir
            rm ${file}
            ln -s ${link} ${file}

        done
    done

    cd $CWD
done

for d in test/*/*/
do      
    (cd "$d" && dos2unix *.comp)
done

#ERROR: parsing input file, currently at --> comp File CLRF -> LF