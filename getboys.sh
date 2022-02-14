#!/bin/bash
# author m.fuse
# exe
cubg=${GAUSS_EXEDIR}"/cubegen"
cubm=${GAUSS_EXEDIR}"/cubman"

### Functions ###
Help()
{
   # Display Help
   echo "Generate centroid from a Boys localized fchk"
   echo "and save them in a xyz file."
   echo "See https://doi.org/10.1039/C4DT00251B "
   echo " or "
   echo "https://www.youtube.com/watch?v=ed3PFYVtha8"
   echo
   echo "Syntax: getboys.sh filename.fchk MO [...]"
   echo 
   echo "filename.fchk   The Gaussian fchk with the orbitals"
   echo "MO              The integer referring to the MO to be localized"
   echo "options:"
   echo "h     Print this Help."
   echo "t     Print a gjf template to get the data required"
   echo
}
check_integer (){
        if [ "$1" -eq "$1" ] 2>/dev/null
        then
                echo 0
        else
                echo 1
        fi
}

print_tmpl () {
        printf "%s\n\n" "$template1" > template_boys.gjf
}

### TEMPLATE ###
template1="$(cat <<EOF
%mem=60GB
%nprocshared=16
%chk=CHK
#P b3lyp/cc-pVDZ int=ultrafine NoSymm

first step

CHR SPIN
COORD

--Link1--
%mem=60GB
%nprocshared=16
%chk=CHK
#P b3lyp/cc-pVDZ int=ultrafine NoSymm
 Geom=Check Guess(NoSymm,Read,Local,Save,Only)
 Pop=Full

Boys localization

CHR SPIN
EOF
)"

### MAIN ###
while getopts ":h:t" option; do
    case $option in
        h) Help
           exit;;
        t) echo "Printing a template in template_boys.gjf"
           print_tmpl
           exit;;
        \?) echo "Error: Invalid option"
           exit;;
    esac
done

argv=("$@")
nargv=${#argv[@]}
if [ $nargv -lt 2 ]
then
        echo "not enough argument"
        exit 1
fi
#Checks if the file exist
FILE=${argv[0]}
if ! [ -f "$FILE" ]; then
    echo "$FILE does not exist."
    exit 1
fi
mol=${FILE//\.fchk/}
#Checks the integers
for i in ${argv[@]:1:nargv}
do
        isint=$( check_integer $i )
        if [[ $isint -eq 0 ]]
        then
                mos+=($i)
        else
                echo $i "Not a number."
                Help
        fi
done

for i in ${mos[@]}
do
        num=$( printf "%4d" $i )
        echo " #### Orbital:${num} ####"
        echo " Making the cube"
        $cubg 1 MO=$i ${mol}.fchk ${mol}_MO${i}.cube -2 h 1>/dev/null 
        echo " Computing the square"
        printf '%s\n' SQ ${mol}_MO${i}.cube y ${mol}_MO${i}_SQ.cube y | $cubm 1>/dev/null
        echo " Getting the centroid"
        printf '%s\n' P ${mol}_MO${i}_SQ.cube y | $cubm |grep "DipAE=" | awk '{ printf " X %12.6f %12.6f %12.6f\n", -$2*0.529177, -$3*0.529177, -$4*0.529177}' >> ${mol}_VS.xyz
        rm ${mol}_MO${i}.cube ${mol}_MO${i}_SQ.cube

done
