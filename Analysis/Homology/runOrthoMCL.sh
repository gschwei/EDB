# created by Gabriel on 23/02/17
# last modified by Gabriel on 23/02/17

# run OrthoMCL with varying inflation values.

step="0.1"

for (( g=11; $g<$(bc<<<"10.1/$step"); g++ )); do
    echo "run OrthoMCL with I =" $g
    mcl ../../../Blast/Blastp/DbMergedSecretomeQuMergedSecretomeShort.csv --abc -I $g -o "Mcl_I"$g"_out.txt"
    echo $g 'is done. Go to next value. Bye.' 
done
