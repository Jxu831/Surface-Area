#!/bin/bash
# Hailey Wallace

rm tmp.csv

# Command line options
while getopts ":hi:l:o:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) # Enter the coords file
         i=$OPTARG;;
      l) # Enter the ligand mol2 file
         l=$OPTARG;;
      o) # Enter the ouput file name
         o=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

echo "Making ligand CG file..."
# Change the ligand coordinate file into a tab delimited tmp.csv file
while read -r p; do
  group=$(echo $p | awk -F'(' '{print $1}')
  echo $p | awk -F'|' '{print $2}' | sed 's/\;/\n/g' | \
  awk -F '\t' -v group2=$group '{print $0","group2}'| \
  sed -e 's/(/\n/g' | sed -e 's/)//' >> tmp.csv #sed 's/,/\t/g'
done < $i

# Read the coordinate file by chunks separated by newline; feed into array
i=1
s=1
declare -a arr
while read -r line
do

    # Increase counter (i) if there is an empty line. Set s to 1.
    [[ $line == "" ]] && ((i++)) && s=1 && continue
    # If s=0, then it is the next block.
    # Value of array os the previous value + current line
    [[ $s == 0 ]] && arr[$i]="${arr[$i]}
$line" || {
            # Otherwise the value of array is current line & s=0
            arr[$i]="$line"
            s=0;
    }
done < tmp.csv

# Reading in corresponding CG atoms
conh2=$(echo "combs_groups/conh2_ASN.txt" "combs_groups/conh2_GLN.txt")
bbcco=$(echo "combs_groups/bbcco_GLY.txt" "combs_groups/bbcco_ALA.txt" "combs_groups/bbcco_PRO.txt")
ph1=$(echo "combs_groups/ph_PHE_1.txt")
ph2=$(echo "combs_groups/ph_PHE_2.txt") # ph numbering in opposite direction
bb_cnh=$(echo "combs_groups/bb_cnh_GLY.txt" "combs_groups/bb_cnh_LYS.txt" "combs_groups/bb_cnh_LYS.txt")
ccn=$(echo "combs_groups/ccn_LYS.txt")
ccoh=$(echo "combs_groups/ccoh_SER.txt" "combs_groups/ccoh_THR.txt")
coh=$(echo "combs_groups/coh_SER.txt" "combs_groups/coh_THR.txt")
coo=$(echo "combs_groups/coo_ASP.txt" "combs_groups/coo_GLU.txt")
csc=$(echo "combs_groups/csc_MET.txt")
csh=$(echo "combs_groups/csh_CYS.txt")
gn=$(echo "combs_groups/gn_ARG.txt")
hid=$(echo "combs_groups/hid_HIS.txt")
hie=$(echo "combs_groups/hie_HIS.txt")
hip=$(echo "combs_groups/hip_HIS.txt")
indole=$(echo "combs_groups/indole_TRP.txt")
phenol1=$(echo "combs_groups/phenol_TYR_1.txt")
phenol2=$(echo "combs_groups/phenol_TYR_2.txt") # phenol numbering in opposite direction
isopropyl1=$(echo "combs_groups/isopropyl_LEU_1.txt" "combs_groups/isopropyl_VAL_1.txt")
isopropyl2=$(echo "combs_groups/isopropyl_LEU_2.txt" "combs_groups/isopropyl_VAL_2.txt")  # isopropyl numbering in opposite direction
pro1=$(echo "combs_groups/pro_PRO_1.txt")
pro2=$(echo "combs_groups/pro_PRO_2.txt") # proline numbering in opposite direction
ch3=$(echo "combs_groups/ch3_ALA.txt" "combs_groups/ch3_ILE.txt")


# Loop through the array and run the ligand support script for corresponding CG atoms
for i in "${arr[@]}"
do
# bb_cco groups
   if [[ "$i" == *"bb_cco"* ]]; then
      echo "################ C -> CA "
      for a in $bbcco; do
         printf '%s\n' $i | paste -d',' - $a | grep -v 'blank'
      done
 # alternate bb_cco groups
      echo "################ alternate N -> CA "
      for a in $bbcco; do
         printf '%s\n' $i | paste -d',' - $a | grep -v 'CA' | sed 's/blank/CA/g'
      done

# conh2 groups
   elif [[ "$i" == *"conh2"* ]]; then
      echo "################"
      for a in $conh2; do
         printf '%s\n' $i | paste -d',' - $a
      done

# ph groups
   elif [[ "$i" == *"ph"* ]]; then
      echo "################"
      for a in $ph1; do
         printf '%s\n' $i | paste -d',' - $a
      done

# ph groups - opposite direction
      echo "################"
      for a in $ph2; do
         printf '%s\n' $i | paste -d',' - $a
      done

# bb_cnh groups
   elif [[ "$i" == *"bb_cnh"* ]]; then
      echo "################"
      for a in $bb_cnh; do
         printf '%s\n' $i | paste -d',' - $a
      done

# ccn groups
   elif [[ "$i" == *"ccn"* ]]; then
      echo "################"
      for a in $ccn; do
         printf '%s\n' $i | paste -d',' - $a
      done

# ccoh groups
   elif [[ "$i" == *"ccoh"* ]]; then
      echo "################"
      for a in $ccoh; do
         printf '%s\n' $i | paste -d',' - $a
      done

# coh groups
   elif [[ "$i" == *"coh"* ]]; then
      echo "################"
      for a in $coh; do
         printf '%s\n' $i | paste -d',' - $a
      done

# coo groups
   elif [[ "$i" == *"coo"* ]]; then
      echo "################"
      for a in $coo; do
         printf '%s\n' $i | paste -d',' - $a
      done

# csc groups
   elif [[ "$i" == *"csc"* ]]; then
      echo "################"
      for a in $csc; do
         printf '%s\n' $i | paste -d',' - $a
      done

# csh groups
   elif [[ "$i" == *"csh"* ]]; then
      echo "################"
      for a in $csh; do
         printf '%s\n' $i | paste -d',' - $a
      done

# gn groups
   elif [[ "$i" == *"gn"* ]]; then
      echo "################"
      for a in $gn; do
         printf '%s\n' $i | paste -d',' - $a
      done

# hid groups
   elif [[ "$i" == *"hid"* ]]; then
      echo "################"
      for a in $hid; do
         printf '%s\n' $i | paste -d',' - $a
      done

# hie groups
   elif [[ "$i" == *"hie"* ]]; then
      echo "################"
      for a in $hie; do
         printf '%s\n' $i | paste -d',' - $a
      done

# hip groups
   elif [[ "$i" == *"hip"* ]]; then
      echo "################"
      for a in $hip; do
         printf '%s\n' $i | paste -d',' - $a
      done

# indole groups
   elif [[ "$i" == *"indole"* ]]; then
      echo "################"
      for a in $indole; do
         printf '%s\n' $i | paste -d',' - $a
      done

# phenol groups
   elif [[ "$i" == *"phenol"* ]]; then
      echo "################"
      for a in $phenol1; do
         printf '%s\n' $i | paste -d',' - $a
      done

# phenol groups - opposite direction
   #elif [[ "$i" == *"phenol"* ]]; then
      echo "################ - opposite direction"
      for a in $phenol2; do
         printf '%s\n' $i | paste -d',' - $a
      done

# isopropyl groups
   elif [[ "$i" == *"isopropyl"* ]]; then
      echo "################"
      for a in $isopropyl1; do
         printf '%s\n' $i | paste -d',' - $a
      done

# isopropyl groups - opposite direction
  # elif [[ "$i" == *"isopropyl"* ]]; then
      echo "################ - opposite direction"
      for a in $isopropyl2; do
         printf '%s\n' $i | paste -d',' - $a
      done

# pro groups
   elif [[ "$i" == *"pro"* ]]; then
      echo "################"
      for a in $pro1; do
         printf '%s\n' $i | paste -d',' - $a | grep -v 'blank'
      done

# pro groups - opposite direction
   #elif [[ "$i" == *"pro"* ]]; then
      echo "################ - opposite direction"
      for a in $pro2; do
         printf '%s\n' $i | paste -d',' - $a | grep -v 'blank'
      done

# ch3 groups
   elif [[ "$i" == *"ch3"* ]]; then
      echo "################"
      for a in $ch3; do
         printf '%s\n' $i | paste -d',' - $a
      done

# default
   else
      echo "################"
      echo "PARSING ERROR: no combs groups in CG coordinate file. Re-run ligand_matcher.py"
   fi

done > tmp_cg_coordinates.csv

# turn this output csv into a second array to work with CG in 'chunks' aka multiple arrays
p=1
k=1
declare -a arr2
while read -r line
do
    [[ $line == *"################"* ]] && ((p++)) && k=1 && continue
    [[ $k == 0 ]] && arr2[$p]="${arr2[$p]}
$line" || {
            arr2[$p]="$line"
            k=0;
    }
done < tmp_cg_coordinates.csv

#rm tmp8.csv
echo " \n" > tmp8.csv
ligand_code=$(grep -A1 'ATOM' $l | grep -v 'ATOM' | awk -F' ' '{print $8}' | uniq)

COUNTER=0
for v in "${arr2[@]}"; do
   vs=$(printf '%s\n' $v)
   let COUNTER=COUNTER+1
   for i in $vs; do
      x_coord=$(echo $i | awk -F',' '{print $1}')
      test=$(grep -e $x_coord $l | awk -F' ' '{print $2}')
      other_info=$(echo $i | awk -F',' '{print $6,$5,$4}')

      # if the coordinate is in the ligand file, then add it to the output file
      # this removes the excess hydrogens atoms from rdkit

      if grep -Fqe $x_coord $l
      then
         echo $ligand_code $test $other_info $COUNTER
      fi
   done
   echo ""

done > tmp8.csv

#rm tmp9.csv
awk '{print $6,$2}' tmp8.csv | sort | uniq | awk '{if(a!=$1) {a=$1; printf "\n%s%s",$0,FS} else {a=$1;$1="";printf $0 }} END {printf "\n"}' | sed '/^$/d' | sed 's/ /,/' > tmp9.csv

#rm tmp10.csv

awk -F, '
{ # first pass
  groups[ $1 ] =  $2
  split( $2, items, " " )
  for (item in items) {
    items_in_groups[ items[item] ] = items_in_groups[ items[item] ] " " $1
  }
}
END {  # second pass
  for ( group in groups ) {
    split( groups[group], items, " " )
    for ( item in items ) {
      split( items_in_groups[ items[item] ] , in_groups, " " )
      if (length( in_groups ) >= 2)  {
        print group,groups[group],in_groups[1]
        break
      }
    }
  }
}
' tmp9.csv > tmp10.csv

if [ -z "$o" ]; then
   if [ -f "ligand.txt" ]; then
      if [ -d "./old_files" ]; then
         echo "Moving previous ligand.txt file to ./old_files/"
         mv ligand.txt ./old_files/
      else
         echo "Creating and moving previous ligand.txt file to ./old_files/"
         mkdir ./old_files
         mv ligand.txt ./old_files/
      fi
   fi
   while read -r line; do
      group_num=$(echo $line | awk -F' ' '{print $1}')
      cluster_num=$(echo $line | awk -F' ' '{print $NF}')
      awk -v a="$group_num" -v b="$cluster_num" 'a == $6 {print $0,b}' tmp8.csv
      echo ""
   done < tmp10.csv > ligand.txt
else
   while read -r line; do
      group_num=$(echo $line | awk -F' ' '{print $1}')
      cluster_num=$(echo $line | awk -F' ' '{print $NF}')
      awk -v a="$group_num" -v b="$cluster_num" 'a == $6 {print $0,b}' tmp8.csv
      echo ""
   done < tmp10.csv > $o
fi


#rm tmp.csv ; rm tmp10.csv ; rm tmp8.csv ; rm tmp9.csv ; rm tmp_cg_coordinates.csv
