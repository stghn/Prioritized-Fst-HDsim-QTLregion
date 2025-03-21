#!/bin/bash

# Instructions for running this Bash script:
# 1. Make sure you have the necessary permissions to execute this script.
# 2. Open a terminal and navigate to the directory where this script is saved.
# 3. Run the script in the backgorund of server by typing the following command:
#    nohup bash Step01_run_first.sh &
# 4. Follow the on-screen instructions and prompts, if any, during the script execution.
# 5. Ensure you have the required dependencies and data files in the same directory as this script.
# 6. Contact the script author for any issues or questions.


# Define variables for different identifiers and values used in the script:
# - 'train' for the training set identifier 'gen9';
# - 'valid' for the validation set identifier 'gen10';
train="gen9";valid="gen10"

## Skip headers of p1_data_001.txt and extract ID animal and phenotype for all animals
tail -n +2 output_rep1/p1_data_001.txt | awk '{print $1,$10}' > pheno

## Add two fixed effects (fix1=100 levels, fix2=4 levels)
gfortran -c R_Fortran_Bash/kind.f90
gfortran R_Fortran_Bash/fixed_effect.f90;./a.out
mv fort.11 data.pheno

## Extract phenotypes for generation 9 ##
tail -n +2 output_rep1/p1_data_001.txt | awk '{print $1,$5}' | awk '$2>8 && $2<10' | awk '{print $1}' > ID.$train

sort  +2 -3 data.pheno > temp1
sort  +0 -1 ID.$train > temp2
join -1 3 -2 1 temp1 temp2 | awk '{print $2,$3,$1,$4}' | column -t -s " " > data.$train

rm pheno temp* ID.$train

## Sort by phenotype column(top file lowest and bottom file highest) 
sort -g -k4 data.$train > sort.data.$train

## Create two subpopulation from generation 9 (5% of low and high phenotypes in generation 5)
awk 'NR==1,NR==750' sort.data.$train > low.pheno.$train 
awk 'NR==14251,NR==15000' sort.data.$train  | sort -g -k4r > high.pheno.$train 

## Replace SNP codes from QMsim simulation to SNP code 0, 1, and 2
tail -n +2 output_rep1/p1_mrk_001.txt | awk '{$1=""}1' > geno.temp
tail -n +2 output_rep1/p1_mrk_001.txt | cut -d " " -f 1 > ID
cat geno.temp | sed "s/ */ /g" > geno

## Compling the replacement SNP code 
gfortran R_Fortran_Bash/ReplaceSNP.f90;./a.out

## Remove extra spaces between genetic code (0,1,2)
cat new.geno | sed "s/ */ /g" > new1.geno

## Join ID animals with genotype file ##
awk 'NR==FNR { a[c=FNR]=$0; next }
   { printf "%-8s\t%s\n", a[FNR], $0 }
 END { for(i=FNR+1;i<=c;i++) print a[i] }' ID new1.geno > new2.geno

rm ID geno geno.temp new1.geno new.geno sort.data.$train 

## Remove tap space between id animal and genetic code(0,1,2)
cat new2.geno | sed "s/\t\t*/ /g" > geno.gen910
rm new2.geno

## Create genotype file for generation 9
awk '{print $3}' data.$train  > ID.$train 

sort  +0 -1 geno.gen910 > temp1
sort  +0 -1 ID.$train  > temp2
join -1 1 -2 1 temp1 temp2 > geno.train.$train 

rm temp* ID.$train 


## Create genotype file for two sup-population from generation 9 (low and high phenotype)
sort  +2 -3 low.pheno.$train  > temp1
sort  +0 -1  geno.train.$train  > temp2
join -1 3 -2 1 temp1 temp2 > sub_low_$train 
rm temp*

sort  +2 -3 high.pheno.$train  > temp1
sort  +0 -1  geno.train.$train  > temp2
join -1 3 -2 1 temp1 temp2 > sub_high_$train 
rm temp*

## Remove columns related to fixed effect levels and phenotype record 
awk 'BEGIN{FS=OFS=" "}{$2=$3=$4=" "}{print}' sub_low_$train  >  S1.low.$train 
awk 'BEGIN{FS=OFS=" "}{$2=$3=$4=" "}{print}' sub_high_$train > S2.high.$train 

cat S1.low.$train  S2.high.$train  > S1S2.$train 

rm sub_*
