#!/bin/bash -i

# Instructions for running this Bash script:
# 1. Make sure you have the necessary permissions to execute this script.
# 2. Open a terminal and navigate to the directory where this script is saved.
# 3. Run the script in the backgorund of server by typing the following command:
#    nohup bash Step01_run_first.sh &
# 4. Follow the on-screen instructions and prompts, if any, during the script execution.
# 5. Ensure you have the required dependencies and data files in the same directory as this script.
# 6. Contact the script author for any issues or questions.

# Enable the expansion of aliases in the script to allow for alias usage.
shopt -s expand_aliases


# Define variables for different identifiers and values used in the script:
# - 'train' for the training set identifier 'gen9';
# - 'valid' for the validation set identifier 'gen10';
# - 'win50', 'win100', 'win200', 'win400' for window sizes '50win', '100win', '200win', '400win' respectively;
# - 'topfst' for the top Fst threshold identifier 'fst99'.
train="gen9"; valid="gen10"
win50="50win"; win100="100win"; win200="200win"; win400="400win";topfst="fst99"


## Calculating allele freq for each sub-population
gfortran R_Fortran_Bash/MAF_S1pop.f90;./a.out
gfortran R_Fortran_Bash/MAF_S2pop.f90;./a.out
gfortran R_Fortran_Bash/MAF_S1S2pop.f90;./a.out

## Calculate Fst scroe based on allele freq from low and high 5% phenotype in generation 9
R4.1.3 CMD BATCH R_Fortran_Bash/calc_fst.R; rm .RData 

## Randomly select 5K animals from generation 10 for testing set
R4.1.3 CMD BATCH R_Fortran_Bash/split_train_test.R ; rm .RData

##calculate genetic percentage for each QTL
bash R_Fortran_Bash/QTLeffect.sh
R4.1.3 CMD BATCH R_Fortran_Bash/additive_QTL.R; rm .RData

## Algorithm for identifying index of SNPs randomly selected within SNP-windows associated with QTLs using FST scores 
R4.1.3 CMD BATCH R_Fortran_Bash/NumQTLsPassed_50SNPwin_12Rnd.R; rm .RData
R4.1.3 CMD BATCH R_Fortran_Bash/NumQTLsPassed_100SNPwin_12Rnd.R; rm .RData
R4.1.3 CMD BATCH R_Fortran_Bash/NumQTLsPassed_200SNPwin_12Rnd.R; rm .RData
R4.1.3 CMD BATCH R_Fortran_Bash/NumQTLsPassed_400SNPwin_12Rnd.R; rm .RData


## list of animals in training set (generation 9)
cut -d " " -f 1 geno.train.$train   > train.${train}.ID


## Create genotype file for randomly selected 5K animals in generation 10 as testing set
sort  +0 -1 geno.gen910 > temp1
sort  +0 -1 test.${valid}.ID > temp2
join -1 1 -2 1 temp1 temp2 > geno.test.$valid
rm temp* 


## Extract selected SNPs based on Fst threshold from the original 600K in generation 9 (training set)
awk '{$1=""}1' geno.train.$train > mrk.temp
cut -d " " -f 1 geno.train.$train > ID.$train

## Prioritize SNPs based on Fst threshold 99% in training set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' ${topfst}.ID  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$train aa > bb

cat bb | sed "s/\t\t*/ /g" > marker${topfst}_$train
rm aa bb


## Prioritize random 12 SNPs within 50 SNP windows based on Fst scroe around each QTL for training set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' Rnd12SNPsFst_$win50  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$train aa > bb

cat bb | sed "s/\t\t*/ /g" > Rnd12SNPsFst_${win50}_$train
rm aa bb


## Prioritize random 12 SNPs within 100 SNP windows based on Fst scroe around each QTL for training set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' Rnd12SNPsFst_$win100  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$train aa > bb

cat bb | sed "s/\t\t*/ /g" > Rnd12SNPsFst_${win100}_$train
rm aa bb


## Prioritize random 12 SNPs within 200 SNP windows based on Fst scroe around each QTL for training set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' Rnd12SNPsFst_$win200  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$train aa > bb

cat bb | sed "s/\t\t*/ /g" > Rnd12SNPsFst_${win200}_$train
rm aa bb


## Prioritize random 12 SNPs within 400 SNP windows based on Fst scroe around each QTL for training set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' Rnd12SNPsFst_$win400  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$train aa > bb

cat bb | sed "s/\t\t*/ /g" > Rnd12SNPsFst_${win400}_$train

rm aa bb mrk.temp ID.$train


##------------------------------------------------------------------------------------------------------------##


## Extract selected SNPs based on Fst threshold from the original 600K in generration 10 (testing set)
awk '{$1=""}1' geno.test.$valid > mrk.temp
cut -d " " -f 1 geno.test.$valid > ID.$valid

## SNPs based on Fst threshold 99% in testing set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' ${topfst}.ID  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$valid aa > bb

cat bb | sed "s/\t\t*/ /g" > marker${topfst}_$valid
rm aa bb

## Prioritize random 12 SNPs within 50 SNP windows based on Fst scroe around each QTL for validation set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' Rnd12SNPsFst_$win50  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$valid aa > bb

cat bb | sed "s/\t\t*/ /g" > Rnd12SNPsFst_${win50}_$valid
rm aa bb


## Prioritize random 12 SNPs within 100 SNP windows based on Fst scroe around each QTL for validation set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' Rnd12SNPsFst_$win100  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$valid aa > bb

cat bb | sed "s/\t\t*/ /g" > Rnd12SNPsFst_${win100}_$valid
rm aa bb


## Prioritize random 12 SNPs within 200 SNP windows based on Fst scroe around each QTL for validation set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' Rnd12SNPsFst_$win200  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$valid aa > bb

cat bb | sed "s/\t\t*/ /g" > Rnd12SNPsFst_${win200}_$valid
rm aa bb


## Prioritize random 12 SNPs within 400 SNP windows based on Fst scroe around each QTL for validation set
awk 'FNR==NR {C[++j]=$1;next} {for (i=1;i<=j;i++) printf "%s ", $C[i]; printf "\n"}' Rnd12SNPsFst_$win400  mrk.temp > aa

awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }'  ID.$valid aa > bb

cat bb | sed "s/\t\t*/ /g" > Rnd12SNPsFst_${win400}_$valid


rm aa bb mrk.temp ID.$valid


