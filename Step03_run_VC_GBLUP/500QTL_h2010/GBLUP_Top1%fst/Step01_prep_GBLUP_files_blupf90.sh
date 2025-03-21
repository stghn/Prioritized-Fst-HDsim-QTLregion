#!/bin/bash -i  

# Instructions for running this Bash script:
# 1. Make sure you have the necessary permissions to execute this script.
# 2. Open a terminal and navigate to the directory where this script is saved.
# 3. Run the script in the backgorund of server by typing the following command:
#    nohup bash Step01_prep_GBLUP_files_blupf90.sh &
# 4. Follow the on-screen instructions and prompts, if any, during the script execution.
# 5. Ensure you have the required dependencies and data files in the same directory as this script.
# 6. Contact the script author for any issues or questions.

## create fake pedigree from training and testing population for renum parameter file
awk '{print $1}' ../data.gen9 > train.gen9.ID
cat train.gen9.ID ../test.gen10.ID > ID.animal
cat ID.animal | sed 's/^/0 /g' > temp1
cat temp1 | sed 's/^/0 /g' > temp2
awk 'BEGIN {FS=OFS=" "} {temp=$3; $3=$1; $1=temp} {print}' temp2 > fake.ped.GBLUP
rm temp* ID*

## create SNP file from training and testing population for renum parameter file 

# SNP data for training population
awk '{print $1}' ../markerfst99_gen9 > ID.gen9
awk '{$1=""}1' ../markerfst99_gen9 > g1
sed -e "s/[ <tab>]*//g" g1 > g2

awk 'NR==FNR { a[c=FNR]=$0; next }
    { printf "%-8s\t%s\n", a[FNR], $0 }
   END { for(i=FNR+1;i<=c;i++) print a[i] }' ID.gen9 g2 > g3

cat g3 | sed "s/\t\t*/ /g" > SNP.train.gen9
rm g* 

# SNP data for testing population
awk '{print $1}' ../markerfst99_gen10 > ID.gen10
awk '{$1=""}1' ../markerfst99_gen10 > g1
sed -e "s/[ <tab>]*//g" g1 > g2

awk 'NR==FNR { a[c=FNR]=$0; next }
    { printf "%-8s\t%s\n", a[FNR], $0 }
   END { for(i=FNR+1;i<=c;i++) print a[i] }' ID.gen10 g2 > g3

cat g3 | sed "s/\t\t*/ /g" > SNP.test.gen10

rm g* 

## combine SNP data of training and testing set
cat SNP.train.gen9 SNP.test.gen10 > SNP.GBLUP
rm SNP.t*

## Extract data related to generation 9 and 10 for the use of TBV
tail -n +2 ../output_rep1/p1_data_001.txt | awk '$5>8' > full.data.geno910

## create datafile for renumf90 parameter
cat ID.gen9 ID.gen10 > ID.gen910

sort  +2 -3 ../data.pheno > temp1
sort  +0 -1  ID.gen910 > temp2
join -1 3 -2 1 temp1 temp2 > data.gen910
rm temp* ID.*

## create map file for using in blupf90
awk 'NR>1' ../output_rep1/lm_mrk_001.txt |  awk '{print $1,$2,$3=sprintf("%10.0f",$3*1000000)}' > mm1 
awk '{print "M"$1}' ../fst99.ID > mm2

sort  +0 -1 mm1 > temp1
sort  +0 -1  mm2 > temp2
join -1 1 -2 1 temp1 temp2 | sort -V > mapfile
rm temp* mm*
sed  -i '1i  SNP_ID  CHR  POS' mapfile

