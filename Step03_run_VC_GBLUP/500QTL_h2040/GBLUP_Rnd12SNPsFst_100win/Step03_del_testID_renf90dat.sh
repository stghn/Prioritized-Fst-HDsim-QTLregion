#!/bin/bash -i

# Instructions for running this Bash script:
# 1. Make sure you have the necessary permissions to execute this script.
# 2. Open a terminal and navigate to the directory where this script is saved.
# 3. Run the script in the backgorund of server by typing the following command:
#    nohup bash Step03_del_testID_renf90dat.sh &
# 4. Follow the on-screen instructions and prompts, if any, during the script execution.
# 5. Ensure you have the required dependencies and data files in the same directory as this script.
# 6. Contact the script author for any issues or questions.

## sort renumber animal IDs
sort -n +3 -4 renf90.dat > sort_data
awk '{print $1 ,$10}' renadd03.ped > renumid

## merge renumber animals with original IDs in traning set
sort +0 -1 ../train.gen9.ID > temp1
sort +1 -2 renumid > temp2
join -1 1 -2 2 temp1 temp2 > renum_train
rm temp*

## subset training data from the renf90.dat generated from renumf90 program
sort +3 -4 sort_data > temp1
sort +1 -2 renum_train > temp2
join -1 4 -2 2 temp1 temp2 > data_train
rm temp*

awk '{print $2,$3,$4,$1}' data_train  | sort -n +3 -4  > data.trn
rm sort_data renum_train data_train

mv data.trn  renf90.dat_train 

## replace the input datafile in renf90.par
sed -i 's/renf90.dat/renf90.dat_train/g' renf90.par


 
