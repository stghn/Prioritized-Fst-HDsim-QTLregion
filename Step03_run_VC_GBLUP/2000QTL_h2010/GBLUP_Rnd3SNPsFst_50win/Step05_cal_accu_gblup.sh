#!/bin/bash -i

# Instructions for running this Bash script:
# 1. Make sure you have the necessary permissions to execute this script.
# 2. Open a terminal and navigate to the directory where this script is saved.
# 3. Run the script in the backgorund of server by typing the following command:
#    nohup bash Step05_cal_accu_gblup.sh &
# 4. Follow the on-screen instructions and prompts, if any, during the script execution.
# 5. Ensure you have the required dependencies and data files in the same directory as this script.
# 6. Contact the script author for any issues or questions.

shopt -s expand_aliases

# Extract the first and tenth columns from 'renadd03.ped', sort by the first column, and save the output to 'out1'
awk '{print $1, $10}' renadd03.ped | sort +0 -1 > out1

# Extract the third and fourth columns from 'solutions' where the second column is equal to 3, sort by the first column, and save to 'sol'
awk '{ if ($2==3) print $3,$4}' solutions | sort +0 -1 > sol

# Join 'out1' and 'sol' based on the first column, reorder the columns, sort by the first column, and save to 'file1'
join -1 +1 -2 +1 out1 sol | awk '{print $2,$1,$3}' | sort +0 -1 > file1

# Extract the first and thirteenth columns from 'full.data.geno910', sort by the first column, and save to 'tbv'
awk '{print $1, $13}' full.data.geno910 | sort +0 -1 > tbv

# Join 'file1' and 'tbv' based on the first column and save to 'sol.tbv'
join -1 +1 -2 +1 file1 tbv > sol.tbv

# Remove temporary files 'out1', 'sol', 'file1', and 'tbv'
rm out1 sol file1 tbv

## Sort 'sol.tbv' and save to 'temp1', sort '../train.gen9.ID' and save to 'temp2',
## join based on the first column and save to 'sol.tbv.train', remove temporary files
sort +0 -1 sol.tbv > temp1
sort +0 -1 ../train.gen9.ID > temp2
join -1 1 -2 1 temp1 temp2 > sol.tbv.train
rm temp*

## Sort 'sol.tbv' and save to 'temp1', sort '../test.gen10.ID' and save to 'temp2', 
##join based on the first column and save to 'sol.tbv.test', remove temporary files and 'sol.tbv'
sort +0 -1 sol.tbv > temp1
sort +0 -1 ../test.gen10.ID > temp2
join -1 1 -2 1 temp1 temp2 > sol.tbv.test
rm temp* sol.tbv

# Execute 'corr_ebv_tbv.R' in R version 4.1.3 in batch mode
R4.1.3 CMD BATCH corr_ebv_tbv.R &
