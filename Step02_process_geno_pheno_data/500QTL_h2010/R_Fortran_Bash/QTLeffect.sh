## Skip headers from effect_qtl_001.txt and remove numbers and : before the numbers in coulmn 3 and 4
tail -n +2 output_rep1/effect_qtl_001.txt | sed 's/[0-9]://g' > temp1

## subtract column 3 (allele1 effect) from column 4 (allele2 effect) to get QTL effect
## and print out QTL id, chromosome and QTL effect 
awk 'BEGIN{FS=OFS=" "}{print $1,$2,$3,$4,$3-$4}' temp1 | awk '{print $1,$2,$5}' > temp2

## convert all of QTL effect to the absolute value 
awk '{for (i=1; i<=NF; i++) if ($i < 0) $i = -$i; print }' temp2 > temp3

mv temp3 QTLeffect
rm temp*

## Skip headers from lm_qtl_001.txt
tail -n +2 output_rep1/lm_qtl_001.txt > skip_lm_qtl

## adding position of QTL inside QTL effect file 
awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }' skip_lm_qtl QTLeffect > temp1

awk 'BEGIN{FS=OFS=" "}{$4=$5=" "}{print}' temp1 > temp2
mv temp2 QTLpos
rm temp1

## Skip headers from p1_freq_qtl_001.txt
tail -n +2 output_rep1/p1_freq_qtl_001.txt  > skip_p1_freq_qtl


### adding allele freq of QTL in generation 5
awk '$2==5 {print}' skip_p1_freq_qtl | sed 's/[0-9]://g' | awk '{print $1, $5, $6}' > temp1
mv temp1 QTLfreq


awk 'NR==FNR { a[c=FNR]=$0; next }
     { printf "%-8s\t%s\n", a[FNR], $0 }
     END { for(i=FNR+1;i<=c;i++) print a[i] }' QTLpos QTLfreq > temp1

awk '{print $1,$2,$3,$4,$6,$7}' temp1 > QTLinfo
sed -i '1iIndex chr pos Alpha Freq1 Freq2' QTLinfo
rm temp1 QTLfreq QTLpos QTLeffect skip*


## adding allelic effect of 0 for the marker file which have postion on chromosome
tail -n +2 output_rep1/lm_mrk_001.txt | sed 's/^/0 /g' > temp1

awk 'BEGIN {FS=OFS=" "} {temp=$1; $1=$5; $5=temp} {print}' temp1 > Markerinfo
sed -i '1iIndex chr pos Alpha' Markerinfo
rm temp1

