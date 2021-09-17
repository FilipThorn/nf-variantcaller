#! /usr/bin/bash -l

input=$1
output_var=$2
output_trans=$3

cat $input | grep -v "#" | awk '{if (length($4)>1||length($5)>1){a="INDEL";b=length($4)-length($5);cnt[b]+=1;} else {a="SNP";b="-";} printf("%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, a, b, $4, $5);}' > $output_var

awk '$3 ~ "SNP" {lbl=$5"->"$6; count[lbl]++} END{for(i in count) print count[i], i}' $output_var | sort -gr |column -t > $output_trans