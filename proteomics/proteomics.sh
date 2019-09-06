# find transcrips expression
awk -F "\t" 'FNR==NR{a[$1]=$1;next} ($1 in a) {print $1"\t"$10"\t"$11"\t"$12"\t"$13}' proteins.idl ALL.isoform.counts.matrix
