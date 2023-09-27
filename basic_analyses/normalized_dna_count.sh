while read j
do
aarpt=`echo $j|awk '{print $1}'`
dnarpt=`echo $j|awk '{print $2}'`
olco=`echo $j|awk '{print $3}'`
dnrptcount=`echo $dnarpt|awk -F, '{print NF}'`
for i in `echo $dnarpt|sed 's/,/ /g'`
do
echo $aarpt $i $olco $dnrptcount|awk '{print $1,$2,$3/$4}' OFS="\t" >> normalized_dna_count.txt
done
done < first_aggregate.txt
