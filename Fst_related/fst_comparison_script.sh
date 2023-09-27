sed -i 's/chr//g' flps_detected_single_chronly_exon_repeats_chr_rpts_combined_compact


##10 10337 0.00129976
##cd /media/workstation/sde3/lokdeep/human_genetics/1000genomes/fst_files/fst_lcr
for  i in `ls *weir_fst.weir.fst_filtered`
do
echo $i
sed '1d' "$i"|sort -k1,1 -k2,2n|awk '{print $1,$2-1,$2,$3}' OFS="\t" > temp
mv temp "$i"
done
########
sed -i 's/chr//1' exon_coordinates
sort -k1,1 -k2,2n exon_coordinates > temp
mv temp exon_coordinates

for f in `ls -1 *fst_filtered`
do
echo $f
n=`echo $f|sed 's/weir_fst.weir.fst_filtered//g'`
bedtools intersect -wa -a "$f" -b exon_coordinates > "$n"exons_fst
done 
##############
###seggregating Fst of LCRs and non-LCRs
###firstly, Fst of LCR regions
for f in `ls *exons_fst`
do
echo $f
n=`echo $f|sed 's/exons_fst//g'`
bedtools intersect -wa -a "$f" -b flps_detected_single_chronly_exon_repeats_chr_rpts_combined_compact > "$n"lcr_fst
done

####Selecting the Fst of non-lcr regions
for f in `ls *exons_fst`
do
echo $f
n=`echo $f|sed 's/exons_fst//g'`
bedtools subtract -A -a "$f" -b flps_detected_single_chronly_exon_repeats_chr_rpts_combined_compact > "$n"non_lcr_fst
done

rm -rf non_lcr_fst_combined
for f in `ls -1 *_non_lcr_fst`
do
echo $f
n=`echo $f|sed 's/_non_lcr_fst//g'`
awk -v n=$n '{print n,"Non-LCR",$4}' OFS="\t" "$f" >> non_lcr_fst_combined
done

rm -rf lcr_fst_combined
for f in `ls -1 *_lcr_fst|grep -v "non_lcr"`
do
echo $f
n=`echo $f|sed 's/_lcr_fst//g'`
awk -v n=$n '{print n,"LCR",$4}' OFS="\t" "$f" >> lcr_fst_combined
done

cat lcr_fst_combined non_lcr_fst_combined > all_fst_combined
