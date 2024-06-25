#we downloaded individual files from https://ptmd.biocuckoo.org/download.php
#and later unzipped and moved to the same directory
mkdir PTMD
##we have unzipped the folders, and then moved them to the same folder named PTMD
mv PTM-Disease\ association Protein\ Information PTM\ Sites PTMD

##We are going to merge two files
##cd $PWD/PTMD/Protein Information
mkdir ../for_flps_ptmd
awk -F"\t" 'NF==7{print $0}' OFS=\t annotation.txt| cut -f1,4,7 > ../for_flps_ptmd/annotation_filtered.txt
cd ../PTM-Disease\ association/
grep Homo data.txt |awk -F"\t" 'NF==9{print $0}' > ../for_flps_ptmd/homo_disease.txt
cd ../for_flps_ptmd/
##removing the header
sed -i '1d' annotation_filtered.txt
##749 annotation_filtered.txt
##making sure that the annotation file has three fields
awk 'NF==3{print $0}' annotation_filtered.txt > temp.txt
mv temp.txt annotation_filtered.txt
##726 annotation_filtered.txt
dos2unix homo_disease.txt 
dos2unix annotation_filtered.txt
rm -rf compiled.txt
while read j
do
unip_id=`echo $j|awk '{print $1}'`
p_id=`echo $j|awk '{print $2}'`
sequence=`echo $j|awk '{print $3}'`
echo $p_id
grep "^$unip_id\b" homo_disease.txt > temp.txt
awk -v p=$p_id -v s=$sequence -F"\t" -v OFS="\t" '{print $0,p,s}' temp.txt >> compiled.txt
done < annotation_filtered.txt
rm temp.txt
##keeping a temp.txt is a bad habit!
##1651 compiled.txt

cut -f2,11 compiled.txt |sort -u|sed -e 's/^/>/g' -e 's/\t/\n/g' > protein_fasta.fa
CompositionMaker -n protein_fasta.fa > composition_out
fLPS2 -c protein_fasta.fa.COMPOSITION protein_fasta.fa > flps_protein.out
##4830 flps_protein.out
sed -i 's/{//g' flps_protein.out
sed -i 's/}//g' flps_protein.out
awk -v OFS="\t" '((($7/(($6-$5)+1)))>0.5) && (length($9) < 5) && ($7 > 3) && ($9!~/X/){print $1,$2,$5,$6,$7,$9}' flps_protein.out >  filtered_flps
##235 filtered_flps

##for bedtools intersect, we gonna shuffle the columns in the right format
awk -v OFS="\t" '{print $1,$3,$4,$2,$6,$7}' filtered_flps > flps.bed
sed -i 's/\t$//g' flps.bed
##removing the lines which do not contain PTM site number
awk -F"\t" '$6 && $7{print $0}' compiled.txt > temp.txt
mv temp.txt compiled.txt
##1131 compiled.txt
awk -F"\t" '$7~/^[0-9]+$/{print $0}' compiled.txt > temp.txt
mv temp.txt compiled.txt
##1111 compiled.txt
awk -v OFS="\t" -F'\t' '{print $2,$7,$7,$1,$3,$4,$5,$6,$11}' compiled.txt > ptmd.bed

bedtools intersect -wa -wb -a flps.bed -b ptmd.bed > flps_ptmd_intersect.out
##32 flps_ptmd_intersect.out

##We gonna split the output file into some files for plotting
##Firstly, to plot the gene
cut -f1,4 flps_ptmd_intersect.out|sort -u > genes_length
#secondly, LCRs
cut -f1,2,3,4,5 flps_ptmd_intersect.out|sort -u > lcrs_info
##TWIST1 shows two LCRs at the same coordinates. We are manually removing the smaller LCR which is completely inside another LCR
##Thirdly, the PTM sites
cut -f4,6,7,11,13 flps_ptmd_intersect.out|sed 's/ /_/g' > ptms_info
