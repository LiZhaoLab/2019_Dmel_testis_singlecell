#identification of de novo transcripts in Dmel testis
#Step 1: Use stringtie to find all transcripts in Drosophila Raleigh line testis and accessory gland sequence files, use stringtie merge to obtain a consensus transcriptome for testis and accessory gland 


for i in Dmel_304_l1.bam Dmel_307_l1.bam Dmel_517_l1.bam 171231_R304_ac.bam 171231_R307_ac.bam 171231_R517_ac.bam Dmel_R360_ac.bam Dmel_R399_ac.bam Dmel_R357_ac.bam Dmel_ed10_ac.bam
do
name=$(basename $i .bam)
stringtie -G /ru-auth/local/home/lzhao/Data_store/witt/Genomes/dmel_r6.15_FB2017_02/gtf/dmel-all-r6.15.gtf -o $name.gtf -A $name.tsv $i
done


stringtie --merge Dmel_304_l1.gtf Dmel_307_l1.gtf Dmel_517_l1.gtf 171231_R304_ac.gtf 171231_R307_ac.gtf 171231_R517_ac.gtf Dmel_R360_ac.gtf Dmel_R399_ac.gtf Dmel_R357_ac.gtf Dmel_ed10_ac.gtf -o dmel_testis_AG_merged.gtf -G ~/Data_store/zhao/Genome/dmel_r6.15_FB2017_02/gtf/dmel-all-r6.15.gtf -l 180601_dmel_testis_AG_merged




#first, make list of every stringtie featurethat doesn't overlap an annotated flybase transcript
comm -13 FBgn.txt <(cat 180601_dmel_testis_AG_merged.gtf| awk '{if ($3=="transcript") print $10, $NF }' | grep -v FB | awk '{print $1}' |sort -u) > nonFBgn.txt

#Our RNA seq data is unstranded, so for single-exon de novo genes, make a plus and minus strand version of each.  This way, the reference is prepared for whatever strand the single-cell data identifies the de novo transcript on.
cat 180601_dmel_testis_AG_merged.gtf | awk '{if ($7==".") print $0}'|awk -F"\t" '{OFS=FS}{ $7="+" ; print   }' | awk '{gsub(/merged./, "mergedplus."); print}' > 180717_dmel_testis_AG_merged_plus
cat 180601_dmel_testis_AG_merged.gtf | awk '{if ($7==".") print $0}'|awk -F"\t" '{OFS=FS}{ $7="-" ; print   }'  | awk '{gsub(/merged./, "mergedminus."); print}' > 180717_dmel_testis_AG_merged_minus
cat 180601_dmel_testis_AG_merged.gtf | awk '{if ($7==".") print $0}'>dmel_testis_AG_merged_unstranded.gtf


#names of unannotated stringtie features with blast hits to segregating and fixed de novo genes from zhao et. al 2014 are in 180731_denovonames2.txt

#make list of non FBgn genes that match zhao et. al 2014 de novo genes
 grep -wFf <(sed 's/"//g' nonFBgn.txt| sed 's/;//g'|sort -u |awk '{print $1}') <(sort -u 180731_denovonames2.txt| awk '{print $1}') |awk '{print $1}'| sed 's/\r//'  > nonFBgn2.txt


#Add stringtie features matching de novo genes to flybase annotation file
rm 180731_dmel_testis_AG_mergedref_notFBgn
cat 180601_dmel_testis_AG_merged_everythingelse |grep -v dmel_testis_AG_merged.14749 |grep -v dmel_testis_AG_merged.14750| grep -v dmel_testis_AG_merged.1475 > 180601_dmel_testis_AG_merged_everythingelse2
for i in `cat nonFBgn2.txt`
do
grep -wF $i dmel_testis_AG_merged_unstranded.gtf| awk -F"\t" '{OFS=FS}{ $7="+" ; print   }' | awk '{gsub(/merged./, "mergedplus."); print}'>> 180731_dmel_testis_AG_mergedref_notFBgn
grep -wF  $i dmel_testis_AG_merged_unstranded.gtf |awk -F"\t" '{OFS=FS}{ $7="-" ; print   }'  | awk '{gsub(/merged./, "mergedminus."); print}' >> 180731_dmel_testis_AG_mergedref_notFBgn
grep -wF $i 180601_dmel_testis_AG_merged_everythingelse2 >> 180731_dmel_testis_AG_mergedref_notFBgn
done



sort -u 180731_dmel_testis_AG_mergedref_notFBgn -o 180731_dmel_testis_AG_mergedref_notFBgn

cat 180731_dmel_testis_AG_mergedref_notFBgn dmel_r6.15.cellranger.gtf |grep -v FBgn0002781 |sort -u  > 180731_dmel_testis_AG_mergedref_notFBgn.gtf
#add plus features, minus features and normal stranded features
rm -rf 180731_dmel_testis_AG_mergedref_FBgn

 
#make reference with combined gtf
../cellranger-2.1.1/cellranger mkref --genome=180731_dmel_testis_AG_mergedref_FBgn --fasta=../Genomes/dmel_r6.15_FB2017_02/fasta/dmel-all-chromosome-r6.15.fasta --genes=180731_dmel_testis_AG_mergedref_notFBgn.gtf

#run cellranger

rm -rf 180731_517_novogene

../cellranger-2.1.1/cellranger count --force-cells=5000  --id=180731_517_novogene --fastqs=180711_novogene/C202SC18061959/raw_data/Dm_R517 --chemistry=SC3Pv2 --transcriptome=180731_dmel_testis_AG_mergedref_FBgn 

echo "finished cellranger"

#Fix bug to add gene short names to TSV file that cellranger produces

cp 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genes.tsv 180731_517_novogene/outs/filtered_gene_bc_matrices/genesbackup.tsv
cat  dmel_r6.15.cellranger.gtf |grep FBgn | awk '{if ($3=="gene") print $10, $12}' | sed 's/"//g' | sed 's/;//g' |sort -u >180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genenames.txt
Rscript merged_and_ordered.R 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genes.tsv 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genenames.txt 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genes2.tsv
i=$(wc -l 180731_517_novogene/outs/filtered_gene_bc_matrices/genesbackup.tsv|awk '{print $1}' )
tail -n $i 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genes2.tsv | sed 's/"//g' | awk -v OFS='\t' '{if ($4 != "NA") print $2, $4; else if ($4 == "NA") print $2, $2 }' > 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genes3.tsv
cp 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genes3.tsv 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genes.tsv
rm 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genenames.txt
rm 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genes2.tsv
rm 180731_517_novogene/outs/filtered_gene_bc_matrices/180731_dmel_testis_AG_mergedref_FBgn/genes3.tsv

echo "finished clean up"


Rscript seurmo.r

echo "finished seurat and monocle"
