# Author: Michelle Tomaselli
# Date: 08/10/2018
# eQTL Hands-On

After log out from the terminal: 
- docker run -v $PWD:$PWD -w $PWD -it dgarrimar/eqtlmapping
- cd teaching/uvic/AdvBI_2018/data/hands-on/eQTL/
- PATH=$PATH:$PWD/bin
- chmod +x file #If says permission denied 

# Task -> 1
#Download the genotype VCFand the corresponding .tbi index into the subdirectory input/unprocessed/1000g/nameofthefiles
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf{.gz,.gz.tbi} --directory-prefix input/unprocessed/1000g


# Task -> 2
cut -f1 input/unprocessed/geuvadis/geuvadis.metadata.txt | sed '1d' | sort | uniq > tmp/geuvadis.samples.txt 
bcftools view -v snps,indels -m 2 -M 2 -q 0.05:minor -S tmp/geuvadis.samples.txt -Ob input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o tmp/genotypes.chr22.vcf.gz
PATH=$PATH:$PWD/bin
filter.genotype.py -t 10 -g <(zcat tmp/genotypes.chr22.vcf.gz) | bgzip > input/processed/genotypes.chr22.vcf.gz
tabix -p vcf input/processed/genotypes.chr22.vcf.gz

# Q --> 1
-f --> Skip sites where filter column does not contain any of the strings listed in list
-v --> Output variant sites only
-m --> Annotate sites which are present ("+") or absent ("-") in the -a file with a new INFO/TAG flag
-M --> Output sites where REF allele is N="2" in this case
-q --> Break phase set if phasing quality is lower than INT
-S --> Subset of samples to annotate. If the samples are named differently in the target VCF and the -a, --annotations VCF, the name mapping can be given as "src_name dst_name\n", separated by whitespaces, each pair on a separate line.
-Ob --> Basic usage
-d --> snps|indels|both|all|none. Output duplicate records of specified type present in multiple files only once. Requires -a, --allow-overlaps.
-Oz --> Convert results into VCF
-o --> When output consists of a single stream, write it to FILE rather than to standard output, where it is written by default.
-t --> Similar as -r, --regions, but the next position is accessed by streaming the whole VCF/BCF rather than using the tbi/csi index. Both -r and -t options can be applied simultaneously: -r uses the index to jump to a region and -t discards positions which are not in the targets. Unlike -r, targets can be prefixed with "^" to request logical complement. For example, "^X,Y,MT" indicates that sequences X, Y and MT should be skipped. Yet another difference between the two is that -r checks both start and end positions of indels, whereas -t checks start positions only. Note that -t cannot be used in combination with -T.
-g --> Output also gVCF blocks of homozygous REF calls. The parameter INT is the minimum per-sample depth required to include a site in the non-variant block.
-p --> Option can be used to indicate how to handle unphased data: vcf

# Q --> 2 answer: 74656
zcat input/processed/genotypes.chr22.vcf.gz | grep -v "#" | wc -l

# Q --> 3 
bcftools stats input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > tmp/stats.before 
less -S tmp/stats.before  answer: 2504
bcftools stats input/processed/genotypes.chr22.vcf.gz > tmp/stats.after 
less -S tmp/stats.after  answer: 445


# Task -> 3
# Q --> 1
v12
# Q --> 2
GRCh37
# Q --> 3
version 29, 83129
# Q --> 4
PATH=$PATH:$PWD/bin


release=12
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release/gencode.v$release.annotation.gtf.gz
mv gencode.v$release.annotation.gtf.gz input/unprocessed/gencode/gencode.annotation.gtf.gz
zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_type \"protein_coding\"\|gene_type \"lincRNA\"" | gtf2bed.sh > tmp/gencode.annotation.bed

# Q --> 5
For obtain the length we have to rest the 3 and the 2 column, to know which is the TSS position we have to check if:
-> +upstream will be the 2n column
-> -downstream will be the 3rd column

# Q --> 6
BED --> chr1 9 20, because the index of it is: 0-index, final is opened

# Q --> 7
Because you can not read and write at the same time, talking about the same name of file

awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$3-$2,$6}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else{print $1,$3-1,$3,$4,$5,$6}}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
sed -i "s/^chr//" tmp/gencode.annotation.bed
join -1 4 -2 1 -t $'\t' <(sort -k4,4 tmp/gencode.annotation.bed) <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | sort -k1,1) > tmp/joint.tsv
awk '$2==22' tmp/joint.tsv > tmp/joint.chr22.tsv
paste <(awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' tmp/joint.chr22.tsv) <(cut -f1-6 --complement tmp/joint.chr22.tsv) | sort -k1,1V -k2,2n > tmp/joint.chr22.bed
cat <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | head -1 | sed "s/TargetID/#chr\tstart\tend\tgene\tlength\tstrand/") tmp/joint.chr22.bed > tmp/genes.chr22.rpkm.bed
less -S tmp/genes.chr22.rpkm.bed #Check result


# Task -> 4
# Q --> 1
lncRNA "Long noncoding RNAs" are less express at the cells.

# Q --> 2 
The expression gene "Y" (lineal regression) needs normalisation of ß compounds, so when we use them to calculate the t-Student values and we would obtain p-value that can be compared with the initial hypothesis for accept one or the other. 

# Q --> 3
Para estandarizar los valores de muestras diferentes realizaríamos un gráfico donde cada una tendría un color diverso obteniendo que todas ellas siguen el mismo patrón. 

# Q --> 4
Porque cada gen tiene una distribución próxima a la estándar normal. 


normalize.R -i tmp/genes.chr22.rpkm.bed -o tmp/genes.chr22.norm.bed
bgzip tmp/genes.chr22.norm.bed
tabix -p bed tmp/genes.chr22.norm.bed.gz
mv tmp/genes.chr22.norm.bed.gz* input/processed

# Task -> 5
- File Before 
check.norm.R -i tmp/genes.chr22.rpkm.bed -o result/plots/noncheck.norm.pdf

- File After
check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf

# Q --> 1
Se realizan dos gráficos (un gráfico indica la expresión del gen en cada muestra y el otro es un Q-Q plot) por cada archivo (uno contiene datos normalizados y el otro no). Una vez hechos los gráficos se observan diferencias claras entre ellos. En el gráfico que contiene los datos normalizados, observamos que el rango que agrupa los valores que indican la expresión del gen van de 3 a -3 y en Q-Q plot se observa una recta perfecta agrupando todos los valores. En cambio, en los gráficos de los datos no normalizados se observa mucha dispersión de valores ya que no hay un rango acotado y en la Q-Q plot los quantiles teóricos son muy diferentes respecto a los de las muestras. 


# Task -> 6
head -1 input/unprocessed/1000g/1000g.phase3_metadata.txt  | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
head -1 input/unprocessed/geuvadis/geuvadis.metadata.txt | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
# Q --> 1
 ->Respecto el archivo 1000g.phase3_metadata.txt analizaríamos como variables: gender, super_pop ya que son las que menos valores tienen. 
 ->Respecto el archivo geuvadis.metadata.txt analizaríamos como variables: Performer"lab name", Factor Value[population] ya que son algunas de las que menos valores tienen. 


# Task -> 7
QTLtools pca --bed input/processed/genes.chr22.norm.bed.gz --scale --center --out result/expression 
QTLtools pca --vcf input/processed/genotypes.chr22.vcf.gz --scale --center --maf 0.05 --distance 50000 --out result/genotypes
# Q --> 1
--bed : es el formato utilizado para fenotipos
--vcf : es el formato utilizado para genotipos
--scale : escala los valores del fenotipo antes de realizar la PCA
--center : centra los valores del fenotipo antes de realizar la PCA
--out : crea un archivo con los datos datos analizados
--maf : solo considera variantes con una maf"minor allel frequency" respecto al valor dado
--distance : solo considera variantes separadas en la distancia dada 

# Q --> 2
Contiene el resultado de la PCA que incluye una matriz, en un caso respecto la expresión y por el otro respecto el genotipo sirviéndonos para detectar valores atípicos antes de realizar el mapping QTL. 

pcaPlot.R -i result/expression -o result/plots/expression.pca.pdf
pcaPlot.R -i result/genotypes -o result/plots/genotypes.pca.pdf

# Q --> 3
Observamos que en el gráfico del archivo de expresión los valores estan dispersos agrupados en un rango de 10 a -10 y, en cambio, en el gráfico del archivo de genotipos observamos dos agrupaciones de los valores diferenciadas entre ellas y en un rango de 5 a -10.

pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color super_pop --out result/plots/genotypes.pca.super_pop.pdf

# Q --> 4
Damos color a las variables del gráfico y observamos que los genotipos de los dos grupos diferenciados son uno el africano y otro el europeo. 

join -j 1 -t $'\t' <(sort -k1,1 input/unprocessed/1000g/1000g.phase3_metadata.txt) <(cut -f1,20 input/unprocessed/geuvadis/geuvadis.metadata.txt | sort -k1,1 | uniq) > tmp/metadata.txt
sed -i '1s/^/sampleID\tpop\tsuper_pop\tgender\tlab\n/' tmp/metadata.txt
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m tmp/metadata.txt --formula "~ (1|gender) + (1|pop) + (1|lab)" -o result/plots/vp.pdf

# Q --> 5
Como se observa en el documento "vp.pdf" el factor que explica mas porcentaje de variación és la variable laboratorio, seguido de la variable población y muy poco vinculado a la variable sexo. Luego vemos que el resto de variables se agrupan en un mismo grupo el cual no tiene ninguna vinculación relacionada con los datos estudiados. 


# Task -> 8   hidden covariates from expression matrix 
peer.R -i input/processed/genes.chr22.norm.bed.gz -p 10 -o tmp/peer.tsv
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m <(paste tmp/peer.tsv tmp/metadata.txt) -f "~ (1|pop) + (1|lab) + PEER1 + PEER2 + PEER3 + PEER4 + PEER5" -o result/plots/vp.peer.pdf
join -j 1 -t $'\t' tmp/metadata.txt tmp/peer.tsv  | Rscript -e 'write.table(t(read.table(file("stdin", open = "r", blocking = T), h = F)), file = "input/processed/covariates.tsv", quote = F, sep = "\t", col.names = F, row.names = F)'
gzip input/processed/covariates.tsv

# Q --> 1
En el documento "vp.peer.pdf" podemos observar que las variables laboratorio y población siguen teniendo influencia en el porcentaje de variación, así como el porcentaje respecto al sexo desaparece. Además, vemos como surgen nuevas covariables ocultas provenientes de la matriz de expresión explicando más el valor obtenido en los diferentes porcentajes. Observándose que la PEER 1 3 y 2 tienen un valor de porcentaje más alto respecto a las variables lab y pop, seguidas de la PEER 4 y 5. Por último, verificamos que la variable sexo se encuentra en la agrupación "Residuals" junto al resto de variables estudiadas. 


# Task -> 9
QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --nominal 1 --out result/nominals.txt
 -- Description: 
1. The phenotype ID
2. The chromosome ID of the phenotype
3. The start position of the phenotype
4. The end position of the phenotype
5. The strand orientation of the phenotype
6. The total number of variants tested in cis
7. The distance between the phenotype and the tested variant (accounting for strand orientation)
8. The ID of the tested variant
9. The chromosome ID of the variant
10. The start position of the variant
11. The end position of the variant
12. The nominal P-value of association between the variant and the phenotype
13. The corresponding regression slope
14. A binary flag equal to 1 is the variant is the top variant in cis

# Q --> 1
Porque hay snips que se heredan en bloques, "linkage disequilibrium".

pvdist.R -i result/nominals.txt --col 12 -o result/plots/pvdist.pdf

# Q --> 2
En el primer gráfico, respecto a la recta en color rojo podemos observar como la base logarítmica en base 10 de los valores p esperados eje(x), va de un rango entre 0 y 4, y respecto al eje(Y) no muestra mucha diferencia ya que los valores se encuentran entre 0 y 5. En esta última, respecto a la recta de color negro podemos ver como hay un incremento bastante notable de los valores p respecto al eje (Y) que van des del 0 al 25, comenzando con mayor fluidez de valores respecto al final donde hay menos. 
En el segundo gráfico, vemos los p valores respecto a la frecuencia en la que se encuentran en el estudio realizado. Podemos observar, como los valores pv rondando el 0'0 se encuentran en mayor frecuencia. En cambio, podemos observar una disminución considerada respecto a la primera barra y la segunda, y que luego a medida que incrementa el pv van aumentado muy poco cada una de estas.  

plink --ld rs200240198 rs200858043 --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/ld2

# Q --> 3
rs200240198 and rs200858043    (pv = 0.825447)
Observas que los alelos presentes en estos snips son GA/AG, tienen una frecuencia de 0.273034 y 0.726966, respectivamente. También vemos que ambos son valores mayores a los "expectation under LE". Además que aquellos que presentan una frecuencia 0 (AA GG), con "expectation under LE" cobran valores iguales entre ellos. 


# Task -> 10
QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --out result/permutations.txt
 -- permutations.txt Description: 
1. The phenotype ID
2. The chromosome ID of the phenotype
3. The start position of the phenotype
4. The end position of the phenotype
5. The strand orientation of the phenotype
6. The total number of variants tested in cis
7. The distance between the phenotype and the tested variant (accounting for strand orientation)
8. The ID of the top variant
9. The chromosome ID of the top variant
10. The start position of the top variant
11. The end position of the top variant
12. The number of degrees of freedom used to compute the P-values
13. Dummy
14. The first parameter value of the fitted beta distribution
15. The second parameter value of the fitted beta distribution (it also gives the effective number of independent tests in the region)
16. The nominal P-value of association between the phenotype and the top variant in cis
17. The corresponding regression slope
18. The P-value of association adjusted for the number of variants tested in cis given by the direct method (i.e. empirircal P-value)
19. The P-value of association adjusted for the number of variants tested in cis given by the fitted beta distribution. We strongly recommend to use this adjusted P-value in any downstream analysis

for j in $(seq 1 16); do
  echo "cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --chunk $j 16 --out result/permutations_$j.txt"
done | xargs -P4 -n14 QTLtools
cat result/permutations_*.txt > result/permutations.txt; rm result/permutations_*.txt

- root@4ae6aa061761:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# R
..... Info R
-Run: 
p <- read.table("result/permutations.txt")                                                      # Read input file
pdf("result/plots/pv-correlation.pdf",  paper = 'a4r', width = 9, height = 6)                   # Open PDF device
plot(p[, 18], p[, 19], xlab = "pv (perm)", ylab = "pv (beta)")                                  # Plot p-values
abline(0, 1, col = "red")                                                                       # Add red line 1=1
plot(-log10(p[, 18]), -log10(p[, 19]), xlab = "-log10 pv (perm)", ylab = "-log10 pv (beta)")    # Repeat in -log10 space to check the behaviour of the small p-values.
abline(0, 1, col = "red")
dev.off()                                                                                       # Close device
quit("no")   										        # Exit R
Pv-correlation.pdf --> Podemos observar como los pvalues "perm" y los "beta" coinciden considerablemente. Esto sucede porque los valores que corresponden a los diferentes genotipos, se encuentran agrupados en el LD=1, por lo que al realizar la regresión lineal son iguales (primer gráfico obtenido). 

# Task -> 11
mtc.R -n result/nominals.txt -p result/permutations.txt --out result/nominalspeq.txt
mtc.R -n result/nominals.txt -p result/permutations.txt --method bonferroni --out tmp/bonferroni.txt
mtc.R -n result/nominals.txt -p result/permutations.txt --method fdr --out tmp/FDR.txt
mtc.R -n result/nominals.txt -p result/permutations.txt —-method perm-fdr --out result/eqtls.tsv      --> Global empirical P-value threshold = 1.20e-02

# Q --> 1
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# awk '{if ($12 < 0.05) print $12}' result/nominals.txt | wc -l
95206
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# awk '{if ($12 < 0.05) print $12}' tmp/FDR.txt | wc -l
17267
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# awk '{if ($12 < 0.05) print $12}' result/eqtls.tsv | wc -l
12044
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# awk '{if ($12 < 0.05) print $12}' tmp/bonferroni.txt | wc -l
6007

Observando los valores de pv>0'05 podemos determinar que cada técnica tiene un rango diferente de exclusion de los valores partiendo de nominals.txt como archivo a analizar. Siendo mucho más estricta la técnica realizada por el método bonferroni y la menos el método que combina permutacion y FDR.  


# Task -> 12 Have a look to the results obtained 
eQTLviewer.R -i <(head -n 10 result/eqtls.tsv) -g input/processed/genotypes.chr22.vcf.gz -e input/processed/genes.chr22.norm.bed.gz -o result/plots/eQTLs_head.pdf --verbose


# Task -> 13
rsync -av rsync://ftp.ensembl.org/ensembl/pub/grch37/release-86/regulation/homo_sapiens/AnnotatedFeatures.gff.gz input/unprocessed/ensembl
zcat input/unprocessed/ensembl/AnnotatedFeatures.gff.gz | awk 'BEGIN{FS=OFS="\t"}{print $1, $4-1, $5, $9}' | sed -r 's/Name=([^;]+);.*/\1/' | grep -v '^GL' | sort -V > tmp/ERB.bed
for feat in $(cut -f4 tmp/ERB.bed | sort | uniq); do bedtools merge -i <(grep -Fw $feat tmp/ERB.bed) -c 4 -o distinct 
done > input/processed/ERB.collapsed.bed 
sed -i "s/^chr//" input/processed/ERB.collapsed.bed
for feat in $(cut -f4 input/processed/ERB.collapsed.bed | sort | uniq); do QTLtools fenrich --qtl <(sed '1d' result/eqtls.tsv | awk '{print $9, $10-1, $10, $8, $1, "."}') --tss tmp/gencode.annotation.bed  --bed <(grep -Fw $feat input/processed/ERB.collapsed.bed) --out tmp/enrich.txt > /dev/null; echo "$(cat tmp/enrich.txt) $feat" 
done | grep -Fwv inf | grep -Fwv nan > result/enrichments.txt
plot.enrich.R -i result/enrichments.txt -o result/plots/enrich.pdf

# Q --> 1
Los valores mas enriquecidos de odds ratio superiores a 1,5 son H3K36me3(asociados a modificaciones sobre las histonas), seguido de PolII, es decir, pv mas pequeños. 
Se observa que EBF1 presenta mayor variabilidad respecto al resto del rango establecido de odds. En cambio la H3K36me3 presenta mayor enriquecimiento pero su variabilidad es muy pequeña. El resto mas o menos tienen una variabilidad similar, siendo E2F6 mayor. En cambio, el valor de enriquecimiento de las restantes es mayor en MAx. 
Se han observado menos EQTL respecto a los valores esperados. 



# Task -> 14
sed '1d' result/eqtls.tsv | cut -f8 | sort | uniq > tmp/eqtls_snps.tsv
head -100 tmp/eqtls_snps.tsv > tmp/eqtls_snps_100.tsv             -->  Analizaremos los 100 primeros

---------------> A mi manera -------------------
# Q --> 1
 ---> STOP_GAINED:
intron_variant: 45%
downstream_gene_variant: 15%
non_coding_transcript_variant: 14%
upstream_gene_variant: 11%
NMD_transcript_variant: 6%
3_prime_UTR_variant: 3%
non_coding_transcript_exon_variant: 2%
regulatory_region_variant: 1%
synonymous_variant: 1%
Others . . . 
--> ACCEPTOR SPLIPCE:
intron_variant: 53%
upstream_gene_variant: 14%
downstream_gene_variant: 11%
non_coding_transcript_variant: 9%
NMD_transcript_variant: 7%
regulatory_region_variant: 2%
3_prime_UTR_variant: 1%
intergenic_variant: 1%
non_coding_transcript_exon_variant: 1%
Others... 

# Q --> 2
En estos 100 SNP encontramos 1 variante respecto al archivo generado como stop_gained.txt, que tiene como consecuencia: stop_gained.  
En estos 100 SNP encontramos 1 variante respecto al archivo generado como aceptorsplice.txt, que tiene como consecuencia: splice_acceptor_variant. 

# Q --> 3
 ---> STOP_GAINED:  documento result/stop_gained.txt
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# wc -l <(awk '$4=="splice_acceptor_variant"' result/stop_gained.txt )
0 /dev/fd/63

--> ACCEPTOR SPLIPCE:  documentoresult/aceptorsplice.txt
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# wc -l <(awk '$4=="splice_acceptor_variant"' result/aceptorsplice.txt )
3 /dev/fd/63

--------------> a partir del link de Ari ---------------
http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?db=core;tl=315nTgfC4grzij8q-4700617

# Q --> 1
intron_variant: 47%
upstream_gene_variant: 15%
downstream_gene_variant: 14%
non_coding_transcript_variant: 11%
NMD_transcript_variant: 6%
regulatory_region_variant: 2%
intergenic_variant: 2%
non_coding_transcript_exon_variant: 1%
3_prime_UTR_variant: 1%
Others

# Q --> 2
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# wc -l eQTL_HandsOn/result/IMPACTHIGH.txt 
24 eQTL_HandsOn/result/IMPACTHIGH.txt

root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# uniq <(cut -f4 eQTL_HandsOn/result/IMPACTHIGH.txt)
Consequence
splice_acceptor_variant,non_coding_transcript_variant
stop_gained
splice_donor_variant,frameshift_variant
frameshift_variant
stop_gained
splice_acceptor_variant
splice_acceptor_variant,NMD_transcript_variant
splice_acceptor_variant
splice_acceptor_variant,non_coding_transcript_variant
stop_gained

# Q --> 3
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# wc -l <(awk '$4=="splice_acceptor_variant"' eQTL_HandsOn/result/IMPACTHIGH.txt)
3 /dev/fd/63



# Task -> 15
cut -f1 result/eqtls.tsv | sed '1d' | sed 's/\..\+//' | sort | uniq > tmp/egenes.txt
awk '{if($1==22) print $4}' tmp/gencode.annotation.bed | sed 's/\..\+//' | sort | uniq > tmp/bg.txt

# Q --> 1
Process: Response to lipopolysaccharide and response to molecule of bacterial origin. 
Function: Ras guanyl-nucleotide exchange factor activity.
Component: endoplasmic reticulum.



# Task -> 16
grep -Fwf <(cut -f1 result/eqtls.tsv ) result/permutations.txt > tmp/rtc_input
cut -f4,7 input/unprocessed/gwas/gwas.catalog.hg19.bed > tmp/gwas_trait
wget http://jungle.unige.ch/QTLtools_examples/hotspots_b37_hg19.bed --directory-prefix tmp
sed -i 's/^chr//' tmp/hotspots_b37_hg19.bed
QTLtools rtc --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --hotspot tmp/hotspots_b37_hg19.bed --gwas-cis tmp/gwas_trait tmp/rtc_input --out result/rtc.txt

  - rtc.txt Description of the columns in the output file:
1. GWAS variant
2. eQTL variant
3. Gene
4. Gene Group
5. GWAS variant chromosome
6. GWAS variant position
7. GWAS variant rank
8. eQTL variant chromosome
9. eQTL variant position
10. eQTL variant rank
11. Gene chromosome
12. Gene position
13. Distance between variants
14. Distance between GWAS and phenotype
15. GWAS variant's region index
16. eQTL variant's region index
17. Region start
18. Region end
19. Number of variants in the region
20. RTC
21. D prime
22. R squared

# Q --> 1
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# awk '{if ($20 > 0.9) print $20}' result/rtc.txt | wc -l
39

# Q --> 2
awk '{if ($20 > 0.999) print $0}' result/rtc.txt
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# awk '{if ($20 > 0.999) print $1,$2,$20}' result/rtc.txt 
other_variant our_variant RTC
rs909685 rs909685 1
rs9611565 rs4820438 0.99902
rs2234052 rs9607799 0.99902

1 --> rs909685
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# grep rs909685 input/unprocessed/gwas/gwas.catalog.hg19.bed
22	39747670	39747671	rs909685	0	+	Rheumatoid arthritis
22	39747670	39747671	rs909685	0	+	Rheumatoid arthritis (ACPA-positive)
Un caso en el que asociemos la enfermedad con el snip, PubMed.
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# awk '{if ($20 > 0.9999) print $0}' result/rtc.txt
other_variant our_variant phenotype phenotype_group other_variant_chr other_variant_start other_variant_rank our_variant_chr our_variant_start our_variant_rank phenotype_chr phenotype_start distance_between_variants distance_between_other_variant_and_pheno other_variant_region_index our_variant_region_index region_start region_end variant_count_in_region RTC D' r^2
rs909685 rs909685 ENSG00000100321.10 ENSG00000100321.10 22 39747671 0 22 39747671 0 22 39745930 0 1741 65647 65647 39660501 39757500 243 1 1 1 --> vemos de nuestra selección que es 1, cual es el nombre  del gen, columna 3  -->  ENSG00000100321.10
https://www.ncbi.nlm.nih.gov/snp/rs909685 Podemos ver que el gen al cual pertenece este SNP es el mismo obtenido a través del documento rtr.txt
https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000100321%20and%20%22Rheumatoid%20arthritis%22 Es el acceso al cual llegamos a través del enlace anterior yendo a General gene information, Gene ontology podemos observar las características que presentan los estudios realizados en este gen como funciones, procesos y componentes.

# Q --> 3
--> SNP: rs909685
   Consequences: 
 - intron_variant 75%
 - NMD_transcript_variant: 13%
 - rnon_coding_transcript_variant: 13%



# Task -> 17
gene=ENSG00000100321.10
compZscore.R --gene $gene --nominal result/nominals.txt -k 50 --output tmp/$gene.rs_z
plink --r square --snps $(cut -f1 tmp/$gene.rs_z) --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/$gene
CAVIAR -z tmp/$gene.rs_z -l tmp/$gene.ld -o result/$gene

root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# head result/ENSG00000100321.10.log 
inf

# Q --> 1
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# wc -l result/ENSG00000100321.10_set 
4 result/ENSG00000100321.10_set
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# cat result/ENSG00000100321.10_set
rs2069235
rs137689
rs137687
rs137688

# Q --> 2
root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# grep -f result/ENSG00000100321.10_set result/ENSG00000100321.10_post 
SNP_ID	Prob_in_pCausalSet	Causal_Post._Prob.
rs2069235	0.5	1
rs137689	0.227911	0.455823
rs137687	0.150071	0.300142
rs137688	0.122017	0.244035

# Q --> 3
rs2069235 --> 1.96196e-74 0.946899
rs137689 --> 8.85543e-28 -0.596472
rs137687 --> 9.57315e-28 -0.598891
rs137688 --> 9.94976e-28 -0.598014
Haciendo un grep ENSG00000100321.10 result/nominals.txt | less podemos observar que a pesar que el primer valor si que es mayor, el resto respecto a los valores que hay dentro del archivo siguen siendo relevantes sobre la expresión del gen. Correspondientes a : 12. The nominal P-value of association between the variant and the phenotype y 13. The corresponding regression slope. 

# Q --> 4
rs2069235 --> intron_variant: 75%
	      NMD_transcript_variant: 13%
	      non_coding_transcript_variant: 13%
rs137689 --> upstream_gene_variant: 100%
rs137687 --> regulatory_region_variant: 50%
	     intergenic_variant: 50%
rs137688 --> upstream_gene_variant: 100%


# Task -> 18
gene=ENSG00000100321.10
cat <(echo "MarkerName P.value") <(grep $gene result/nominals.txt | cut -d " " -f8,12) > tmp/metal.$gene         SNP --> rs2069235

root@7fb12de33375:/Users/michelletomaselli/teaching/uvic/AdvBI_2018/data/hands-on/eQTL# head -5 tmp/metal.ENSG00000100321.10 
MarkerName P.value
rs57485857 0.0835219
rs3747169 0.736051
rs6001099 0.581104
rs111238998 0.0660701

En el plot metal.ENSG00000100321.10 podemos ver como efectivamente uno de nuestros SNP el cual supusimos que influía en el gen, se encuentra en la parte superior "rs2069235" siendo significativo en nuestro estudio realizado. 



# Task -> 19
mkdir eQTL_HandsOn
mv result/ eQTL_HandsOn
gzip eQTL_HandsOn/result/nominals.txt 
cd eQTL_HandsOn
git init
git add *
git commit -m "eQTL_HandsOn"
git remote add origin https://github.com/michelletomaselli/eQTL_HandsOn.git
git push -u origin master


