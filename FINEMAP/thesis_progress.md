### FINEMAP
<!-- toc -->

- [April 16th Yep, still annotating](#April-16th-Yep-still-annotating)
- [April 12th still annotating](#April-12th-still-annotating)
- [April 5th annotating](#April-5th-annotating)
  * [use self-annotated promoter, chipseeker output is wrong](#use-self-annotated-promoter-chipseeker-output-is-wrong)
- [March 22nd](#March-22nd)
  * [run1 all with 1Mb sep except for SCZ](#run1-all-with-1Mb-sep-except-for-SCZ)
  * [AD PASS](#AD-PASS)
  * [ASD PASS](#ASD-PASS)
  * [ADHD PASS](#ADHD-PASS)
  * [BP PASS](#BP-PASS)
  * [SCZ run1 3Mb, LD>0.6](#SCZ-run1-3Mb-LD06)
- [March 21st, at New York Public library :D](#March-21st-at-New-York-Public-library-D)
  * [MAF ALL BASED ON 1KG PHASE3 EUR SUBSET](#MAF-ALL-BASED-ON-1KG-PHASE3-EUR-SUBSET)
  * [switch MDD source to MDD2018_ex23andMe -- same paper, larger dataset](#switch-MDD-source-to-MDD2018_ex23andMe----same-paper-larger-dataset)
- [March 15th -- Working toward ddl :P](#March-15th----Working-toward-ddl-P)
  * [aquiring data](#aquiring-data)
  * [Sources of summary statistics](#Sources-of-summary-statistics)
- [Feb 22](#Feb-22)
- [Original Notes](#Original-Notes)
  * [12.3](#123)
  * [12.8 AD](#128-AD)
  * [12.19 BP](#1219-BP)
  * [12.19 ADHD](#1219-ADHD)
  * [12.19 MDD](#1219-MDD)
  * [12.19 ASD](#1219-ASD)
- [Random Notes](#Random-Notes)

<!-- tocstop -->
#### April 16th Yep, still annotating
* SCZ TCF4
* ASD -- too few coding variants to do GO :P
* SCZ chr6/chr2 ld did not finish calculting
* ADHD
no GO enrichment at all
* AD

``````bash
for i in *.snp; do sed 1d $i | cat >> AD_all.snp; done
``````
* BP -- metabolism??????????
* Segmentation fault -- delete for now, will need to address in the future
* AD chr6 block3 what is wrong with it.... delete for now, need to address in the future
#### April 12th still annotating
* add in SNPs not avalible for fine-mapping
```bash
cat single* double* | awk '{print $2"\t"$1"\t"$3}'> AD_noteli.set
```
#### April 5th annotating
##### use self-annotated promoter, chipseeker output is wrong
* S1. prepare chipseeker input with python
* S2 liftover by R
* S3 chipseeker -- UpSet plot
* S4 data cleaning - bar plot with python
* S5 prepare peak file with python
* S6 intersect with bedtools
* S7 annotation with chromatin accessability and EDA with python
* 

#### March 22nd
##### run1 all with 1Mb sep except for SCZ
##### AD PASS
##### ASD PASS
##### ADHD PASS
##### BP PASS
##### SCZ run1 3Mb, LD>0.6
* Killed list: 1, 3, 4, 8, 10, 11

#### March 21st, at New York Public library :D 
##### MAF ALL BASED ON 1KG PHASE3 EUR SUBSET
* Re-calc MAF with european vcf data
* Re-calc block with 3MB size
* WTF why are there duplicates in lst files ?????????????????
##### switch MDD source to MDD2018_ex23andMe -- same paper, larger dataset

#### March 15th -- Working toward ddl :P
##### aquiring data
* download 1KG p3
* extract EUR
``````bash
FILE=ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
bcftools view -S EUR.id -O z -o ALL.chr$CHR.phase3.EUR.vcf.gz $FILE
``````
##### Sources of summary statistics
* AD
Genome-wide meta-analysis identifies new loci and functional pathways influencing Alzheimer’s disease risk
**imputation: 1KG Phase unknow
Cohort: mostly Europe**
[GWAS Summary Statistics \| CTG](https://ctg.cncr.nl/software/summary_statistics)
* MDD
http://www.med.unc.edu/pgc/files/resultfiles/pgc-mdd-2018-readme-v.3
**imputation: 1KG Phase1
Cohort: Europe** 
* ADHD
Discovery of the first genome-wide significant risk loci for attention deficit/hyperactivity disorder
**imputation: 1KG Phase3
Cohort: mostly Europe**

* BP & SCZ
Genomic Dissection of Bipolar Disorder and Schizophrenia, Including 28 Subphenotypes
**imputation: 1KG Phase3
Cohort: Europe**
* ASD 
Identification of common genetic risk variants for autism spectrum disorder
**imputation: 1KG Phase3
Cohort: Europe** 

#### Feb 22
* update maf2finemap -- running one chromosome a time for parallization
* What to do with blocks with 1 | 2 SNPs? Does not make sense to throw them away -- maybe filter again with more stringent threshold?


#### Original Notes

##### 12.3
* AD new sumstats
cutoff : 5e-6 : 4360 remained.

* ASD (2017 PGC)
*cutoff: 1e-5 : 4028 / 9106201

##### 12.8 AD

**ERROR IN FINEMAP -- NEED TO FIX IT LATER (SEGMENTATION FAULT)**

* split-by-chr.sh
``````bash
* for i in `seq 1 22`; do awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"12"\t"$13}' AD_final_hg38_Jansen_chr$i.bed > AD_final_hg38_Jansen_chr$i.sst; done
``````
* calc-freq.sh
``````bash
awk -F "\t" '{if ($10 <= 5e-6) {print}}' AD_Jansen_hg38.bed > AD_Jansen_hg38_filtered.bed
* for i in `seq 1 22`; do grep -P "^chr$i\t" AD_Jansen_hg38_filtered.bed | awk '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' > AD_filtered_final_hg38_Jansen_chr$i.sst; done
``````
* s1_chr_sst_to_block_sst.sh
* seqnames        start   end     width   strand  A1      A2      SNP     Z       P       MAF     BETA    SE
*  s2_block_to-zfile.sh

##### 12.19 BP
BP

* CHR     SNP     BP      A1      A2      FRQ_A_20129     FRQ_U_21524     INFO    OR      SE      P       Direction       HetPVa
1. reformat
``````bash
sed 1d BDvsCONT.sumstats | awk '{print "chr"$1"\t"$2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$9"\t"$10"\t"$11}' > BD_hg19.bed
``````
3. liftover -> BD_hg38.bed
4. split-by-chr.sh
CHR   BP  SNP   A1   A2    OR    SE   P
4. calc MAF
``````bash
for i in `seq 1 22`; do awk '{print $3}' BD_hg38_chr$i.sst > BD_hg38_chr$i.lst; done
``````
* calc-freq.sh
5.  filter SST
``````bash
awk -F "\t" '{if ($11 <= 5e-6){print}}' BD_hg38.bed > BD_hg38_filtered.bed
``````
``````bash
for i in `seq 1 22`; do grep -P "^chr$i\t" BD_hg38_filtered.bed | awk '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' >  BD_hg38_filtered_chr$i.sst; done
```
```
6. 
* s1_chr_sst_to_block_sst.sh
* s2_block_to-zfile.sh
* s3_gen_ld.sh
``````

##### 12.19 ADHD
ADHD
**!!! SOMETHING WRONG WITH CERTAIN SNPS**
**chr2    110298586       rs151147303**    
**chr2    110298606       rs140127414**   


* CHR     SNP     BP      A1      A2      INFO    OR      SE      P
``````bash
sed 1d adhd_eur_jun2017 | awk '{print "chr"$1"\t"$3-1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' > adhd_hg19.bed
``````
* hg38
seqnames        start   end     width   strand  SNP     A1      A2      OR      SE      P
``````bash
for i in `seq 1 22`; do grep -P "^chr$i\t" adhd_hg38_filtered.bed | awk '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' >  adhd_filtered_hg38_chr$i.sst; done
``````

##### 12.19 MDD
MDD

seqnames        start   end     width   strand  SNP     A1      A2      F1      F2      INFO    OR      SE      P       REPORT

``````bash
for i in `seq 1 22`; do grep -P "^chr$i\t" MDD_hg38.bed | awk '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$12"\t"$13"\t"$14}' > MDD_hg38_chr$i.sst; done``````
*  awk -F "\t" '{if($14 <= 5e-6){print}}' MDD_hg38.bed > MDD_hg38_filtered.bed
``````bash
for i in `seq 1 22`; do grep -P "^chr$i\t" MDD_hg38_filtered.bed | awk '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$12"\t"$13"\t"$14}' > MDD_hg38_filtered_chr$i.sst; done
``````
* calc-freq
* s1_chr_sst_to_block_sst.sh

``````bash
for i in *-*.set; do sed 1d $i | cat >> multi_snps.set; done
sed -i '1 i\BLOCK   SNP     PROB    BP      CHR     A1      A2      OR      SE      P' multi_snps.set
``````

##### 12.19 ASD

* Common risk variants identified in autism spectrum disorder’.
https://www.biorxiv.org/content/early/2017/11/27/224774 (2017)
* 18,381 ASD cases and 27,969 controls

Problematic blocks
* chr15_block_1
---

#### Random Notes
? stepwise conditional analysis

* Imputation: combine date generated in multiple different platforms
* Meta analysis (required summary statistics at each variant)
  * Benefit: increase sample size/statistical power
  
* Phasing: statistical estimation of haplotypes from genotype data.
  * Necessary for allele imputation from reference
  * Most widely used algorithms now: HMM