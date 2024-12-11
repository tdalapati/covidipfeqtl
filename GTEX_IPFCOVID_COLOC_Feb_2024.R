library(coloc)
library(remotes)
library(devtools)
require(data.table)
require(coloc)
require(cowplot)
require(ggrepel)
require(dplyr)
require(tidyr)
library(tidyverse)
require(LDlinkR)
library(locuszoomr)
require(ensembldb)
library(EnsDb.Hsapiens.v86)

setwd("~/Documents/")
#setwd("~/Documents/KoLab/Computation/IPFCOVID")
ipf = fread("meta_gwas_5way_summary_stats_HG38.txt")
covid = fread("COVID19_HGI_A2_ALL_20220403.10k copy.txt")

rs1631350 <- subset(ipf, rsid =="rs1631350")
rs2897075 <- subset(ipf, rsid =="rs2897075")

#gtex, regression in reference to alternative allele relative to reference allele, effect alelle

#TRIM4 LUNG IPF
lung = fread("gtex_eqtl_lung_allpairs_chr7_HG38.txt")
trim4 <- subset(lung, gene_id == 'ENSG00000146833.15')
trim4 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> trim4
trim4_exact <- subset(trim4, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
trim4$chr = gsub("chr", "", trim4$chr) %>% as.integer()
trim4$bp <- as.numeric(trim4$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,trim4, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.74e-77  3.06e-62  8.96e-16  1.00e+00  9.30e-16 
#[1] "PP abf for shared variant: 9.3e-14%"


#ZKSCAN1 in LUNG IPF
lung = fread("gtex_eqtl_lung_allpairs_chr7_HG38.txt")
zkscan1 <- subset(lung, gene_id == 'ENSG00000106261.16')
zkscan1 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> zkscan1
zkscan1_exact <- subset(zkscan1, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
zkscan1$chr = gsub("chr", "", zkscan1$chr) %>% as.integer()
zkscan1$bp <- as.numeric(zkscan1$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,zkscan1, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#5.28e-16  5.89e-01  9.24e-17  1.03e-01  3.08e-01 
#[1] "PP abf for shared variant: 30.8%"

#COPS6 in LUNG IPF
lung = fread("gtex_eqtl_lung_allpairs_chr7_HG38.txt")
COPS6 <- subset(lung, gene_id == 'ENSG00000168090.9')
COPS6 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> COPS6
COPS6_exact <- subset(COPS6, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
COPS6$chr = gsub("chr", "", COPS6$chr) %>% as.integer()
COPS6$bp <- as.numeric(COPS6$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,COPS6, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.60e-16  8.48e-01  8.71e-17  9.72e-02  5.50e-02 
#[1] "PP abf for shared variant: 5.5%"

#TRIM4 in LUNG COVID
subcovid7 <- subset(covid, `#CHR`==7)
mysnp <- subset(subcovid7, rsid=="rs2897075")

mergedcovidgtex7 = merge(subcovid7,trim4, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#8.37e-42  4.44e-38  1.88e-04  1.00e+00  1.94e-04 
#[1] "PP abf for shared variant: 0.0194%"

#ZKSCAN1 in LUNG COVID
subcovid7 <- subset(covid, `#CHR`==7)
mysnp <- subset(subcovid7, rsid=="rs2897075")

mergedcovidgtex7 = merge(subcovid7,zkscan1, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.23e-04  6.54e-01  6.48e-06  3.41e-02  3.12e-01 
#[1] "PP abf for shared variant: 31.2%

#COPS6 in LUNG COVID
subcovid7 <- subset(covid, `#CHR`==7)
mysnp <- subset(subcovid7, rsid=="rs2897075")

mergedcovidgtex7 = merge(subcovid7,COPS6, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.75e-04  9.30e-01  2.50e-06  1.32e-02  5.66e-02 
#[1] "PP abf for shared variant: 5.66%"

rm(lung)

#WHOLE BLOOD
wb = fread("gtex_wholeblood_eqtl_allpairs_chr7_HG38.txt")
trim4 <- subset(wb, gene_id == 'ENSG00000146833.15')
trim4 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> trim4
trim4_exact <- subset(trim4, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
trim4$chr = gsub("chr", "", trim4$chr) %>% as.integer()
trim4$bp <- as.numeric(trim4$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,trim4, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#9.95e-55  1.11e-39  8.96e-16  1.00e+00  9.08e-16 
#[1] "PP abf for shared variant: 9.08e-14%"

#ZKSCAN1 in WB IPF
zkscan1 <- subset(wb, gene_id == 'ENSG00000106261.16')
zkscan1 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> zkscan1
zkscan1_exact <- subset(zkscan1, variant_id == 'chr7_100032719_C_T_b38')
zkscan1$chr = gsub("chr", "", zkscan1$chr) %>% as.integer()
zkscan1$bp <- as.numeric(zkscan1$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,zkscan1, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.23e-16  8.07e-01  1.36e-16  1.52e-01  4.17e-02 
#[1] "PP abf for shared variant: 4.17%"

#COPS6 in WB IPF
COPS6 <- subset(wb, gene_id == 'ENSG00000168090.9')
COPS6 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> COPS6
COPS6_exact <- subset(COPS6, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
COPS6$chr = gsub("chr", "", COPS6$chr) %>% as.integer()
COPS6$bp <- as.numeric(COPS6$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,COPS6, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#5.27e-22  5.88e-07  8.96e-16  1.00e+00  6.16e-08 
#[1] "PP abf for shared variant: 6.16e-06%"


#TRIM4 in WB COVID
subcovid7 <- subset(covid, `#CHR`==7)
mysnp <- subset(subcovid7, rsid=="rs2897075")

mergedcovidgtex7 = merge(subcovid7,trim4, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.93e-30  1.03e-26  1.88e-04  9.99e-01  5.58e-04 
#[1] "PP abf for shared variant: 0.0558%"

#ZKSCAN1 in WB COVID
mergedcovidgtex7 = merge(subcovid7,zkscan1, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.77e-04  9.43e-01  1.69e-06  8.94e-03  4.83e-02 
#[1] "PP abf for shared variant: 4.83%"

#COPS6 in WB COVID
mergedcovidgtex7 = merge(subcovid7,COPS6, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#6.52e-10  3.46e-06  1.88e-04  1.00e+00  1.05e-04 
#[1] "PP abf for shared variant: 0.0105%"

rm(wb)

#FIBROBLAST
fb = fread("gtex_fibroblasts_eqtl_allpairs_chr7_HG38.txt")
trim4 <- subset(fb, gene_id == 'ENSG00000146833.15')
trim4 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> trim4
trim4_exact <- subset(trim4, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
trim4$chr = gsub("chr", "", trim4$chr) %>% as.integer()
trim4$bp <- as.numeric(trim4$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,trim4, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.01e-91  2.25e-76  8.96e-16  1.00e+00  8.99e-16 
#[1] "PP abf for shared variant: 8.99e-14%"

#ZKSCAN1 in FB IPF
zkscan1 <- subset(fb, gene_id == 'ENSG00000106261.16')
zkscan1 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> zkscan1
zkscan1_exact <- subset(zkscan1, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
zkscan1$chr = gsub("chr", "", zkscan1$chr) %>% as.integer()
zkscan1$bp <- as.numeric(zkscan1$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,zkscan1, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#8.09e-30  9.02e-15  8.96e-16  1.00e+00  2.70e-08 
#[1] "PP abf for shared variant: 2.7e-06%"

#COPS6 in fb IPF
COPS6 <- subset(fb, gene_id == 'ENSG00000168090.9')
COPS6 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> COPS6
COPS6_exact <- subset(COPS6, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
COPS6$chr = gsub("chr", "", COPS6$chr) %>% as.integer()
COPS6$bp <- as.numeric(COPS6$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,COPS6, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#6.18e-16  6.90e-01  9.38e-17  1.04e-01  2.05e-01 
#[1] "PP abf for shared variant: 20.5%"

#TRIM4 in FB COVID
subcovid7 <- subset(covid, `#CHR`==7)
mysnp <- subset(subcovid7, rsid=="rs2897075")

mergedcovidgtex7 = merge(subcovid7,trim4, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.11e-48  1.11e-44  1.89e-04  9.99e-01  1.23e-03 
#[1] "PP abf for shared variant: 0.123%"  

#ZKSCAN1 in FB COVID
mergedcovidgtex7 = merge(subcovid7,zkscan1, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.36e-18  7.18e-15  1.89e-04  9.99e-01  4.32e-04 
#[1] "PP abf for shared variant: 0.0432%"

#COPS6 in FB COVID
mergedcovidgtex7 = merge(subcovid7,COPS6, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.45e-04  7.66e-01  5.51e-06  2.89e-02  2.05e-01 
#[1] "PP abf for shared variant: 20.5%"

rm(fb)

#LCL
lcl = fread("gtex_lcl_eqtl_allpairs_chr7_HG38.txt")
trim4 <- subset(lcl, gene_id == 'ENSG00000146833.15')
trim4 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> trim4
trim4_exact <- subset(trim4, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
trim4$chr = gsub("chr", "", trim4$chr) %>% as.integer()
trim4$bp <- as.numeric(trim4$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,trim4, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.88e-30  2.10e-15  8.96e-16  1.00e+00  1.95e-14 
#[1] "PP abf for shared variant: 1.95e-12%"

#ZKSCAN1 in LCL IPF
zkscan1 <- subset(lcl, gene_id == 'ENSG00000106261.16')
zkscan1 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> zkscan1
zkscan1_exact <- subset(zkscan1, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
zkscan1$chr = gsub("chr", "", zkscan1$chr) %>% as.integer()
zkscan1$bp <- as.numeric(zkscan1$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,zkscan1, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.42e-16  8.28e-01  9.78e-17  1.09e-01  6.28e-02 
#[1] "PP abf for shared variant: 6.28%"

#COPS6 in LCL IPF
COPS6 <- subset(lcl, gene_id == 'ENSG00000168090.9')
COPS6 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> COPS6
COPS6_exact <- subset(COPS6, variant_id == 'chr7_100032719_C_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
COPS6$chr = gsub("chr", "", COPS6$chr) %>% as.integer()
COPS6$bp <- as.numeric(COPS6$bp)

subipf7 <- subset(ipf, chromosome==7)
mysnp <- subset(subipf7, rsid=="rs2897075")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex7 = merge(subipf7,COPS6, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedipfgtex7$effect_allele != mergedipfgtex7$alt)
print(Idx)
mergedipfgtex7[Idx,]$slope = - mergedipfgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

names(mergedipfgtex7)[names(mergedipfgtex7) == 'start'] <- 'position'
names(mergedipfgtex7)[names(mergedipfgtex7) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex7, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex7, position < end)
tmp2

mergedipfgtex7[1,]
tmp3 = subset(mergedipfgtex7, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex7$p), y=-log10(mergedipfgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#6.94e-16  7.74e-01  1.57e-16  1.75e-01  5.05e-02 
#[1] "PP abf for shared variant: 5.05%"

#TRIM4 in LCL COVID
subcovid7 <- subset(covid, `#CHR`==7)
mysnp <- subset(subcovid7, rsid=="rs2897075")

mergedcovidgtex7 = merge(subcovid7,trim4, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#9.85e-11  5.34e-07  1.84e-04  9.99e-01  1.08e-03 
#[1] "PP abf for shared variant: 0.108%"

#ZKSCAN1 in LCL COVID
mergedcovidgtex7 = merge(subcovid7,zkscan1, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.71e-04  9.28e-01  1.76e-06  9.46e-03  6.27e-02 
#[1] "PP abf for shared variant: 6.27%"

#COPS6 in LCL COVID
mergedcovidgtex7 = merge(subcovid7,COPS6, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex7, rsid=="rs2897075")
#othersnp <- subset(mergedipfgtex, rsid =="rs2525542")

Idx = which(mergedcovidgtex7$alt != mergedipfgtex7$ALT)
print(Idx)
mergedcovidgtex7[Idx,]$slope = - mergedcovidgtex7[Idx,]$slope

start = 100032719 - 500000
end = 100032719 + 500000

tmp1 = subset(mergedcovidgtex7, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex7, POS < end)
tmp2

mergedcovidgtex7[1,]
tmp3 = subset(mergedcovidgtex7, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex7$all_inv_var_meta_p), y=-log10(mergedcovidgtex7$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs2897075')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]


myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.70e-04  9.22e-01  2.78e-06  1.50e-02  6.26e-02 
#[1] "PP abf for shared variant: 6.26%"

#MUC5B
#Lung
lung = fread("gtex_eqtl_lung_allpairs_chr11_HG38.txt")
muc5b <- subset(lung, gene_id == 'ENSG00000117983.17')
muc5b %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> muc5b
muc5b_exact <- subset(muc5b, variant_id == 'chr11_1219991_G_T_b38')
#subgtex <- subset(lung, chr=="chr7")
#subgtex$variantloc <- subgtex$variant_id
muc5b$chr = gsub("chr", "", muc5b$chr) %>% as.integer()
muc5b$bp <- as.numeric(muc5b$bp)

subipf11 <- subset(ipf, chromosome==11)
subipf11[2456, "p"] <- "1.0e-203"
mysnp <- subset(subipf11, rsid=="rs35705950")
#othersnp <- subset(subipf7, rsid=="rs2525542")

mergedipfgtex11 = merge(subipf11,muc5b, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex11, rsid=="rs35705950")

Idx = which(mergedipfgtex11$effect_allele != mergedipfgtex11$alt)
print(Idx)
mergedipfgtex11[Idx,]$slope = - mergedipfgtex11[Idx,]$slope

start = 1219991 - 500000
end = 1219991 + 500000

names(mergedipfgtex11)[names(mergedipfgtex11) == 'start'] <- 'position'
names(mergedipfgtex11)[names(mergedipfgtex11) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex11, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex11, position < end)
tmp2

mergedipfgtex11[1,]
tmp3 = subset(mergedipfgtex11, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex11$p), y=-log10(mergedipfgtex11$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs35705950')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.59e-201  6.10e-10 1.24e-194  4.12e-10  1.00e+00 
#[1] "PP abf for shared variant: 100%"

#myres <- coloc.abf(dataset1=list(pvalues = tmp4$pval_nominal, N= 515,
#                                 type="quant", snp=tmp4$rsid),
#                   dataset2=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
#                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

subcovid11 <- subset(covid, `#CHR`==11)
mysnp <- subset(subcovid11, rsid=="rs35705950")

mergedcovidgtex11 = merge(subcovid11,muc5b, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex11, rsid=="rs35705950")

Idx = which(mergedcovidgtex11$alt != mergedipfgtex11$ALT)
print(Idx)
mergedcovidgtex11[Idx,]$slope = - mergedcovidgtex11[Idx,]$slope

start = 1219991 - 500000
end = 1219991 + 500000

tmp1 = subset(mergedcovidgtex11, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex11, POS < end)
tmp2

mergedcovidgtex11[1,]
tmp3 = subset(mergedcovidgtex11, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex11$all_inv_var_meta_p), y=-log10(mergedcovidgtex11$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs35705950')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.34e-16  6.80e-09  4.91e-11  5.26e-07  1.00e+00 
#[1] "PP abf for shared variant: 100%"


lcl = fread("gtex_lcl_eqtl_allpairs_chr11_HG38.txt")
muc5b <- subset(lcl, gene_id == 'ENSG00000117983.17')

#MUC5B - pvalue of gtex vs pvalue of IPF GWAS 
mergedipfgtex11$pval_nominal_transformed <- as.numeric(-log10(mergedipfgtex11$pval_nominal))
mergedipfgtex11$p_transformed <- as.numeric(-log10(mergedipfgtex11$p))
mergedipfgtex11$pc <- predict(prcomp(~pval_nominal_transformed+p_transformed, mergedipfgtex11))[,1]
mysnp <- subset(mergedipfgtex11, rsid == 'rs35705950')
p = ggplot(mergedipfgtex11, aes(p_transformed, pval_nominal_transformed, color = pc)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  theme_bw() +
  scale_color_gradient(low = "blue", high = "orangered")+
  xlab("-log10(P-value, IPF GWAS)")+
  ylab("-log10(P-value, GTEx Lung)")
p
p+ geom_text_repel(
  data = mysnp,
  aes(label = rsid),
  nudge_x = 0,
  nudge_y = -2,
  size = 5,
  min.segment.length = 0,
  color = "black")

ggsave("rs35705950_IPF_GTEx_pvalue_comparison.pdf", width = 6, height = 6)

#ATP11A
#Lung
lung = fread("gtex_eqtl_lung_allpairs_chr13_HG38.txt")
atp11 <- subset(lung, gene_id == 'ENSG00000068650.18')
atp11 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> atp11
atp11_exact <- subset(atp11, variant_id == 'chr13_112881427_C_T_b38')
atp11$chr = gsub("chr", "", atp11$chr) %>% as.integer()
atp11$bp <- as.numeric(atp11$bp)

subipf13 <- subset(ipf, chromosome==13)
mysnp <- subset(subipf13, rsid=="rs12585036")

mergedipfgtex13 = merge(subipf13,atp11, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex13, rsid=="rs12585036")

Idx = which(mergedipfgtex13$effect_allele != mergedipfgtex13$alt)
print(Idx)
mergedipfgtex13[Idx,]$slope = - mergedipfgtex13[Idx,]$slope

start = 112881427 - 500000
end = 112881427 + 500000

names(mergedipfgtex13)[names(mergedipfgtex13) == 'start'] <- 'position'
names(mergedipfgtex13)[names(mergedipfgtex13) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex13, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex13, position < end)
tmp2

mergedipfgtex13[1,]
tmp3 = subset(mergedipfgtex13, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex13$p), y=-log10(mergedipfgtex13$pval_nominal))
mysnp <- subset(tmp3, rsid == 'rs12585036')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

plot(x = mergedipfgtex13$position, y=-log10(mergedipfgtex13$pval_nominal))
mysnp <- subset(mergedipfgtex13, rsid == 'rs12585036')
points(x=(mysnp$position), y=-log10(mysnp$pval_nominal),col='red')


tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#9.97e-09  5.26e-01  2.41e-09  1.27e-01  3.47e-01 
#[1] "PP abf for shared variant: 34.7%"

subcovid13 <- subset(covid, `#CHR`==13)
mysnp <- subset(subcovid13, rsid=="rs12585036")

mergedcovidgtex13 = merge(subcovid13,atp11, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex13, rsid=="rs12585036")

Idx = which(mergedcovidgtex13$alt != mergedipfgtex13$ALT)
print(Idx)
mergedcovidgtex13[Idx,]$slope = - mergedcovidgtex13[Idx,]$slope

start = 112881427 - 500000
end = 112881427 + 500000

tmp1 = subset(mergedcovidgtex13, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex13, POS < end)
tmp2

mergedcovidgtex13[1,]
tmp3 = subset(mergedcovidgtex13, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex13$all_inv_var_meta_p), y=-log10(mergedcovidgtex13$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12585036')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

plot(x = mergedcovidgtex13$POS, y=-log10(mergedcovidgtex13$pval_nominal))
mysnp <- subset(mergedcovidgtex13, rsid == 'rs12585036')
points(x=(mysnp$POS), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.98e-13  5.85e-01  9.57e-15  6.62e-03  4.08e-01 
#[1] "PP abf for shared variant: 40.8%"


#wb
wb = fread("gtex_wholeblood_eqtl_allpairs_chr13_HG38.txt")
atp11 <- subset(wb, gene_id == 'ENSG00000068650.18')
atp11 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> atp11
atp11_exact <- subset(atp11, variant_id == 'chr13_112881427_C_T_b38')
atp11$chr = gsub("chr", "", atp11$chr) %>% as.integer()
atp11$bp <- as.numeric(atp11$bp)

subipf13 <- subset(ipf, chromosome==13)
mysnp <- subset(subipf13, rsid=="rs12585036")

mergedipfgtex13 = merge(subipf13,atp11, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex13, rsid=="rs12585036")

Idx = which(mergedipfgtex13$effect_allele != mergedipfgtex13$alt)
print(Idx)
mergedipfgtex13[Idx,]$slope = - mergedipfgtex13[Idx,]$slope

start = 112881427 - 500000
end = 112881427 + 500000

names(mergedipfgtex13)[names(mergedipfgtex13) == 'start'] <- 'position'
names(mergedipfgtex13)[names(mergedipfgtex13) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex13, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex13, position < end)
tmp2

mergedipfgtex13[1,]
tmp3 = subset(mergedipfgtex13, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex13$p), y=-log10(mergedipfgtex13$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12585036')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#4.54e-09  2.40e-01  2.65e-09  1.39e-01  6.21e-01 
#[1] "PP abf for shared variant: 62.1%"

#ATP11A - pvalue of gtex WB vs pvalue of IPF GWAS 
mergedipfgtex13$pval_nominal_transformed <- as.numeric(-log10(mergedipfgtex13$pval_nominal))
mergedipfgtex13$p_transformed <- as.numeric(-log10(mergedipfgtex13$p))
mergedipfgtex13$pc <- predict(prcomp(~pval_nominal_transformed+p_transformed, mergedipfgtex13))[,1]
mysnp <- subset(mergedipfgtex13, rsid == 'rs12585036')
p = ggplot(mergedipfgtex13, aes(p_transformed, pval_nominal_transformed, color = pc)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  theme_bw() +
  scale_color_gradient(low = "blue", high = "orangered")+
  xlab("-log10(P-value, IPF GWAS)")+
  ylab("-log10(P-value, GTEx Whole Blood)")
p
p+ geom_text_repel(
  data = mysnp,
  aes(label = rsid),
  nudge_x = 0,
  nudge_y = 1,
  size = 5,
  min.segment.length = 0,
  color = "black")

ggsave("rs12585036_IPF_GTEx_WB_pvalue_comparison.pdf", width = 6, height = 6)


subcovid13 <- subset(covid, `#CHR`==13)
mysnp <- subset(subcovid13, rsid=="rs12585036")

mergedcovidgtex13 = merge(subcovid13,atp11, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex13, rsid=="rs12585036")

Idx = which(mergedcovidgtex13$alt != mergedipfgtex13$ALT)
print(Idx)
mergedcovidgtex13[Idx,]$slope = - mergedcovidgtex13[Idx,]$slope

start = 112881427 - 500000
end = 112881427 + 500000

tmp1 = subset(mergedcovidgtex13, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex13, POS < end)
tmp2

mergedcovidgtex13[1,]
tmp3 = subset(mergedcovidgtex13, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex13$all_inv_var_meta_p), y=-log10(mergedcovidgtex13$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12585036')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.25e-13  2.41e-01  1.41e-14  9.70e-03  7.49e-01 
#[1] "PP abf for shared variant: 74.9%"

#FIBROBLASTS
fb = fread("gtex_fibroblasts_eqtl_allpairs_chr13_HG38.txt")
atp11 <- subset(fb, gene_id == 'ENSG00000068650.18')
atp11 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> atp11
atp11_exact <- subset(atp11, variant_id == 'chr13_112881427_C_T_b38')
atp11$chr = gsub("chr", "", atp11$chr) %>% as.integer()
atp11$bp <- as.numeric(atp11$bp)

subipf13 <- subset(ipf, chromosome==13)
mysnp <- subset(subipf13, rsid=="rs12585036")

mergedipfgtex13 = merge(subipf13,atp11, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex13, rsid=="rs12585036")

Idx = which(mergedipfgtex13$effect_allele != mergedipfgtex13$alt)
print(Idx)
mergedipfgtex13[Idx,]$slope = - mergedipfgtex13[Idx,]$slope

start = 112881427 - 500000
end = 112881427 + 500000

names(mergedipfgtex13)[names(mergedipfgtex13) == 'start'] <- 'position'
names(mergedipfgtex13)[names(mergedipfgtex13) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex13, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex13, position < end)
tmp2

mergedipfgtex13[1,]
tmp3 = subset(mergedipfgtex13, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex13$p), y=-log10(mergedipfgtex13$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12585036')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

plot(x=(mergedipfgtex13$position), y=-log10(mergedipfgtex13$pval_nominal))
mysnp <- subset(mergedipfgtex13, rsid == 'rs12585036')
points(x=-log10(mysnp$position), y=-log10(mysnp$pval_nominal),col='red')


tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.69e-10  8.93e-03  1.88e-08  9.91e-01  5.46e-04 
#[1] "PP abf for shared variant: 0.0546%"

subcovid13 <- subset(covid, `#CHR`==13)
mysnp <- subset(subcovid13, rsid=="rs12585036")

mergedcovidgtex13 = merge(subcovid13,atp11, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex13, rsid=="rs12585036")

Idx = which(mergedcovidgtex13$alt != mergedipfgtex13$ALT)
print(Idx)
mergedcovidgtex13[Idx,]$slope = - mergedcovidgtex13[Idx,]$slope

start = 112881427 - 500000
end = 112881427 + 500000

tmp1 = subset(mergedcovidgtex13, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex13, POS < end)
tmp2

mergedcovidgtex13[1,]
tmp3 = subset(mergedcovidgtex13, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex13$all_inv_var_meta_p), y=-log10(mergedcovidgtex13$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12585036')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.11e-12  8.20e-01  1.78e-13  1.31e-01  4.89e-02 
#[1] "PP abf for shared variant: 4.89%"

#LCL
lcl = fread("gtex_lcl_eqtl_allpairs_chr13_HG38.txt")
atp11 <- subset(lcl, gene_id == 'ENSG00000068650.18')
atp11 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> atp11
atp11_exact <- subset(atp11, variant_id == 'chr13_112881427_C_T_b38')
atp11$chr = gsub("chr", "", atp11$chr) %>% as.integer()
atp11$bp <- as.numeric(atp11$bp)

subipf13 <- subset(ipf, chromosome==13)
mysnp <- subset(subipf13, rsid=="rs12585036")

mergedipfgtex13 = merge(subipf13,atp11, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex13, rsid=="rs12585036")

Idx = which(mergedipfgtex13$effect_allele != mergedipfgtex13$alt)
print(Idx)
mergedipfgtex13[Idx,]$slope = - mergedipfgtex13[Idx,]$slope

start = 112881427 - 500000
end = 112881427 + 500000

names(mergedipfgtex13)[names(mergedipfgtex13) == 'start'] <- 'position'
names(mergedipfgtex13)[names(mergedipfgtex13) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex13, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex13, position < end)
tmp2

mergedipfgtex13[1,]
tmp3 = subset(mergedipfgtex13, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex13$p), y=-log10(mergedipfgtex13$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12585036')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.36e-08  7.16e-01  4.26e-09  2.25e-01  5.88e-02 
#[1] "PP abf for shared variant: 5.88%"

subcovid13 <- subset(covid, `#CHR`==13)
mysnp <- subset(subcovid13, rsid=="rs12585036")

mergedcovidgtex13 = merge(subcovid13,atp11, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex13, rsid=="rs12585036")

Idx = which(mergedcovidgtex13$alt != mergedipfgtex13$ALT)
print(Idx)
mergedcovidgtex13[Idx,]$slope = - mergedcovidgtex13[Idx,]$slope

start = 112881427 - 500000
end = 112881427 + 500000

tmp1 = subset(mergedcovidgtex13, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex13, POS < end)
tmp2

mergedcovidgtex13[1,]
tmp3 = subset(mergedcovidgtex13, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex13$all_inv_var_meta_p), y=-log10(mergedcovidgtex13$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12585036')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.21e-12  9.05e-01  2.52e-14  1.88e-02  7.65e-02 
#[1] "PP abf for shared variant: 7.65%"

#DPP9, rs12610495
#lung
lung = fread("gtex_eqtl_lung_allpairs_chr19_HG38.txt")
dpp9 <- subset(lung, gene_id == 'ENSG00000142002.16')
dpp9 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> dpp9
dpp9_exact <- subset(dpp9, variant_id == 'chr19_4717660_A_G_b38')
dpp9$chr = gsub("chr", "", dpp9$chr) %>% as.integer()
dpp9$bp <- as.numeric(dpp9$bp)

subipf19 <- subset(ipf, chromosome==19)
mysnp <- subset(subipf19, rsid=="rs12610495")

mergedipfgtex19 = merge(subipf19,dpp9, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex19, rsid=="rs12610495")

Idx = which(mergedipfgtex19$effect_allele != mergedipfgtex19$alt)
print(Idx)
mergedipfgtex19[Idx,]$slope = - mergedipfgtex19[Idx,]$slope

start = 4717660 - 500000
end = 4717660 + 500000

names(mergedipfgtex19)[names(mergedipfgtex19) == 'start'] <- 'position'
names(mergedipfgtex19)[names(mergedipfgtex19) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex19, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex19, position < end)
tmp2

mergedipfgtex19[1,]
tmp3 = subset(mergedipfgtex19, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex19$p), y=-log10(mergedipfgtex19$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12610495')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#5.25e-22  1.51e-12  3.48e-10  1.00e+00  1.06e-07 
#[1] "PP abf for shared variant: 1.06e-05%"

subcovid19 <- subset(covid, `#CHR`==19)
mysnp <- subset(subcovid19, rsid=="rs12610495")

mergedcovidgtex19 = merge(subcovid19,dpp9, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex19, rsid=="rs12610495")

Idx = which(mergedcovidgtex19$alt != mergedipfgtex19$ALT)
print(Idx)
mergedcovidgtex19[Idx,]$slope = - mergedcovidgtex19[Idx,]$slope

start = 4717660 - 500000
end = 4717660 + 500000

tmp1 = subset(mergedcovidgtex19, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex19, POS < end)
tmp2

mergedcovidgtex19[1,]
tmp3 = subset(mergedcovidgtex19, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex19$all_inv_var_meta_p), y=-log10(mergedcovidgtex19$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12610495')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.71e-64  2.57e-12  2.99e-52  1.00e+00  1.53e-07 
#[1] "PP abf for shared variant: 1.53e-05%"

#Whole Blood
wb = fread("gtex_wholeblood_eqtl_allpairs_chr19_HG38.txt")
dpp9 <- subset(wb, gene_id == 'ENSG00000142002.16')
dpp9 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> dpp9
dpp9_exact <- subset(dpp9, variant_id == 'chr19_4717660_A_G_b38')
dpp9$chr = gsub("chr", "", dpp9$chr) %>% as.integer()
dpp9$bp <- as.numeric(dpp9$bp)

subipf19 <- subset(ipf, chromosome==19)
mysnp <- subset(subipf19, rsid=="rs12610495")

mergedipfgtex19 = merge(subipf19,dpp9, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex19, rsid=="rs12610495")

Idx = which(mergedipfgtex19$effect_allele != mergedipfgtex19$alt)
print(Idx)
mergedipfgtex19[Idx,]$slope = - mergedipfgtex19[Idx,]$slope

start = 4717660 - 500000
end = 4717660 + 500000

names(mergedipfgtex19)[names(mergedipfgtex19) == 'start'] <- 'position'
names(mergedipfgtex19)[names(mergedipfgtex19) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex19, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex19, position < end)
tmp2

mergedipfgtex19[1,]
tmp3 = subset(mergedipfgtex19, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex19$p), y=-log10(mergedipfgtex19$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12610495')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.59e-27  2.18e-17  3.48e-10  1.00e+00  1.34e-11 
#[1] "PP abf for shared variant: 1.34e-09%"

subcovid19 <- subset(covid, `#CHR`==19)
mysnp <- subset(subcovid19, rsid=="rs12610495")

mergedcovidgtex19 = merge(subcovid19,dpp9, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex19, rsid=="rs12610495")

Idx = which(mergedcovidgtex19$alt != mergedipfgtex19$ALT)
print(Idx)
mergedcovidgtex19[Idx,]$slope = - mergedcovidgtex19[Idx,]$slope

start = 4717660 - 500000
end = 4717660 + 500000

tmp1 = subset(mergedcovidgtex19, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex19, POS < end)
tmp2

mergedcovidgtex19[1,]
tmp3 = subset(mergedcovidgtex19, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex19$all_inv_var_meta_p), y=-log10(mergedcovidgtex19$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12610495')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.16e-68  3.88e-17  2.99e-52  1.00e+00  1.26e-16 
#[1] "PP abf for shared variant: 1.26e-14%"

#Fibroblasts
fb = fread("gtex_fibroblasts_eqtl_allpairs_chr19_HG38.txt")
dpp9 <- subset(fb, gene_id == 'ENSG00000142002.16')
dpp9 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> dpp9
dpp9_exact <- subset(dpp9, variant_id == 'chr19_4717660_A_G_b38')
dpp9$chr = gsub("chr", "", dpp9$chr) %>% as.integer()
dpp9$bp <- as.numeric(dpp9$bp)

subipf19 <- subset(ipf, chromosome==19)
mysnp <- subset(subipf19, rsid=="rs12610495")

mergedipfgtex19 = merge(subipf19,dpp9, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex19, rsid=="rs12610495")

Idx = which(mergedipfgtex19$effect_allele != mergedipfgtex19$alt)
print(Idx)
mergedipfgtex19[Idx,]$slope = - mergedipfgtex19[Idx,]$slope

start = 4717660 - 500000
end = 4717660 + 500000

names(mergedipfgtex19)[names(mergedipfgtex19) == 'start'] <- 'position'
names(mergedipfgtex19)[names(mergedipfgtex19) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex19, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex19, position < end)
tmp2

mergedipfgtex19[1,]
tmp3 = subset(mergedipfgtex19, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex19$p), y=-log10(mergedipfgtex19$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12610495')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#8.19e-12  2.36e-02  3.21e-12  8.27e-03  9.68e-01 
#[1] "PP abf for shared variant: 96.8%"

#DPP9 - pvalue of gtex vs pvalue of IPF GWAS 
mergedipfgtex19$pval_nominal_transformed <- as.numeric(-log10(mergedipfgtex19$pval_nominal))
mergedipfgtex19$p_transformed <- as.numeric(-log10(mergedipfgtex19$p))
mergedipfgtex19$pc <- predict(prcomp(~pval_nominal_transformed+p_transformed, mergedipfgtex19))[,1]
mysnp <- subset(mergedipfgtex19, rsid == 'rs12610495')
p = ggplot(mergedipfgtex19, aes(p_transformed, pval_nominal_transformed, color = pc)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  theme_bw() +
  scale_color_gradient(low = "orangered", high = "blue")+
  xlab("-log10(P-value, IPF GWAS)")+
  ylab("-log10(P-value, GTEx Fibroblast)")
p
p+ geom_text_repel(
  data = mysnp,
  aes(label = rsid),
  nudge_x = 0,
  nudge_y = -2,
  size = 5,
  min.segment.length = 0,
  color = "black")

ggsave("rs12610495_IPF_GTEx_pvalue_comparison.pdf", width = 6, height = 6)



subcovid19 <- subset(covid, `#CHR`==19)
mysnp <- subset(subcovid19, rsid=="rs12610495")

mergedcovidgtex19 = merge(subcovid19,dpp9, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex19, rsid=="rs12610495")

Idx = which(mergedcovidgtex19$alt != mergedipfgtex19$ALT)
print(Idx)
mergedcovidgtex19[Idx,]$slope = - mergedcovidgtex19[Idx,]$slope

start = 4717660 - 500000
end = 4717660 + 500000

tmp1 = subset(mergedcovidgtex19, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex19, POS < end)
tmp2

mergedcovidgtex19[1,]
tmp3 = subset(mergedcovidgtex19, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex19$all_inv_var_meta_p), y=-log10(mergedcovidgtex19$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12610495')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.15e-54  2.39e-02  2.64e-54  7.85e-03  9.68e-01 
#[1] "PP abf for shared variant: 96.8%"

#LCL
lcl = fread("gtex_lcl_eqtl_allpairs_chr19_HG38.txt")
dpp9 <- subset(lcl, gene_id == 'ENSG00000142002.16')
dpp9 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> dpp9
dpp9_exact <- subset(dpp9, variant_id == 'chr19_4717660_A_G_b38')
dpp9$chr = gsub("chr", "", dpp9$chr) %>% as.integer()
dpp9$bp <- as.numeric(dpp9$bp)

subipf19 <- subset(ipf, chromosome==19)
mysnp <- subset(subipf19, rsid=="rs12610495")

mergedipfgtex19 = merge(subipf19,dpp9, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex19, rsid=="rs12610495")

Idx = which(mergedipfgtex19$effect_allele != mergedipfgtex19$alt)
print(Idx)
mergedipfgtex19[Idx,]$slope = - mergedipfgtex19[Idx,]$slope

start = 4717660 - 500000
end = 4717660 + 500000

names(mergedipfgtex19)[names(mergedipfgtex19) == 'start'] <- 'position'
names(mergedipfgtex19)[names(mergedipfgtex19) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex19, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex19, position < end)
tmp2

mergedipfgtex19[1,]
tmp3 = subset(mergedipfgtex19, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex19$p), y=-log10(mergedipfgtex19$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12610495')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.05e-10  5.89e-01  8.06e-11  2.32e-01  1.80e-01 
#[1] "PP abf for shared variant: 18%"

subcovid19 <- subset(covid, `#CHR`==19)
mysnp <- subset(subcovid19, rsid=="rs12610495")

mergedcovidgtex19 = merge(subcovid19,dpp9, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex19, rsid=="rs12610495")

Idx = which(mergedcovidgtex19$alt != mergedipfgtex19$ALT)
print(Idx)
mergedcovidgtex19[Idx,]$slope = - mergedcovidgtex19[Idx,]$slope

start = 4717660 - 500000
end = 4717660 + 500000

tmp1 = subset(mergedcovidgtex19, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex19, POS < end)
tmp2

mergedcovidgtex19[1,]
tmp3 = subset(mergedcovidgtex19, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex19$all_inv_var_meta_p), y=-log10(mergedcovidgtex19$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs12610495')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]
tmp5 <- subset(tmp4, maf > 0)

myres <- coloc.abf(dataset1=list(pvalues = tmp5$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp5$rsid),
                   dataset2=list(pvalues = tmp5$pval_nominal, N= 147,
                                 type="quant", snp=tmp5$rsid),
                   MAF=tmp5$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.25e-52  7.54e-01  6.33e-54  2.09e-02  2.25e-01 
#[1] "PP abf for shared variant: 22.5%"


#rs1105569
#Lung
# lung = fread("gtex_eqtl_lung_allpairs_chr17_HG38.txt")
# linc <- subset(lung, gene_id == 'ENSG00000204650.14')
# linc %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> linc
# linc_exact <- subset(linc, variant_id == 'chr17_45716022_C_T_b38')
# linc$chr = gsub("chr", "", linc$chr) %>% as.integer()
# linc$bp <- as.numeric(linc$bp)
# 
# subipf17 <- subset(ipf, chromosome==17)
# mysnp <- subset(subipf17, rsid=="rs1105569")
# 
# mergedipfgtex17 = merge(subipf17,linc, by.x = c("chromosome","start"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
# print(Idx)
# mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'
# 
# tmp1 = subset(mergedipfgtex17, position > start)
# tmp1
# plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedipfgtex17, position < end)
# tmp2
# 
# mergedipfgtex17[1,]
# tmp3 = subset(mergedipfgtex17, position > start & position < end)
# tmp3
# plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #4.63e-146 2.41e-130  1.25e-16  6.49e-01  3.51e-01 
# #[1] "PP abf for shared variant: 35.1%"
# 
# subcovid17 <- subset(covid, `#CHR`==17)
# mysnp <- subset(subcovid17, rsid=="rs1105569")
# 
# mergedcovidgtex17 = merge(subcovid17,linc, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
# print(Idx)
# mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# tmp1 = subset(mergedcovidgtex17, POS > start)
# tmp1
# plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedcovidgtex17, POS < end)
# tmp2
# 
# mergedcovidgtex17[1,]
# tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
# tmp3
# plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #3.16e-137 1.02e-127  9.72e-11  3.12e-01  6.88e-01 
# #[1] "PP abf for shared variant: 68.8%"
# 
# #Whole Blood
# wb = fread("gtex_wholeblood_eqtl_allpairs_chr17_HG38.txt")
# linc <- subset(wb, gene_id == 'ENSG00000204650.14')
# linc %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> linc
# linc_exact <- subset(linc, variant_id == 'chr17_45716022_C_T_b38')
# linc$chr = gsub("chr", "", linc$chr) %>% as.integer()
# linc$bp <- as.numeric(linc$bp)
# 
# subipf17 <- subset(ipf, chromosome==17)
# mysnp <- subset(subipf17, rsid=="rs1105569")
# 
# mergedipfgtex17 = merge(subipf17,linc, by.x = c("chromosome","start"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
# print(Idx)
# mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'
# 
# tmp1 = subset(mergedipfgtex17, position > start)
# tmp1
# plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedipfgtex17, position < end)
# tmp2
# 
# mergedipfgtex17[1,]
# tmp3 = subset(mergedipfgtex17, position > start & position < end)
# tmp3
# plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #7.45e-100  3.88e-84  1.30e-16  6.75e-01  3.25e-01 
# #[1] "PP abf for shared variant: 32.5%"
# 
# subcovid17 <- subset(covid, `#CHR`==17)
# mysnp <- subset(subcovid17, rsid=="rs1105569")
# 
# mergedcovidgtex17 = merge(subcovid17,linc, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
# print(Idx)
# mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# tmp1 = subset(mergedcovidgtex17, POS > start)
# tmp1
# plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedcovidgtex17, POS < end)
# tmp2
# 
# mergedcovidgtex17[1,]
# tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
# tmp3
# plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #4.35e-92  1.39e-82  9.77e-11  3.12e-01  6.88e-01 
# #[1] "PP abf for shared variant: 68.8%"
# 
# #Fibroblast
# fb = fread("gtex_fibroblasts_eqtl_allpairs_chr17_HG38.txt")
# linc <- subset(fb, gene_id == 'ENSG00000204650.14')
# linc %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> linc
# linc_exact <- subset(linc, variant_id == 'chr17_45716022_C_T_b38')
# linc$chr = gsub("chr", "", linc$chr) %>% as.integer()
# linc$bp <- as.numeric(linc$bp)
# 
# subipf17 <- subset(ipf, chromosome==17)
# mysnp <- subset(subipf17, rsid=="rs1105569")
# 
# mergedipfgtex17 = merge(subipf17,linc, by.x = c("chromosome","start"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
# print(Idx)
# mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'
# 
# tmp1 = subset(mergedipfgtex17, position > start)
# tmp1
# plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedipfgtex17, position < end)
# tmp2
# 
# mergedipfgtex17[1,]
# tmp3 = subset(mergedipfgtex17, position > start & position < end)
# tmp3
# plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #1.77e-148 9.23e-133  1.30e-16  6.78e-01  3.22e-01 
# #[1] "PP abf for shared variant: 32.2%"
# 
# subcovid17 <- subset(covid, `#CHR`==17)
# mysnp <- subset(subcovid17, rsid=="rs1105569")
# 
# mergedcovidgtex17 = merge(subcovid17,linc, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
# print(Idx)
# mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# tmp1 = subset(mergedcovidgtex17, POS > start)
# tmp1
# plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedcovidgtex17, POS < end)
# tmp2
# 
# mergedcovidgtex17[1,]
# tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
# tmp3
# plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #9.89e-142 3.09e-132  9.86e-11  3.08e-01  6.92e-01 
# #[1] "PP abf for shared variant: 69.2%"
# 
# #LCL
# lcl = fread("gtex_lcl_eqtl_allpairs_chr17_HG38.txt")
# linc <- subset(lcl, gene_id == 'ENSG00000204650.14')
# linc %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> linc
# linc_exact <- subset(linc, variant_id == 'chr17_45716022_C_T_b38')
# linc$chr = gsub("chr", "", linc$chr) %>% as.integer()
# linc$bp <- as.numeric(linc$bp)
# 
# subipf17 <- subset(ipf, chromosome==17)
# mysnp <- subset(subipf17, rsid=="rs1105569")
# 
# mergedipfgtex17 = merge(subipf17,linc, by.x = c("chromosome","start"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
# print(Idx)
# mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'
# 
# tmp1 = subset(mergedipfgtex17, position > start)
# tmp1
# plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedipfgtex17, position < end)
# tmp2
# 
# mergedipfgtex17[1,]
# tmp3 = subset(mergedipfgtex17, position > start & position < end)
# tmp3
# plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #2.97e-28  1.54e-12  1.30e-16  6.78e-01  3.22e-01 
# #[1] "PP abf for shared variant: 32.2%"
# 
# subcovid17 <- subset(covid, `#CHR`==17)
# mysnp <- subset(subcovid17, rsid=="rs1105569")
# 
# mergedcovidgtex17 = merge(subcovid17,linc, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
# print(Idx)
# mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# tmp1 = subset(mergedcovidgtex17, POS > start)
# tmp1
# plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedcovidgtex17, POS < end)
# tmp2
# 
# mergedcovidgtex17[1,]
# tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
# tmp3
# plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #1.23e-21  3.88e-12  1.01e-10  3.18e-01  6.82e-01 
# #[1] "PP abf for shared variant: 68.2%"

#Lung
lung = fread("gtex_eqtl_lung_allpairs_chr17_HG38.txt")
kansl <- subset(lung, gene_id == 'ENSG00000120071.13')
kansl %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> kansl
kansl_exact <- subset(kansl, variant_id == 'chr17_45716022_C_T_b38')
kansl$chr = gsub("chr", "", kansl$chr) %>% as.integer()
kansl$bp <- as.numeric(kansl$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,kansl, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.27e-36  1.70e-20  1.92e-16  1.00e+00  2.78e-16 
#[1] "PP abf for shared variant: 2.78e-14%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,kansl, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.79e-13  2.50e-03  9.62e-11  3.09e-01  6.89e-01 
#[1] "PP abf for shared variant: 68.9%"


#Whole Blood
wb = fread("gtex_wholeblood_eqtl_allpairs_chr17_HG38.txt")
kansl <- subset(wb, gene_id == 'ENSG00000120071.13')
kansl %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> kansl
kansl_exact <- subset(kansl, variant_id == 'chr17_45716022_C_T_b38')
kansl$chr = gsub("chr", "", kansl$chr) %>% as.integer()
kansl$bp <- as.numeric(kansl$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,kansl, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.05e-40  5.47e-25  1.92e-16  1.00e+00  2.76e-16 
#[1] "PP abf for shared variant: 2.76e-14%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,kansl, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#4.97e-17  1.59e-07  9.72e-11  3.10e-01  6.90e-01 
#[1] "PP abf for shared variant: 69%"


#Fibroblast
fb = fread("gtex_fibroblasts_eqtl_allpairs_chr17_HG38.txt")
kansl <- subset(fb, gene_id == 'ENSG00000120071.13')
kansl %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> kansl
kansl_exact <- subset(kansl, variant_id == 'chr17_45716022_C_T_b38')
kansl$chr = gsub("chr", "", kansl$chr) %>% as.integer()
kansl$bp <- as.numeric(kansl$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,kansl, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#5.06e-27  2.63e-11  1.33e-16  6.90e-01  3.10e-01 
#[1] "PP abf for shared variant: 31%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,kansl, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.54e-20  7.95e-11  9.74e-11  3.04e-01  6.96e-01 
#[1] "PP abf for shared variant: 69.6%"

plot(x=-(mergedcovidgtex17$all_inv_var_meta_beta), y=(mergedcovidgtex17$slope))
mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-(mysnp$all_inv_var_meta_beta), y=(mysnp$slope), col='red')
plot(x=-(mergedipfgtex17$beta), y=(mergedipfgtex17$slope))



#LCL
lcl = fread("gtex_lcl_eqtl_allpairs_chr17_HG38.txt")
kansl <- subset(lcl, gene_id == 'ENSG00000120071.13')
kansl %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> kansl
kansl_exact <- subset(kansl, variant_id == 'chr17_45716022_C_T_b38')
kansl$chr = gsub("chr", "", kansl$chr) %>% as.integer()
kansl$bp <- as.numeric(kansl$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,kansl, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.05e-20  1.59e-04  1.29e-16  6.69e-01  3.31e-01 
#[1] "PP abf for shared variant: 33.1%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,kansl, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.27e-13  4.01e-04  9.84e-11  3.09e-01  6.91e-01 
#[1] "PP abf for shared variant: 69.1%"

#CRHR1
#rs1105569 is not an eqtl for CRHR1 in the lung, whole blood, fb, or lcl

#MAPT
lung = fread("gtex_eqtl_lung_allpairs_chr17_HG38.txt")
mapt <- subset(lung, gene_id == 'ENSG00000186868.15')
mapt %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> mapt
mapt_exact <- subset(mapt, variant_id == 'chr17_45716022_C_T_b38')
mapt$chr = gsub("chr", "", mapt$chr) %>% as.integer()
mapt$bp <- as.numeric(mapt$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,mapt, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.94e-40  2.05e-24  1.28e-16  6.68e-01  3.32e-01 
#[1] "PP abf for shared variant: 33.2%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,mapt, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.16e-33  1.02e-23  9.65e-11  3.10e-01  6.90e-01 
#[1] "PP abf for shared variant: 69%"

#whole blood
wb = fread("gtex_wholeblood_eqtl_allpairs_chr17_HG38.txt")
mapt <- subset(wb, gene_id == 'ENSG00000186868.15')
mapt %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> mapt
mapt_exact <- subset(mapt, variant_id == 'chr17_45716022_C_T_b38')
mapt$chr = gsub("chr", "", mapt$chr) %>% as.integer()
mapt$bp <- as.numeric(mapt$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,mapt, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.23e-16  6.41e-01  5.01e-17  2.61e-01  9.87e-02 
#[1] "PP abf for shared variant: 9.87%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,mapt, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.52e-10  8.06e-01  1.99e-11  6.35e-02  1.30e-01 
#[1] "PP abf for shared variant: 13%"

#Fibroblast
fb = fread("gtex_fibroblasts_eqtl_allpairs_chr17_HG38.txt")
mapt <- subset(fb, gene_id == 'ENSG00000186868.15')
mapt %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> mapt
mapt_exact <- subset(mapt, variant_id == 'chr17_45716022_C_T_b38')
mapt$chr = gsub("chr", "", mapt$chr) %>% as.integer()
mapt$bp <- as.numeric(mapt$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,mapt, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#6.25e-17  3.25e-01  9.07e-17  4.72e-01  2.03e-01 
#[1] "PP abf for shared variant: 20.3%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,mapt, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.71e-10  5.35e-01  4.74e-11  1.48e-01  3.17e-01 
#[1] "PP abf for shared variant: 31.7%"


#LCL
lcl = fread("gtex_lcl_eqtl_allpairs_chr17_HG38.txt")
mapt <- subset(lcl, gene_id == 'ENSG00000186868.15')
mapt %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> mapt
mapt_exact <- subset(mapt, variant_id == 'chr17_45716022_C_T_b38')
mapt$chr = gsub("chr", "", mapt$chr) %>% as.integer()
mapt$bp <- as.numeric(mapt$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,mapt, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.37e-16  7.12e-01  4.22e-17  2.20e-01  6.84e-02 
#[1] "PP abf for shared variant: 6.84%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,mapt, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.77e-10  8.71e-01  1.36e-11  4.26e-02  8.62e-02 
#[1] "PP abf for shared variant: 8.62%"


# #LINC02210
# lung = fread("gtex_eqtl_lung_allpairs_chr17_HG38.txt")
# linc <- subset(lung, gene_id == 'ENSG00000204650.14')
# linc %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> linc
# linc_exact <- subset(linc, variant_id == 'chr17_45716022_C_T_b38')
# linc$chr = gsub("chr", "", linc$chr) %>% as.integer()
# linc$bp <- as.numeric(linc$bp)
# 
# subipf17 <- subset(ipf, chromosome==17)
# mysnp <- subset(subipf17, rsid=="rs1105569")
# 
# mergedipfgtex17 = merge(subipf17,linc, by.x = c("chromosome","start"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
# print(Idx)
# mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'
# 
# tmp1 = subset(mergedipfgtex17, position > start)
# tmp1
# plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedipfgtex17, position < end)
# tmp2
# 
# mergedipfgtex17[1,]
# tmp3 = subset(mergedipfgtex17, position > start & position < end)
# tmp3
# plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #4.63e-146 2.41e-130  1.25e-16  6.49e-01  3.51e-01 
# #[1] "PP abf for shared variant: 35.1%"
# 
# subcovid17 <- subset(covid, `#CHR`==17)
# mysnp <- subset(subcovid17, rsid=="rs1105569")
# 
# mergedcovidgtex17 = merge(subcovid17,linc, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
# print(Idx)
# mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# tmp1 = subset(mergedcovidgtex17, POS > start)
# tmp1
# plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedcovidgtex17, POS < end)
# tmp2
# 
# mergedcovidgtex17[1,]
# tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
# tmp3
# plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #3.16e-137 1.02e-127  9.72e-11  3.12e-01  6.88e-01 
# #[1] "PP abf for shared variant: 68.8%"

#whole blood
wb = fread("gtex_wholeblood_eqtl_allpairs_chr17_HG38.txt")
linc <- subset(wb, gene_id == 'ENSG00000204650.14')
linc %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> linc
linc_exact <- subset(linc, variant_id == 'chr17_45716022_C_T_b38')
linc$chr = gsub("chr", "", linc$chr) %>% as.integer()
linc$bp <- as.numeric(linc$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,linc, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.45e-100  3.88e-84  1.30e-16  6.75e-01  3.25e-01 
#[1] "PP abf for shared variant: 32.5%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,linc, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#4.35e-92  1.39e-82  9.77e-11  3.12e-01  6.88e-01 
#[1] "PP abf for shared variant: 68.8%"

#fibroblasts
fb = fread("gtex_fibroblasts_eqtl_allpairs_chr17_HG38.txt")
linc <- subset(fb, gene_id == 'ENSG00000204650.14')
linc %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> linc
linc_exact <- subset(linc, variant_id == 'chr17_45716022_C_T_b38')
linc$chr = gsub("chr", "", linc$chr) %>% as.integer()
linc$bp <- as.numeric(linc$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,linc, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.77e-148 9.23e-133  1.30e-16  6.78e-01  3.22e-01 
#[1] "PP abf for shared variant: 32.2%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,linc, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#9.89e-142 3.09e-132  9.86e-11  3.08e-01  6.92e-01 
#[1] "PP abf for shared variant: 69.2%"

# #LCL
# lcl = fread("gtex_lcl_eqtl_allpairs_chr17_HG38.txt")
# linc <- subset(lcl, gene_id == 'ENSG00000204650.14')
# linc %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> linc
# linc_exact <- subset(linc, variant_id == 'chr17_45716022_C_T_b38')
# linc$chr = gsub("chr", "", linc$chr) %>% as.integer()
# linc$bp <- as.numeric(linc$bp)
# 
# subipf17 <- subset(ipf, chromosome==17)
# mysnp <- subset(subipf17, rsid=="rs1105569")
# 
# mergedipfgtex17 = merge(subipf17,linc, by.x = c("chromosome","start"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
# print(Idx)
# mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'
# 
# tmp1 = subset(mergedipfgtex17, position > start)
# tmp1
# plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedipfgtex17, position < end)
# tmp2
# 
# mergedipfgtex17[1,]
# tmp3 = subset(mergedipfgtex17, position > start & position < end)
# tmp3
# plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #2.97e-28  1.54e-12  1.30e-16  6.78e-01  3.22e-01 
# #[1] "PP abf for shared variant: 32.2%"
# 
# subcovid17 <- subset(covid, `#CHR`==17)
# mysnp <- subset(subcovid17, rsid=="rs1105569")
# 
# mergedcovidgtex17 = merge(subcovid17,linc, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
# print(Idx)
# mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# tmp1 = subset(mergedcovidgtex17, POS > start)
# tmp1
# plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedcovidgtex17, POS < end)
# tmp2
# 
# mergedcovidgtex17[1,]
# tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
# tmp3
# plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #1.23e-21  3.88e-12  1.01e-10  3.18e-01  6.82e-01 
# #[1] "PP abf for shared variant: 68.2%"


#sppl2c
#rs1105569 is not an eqtl for sppl2c in the lung, whole blood, fb, or lcl

#ARHGAP27 in LUNG
ARHGAP27 <- subset(lung, gene_id == 'ENSG00000159314.11')
ARHGAP27 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARHGAP27
ARHGAP27_exact <- subset(ARHGAP27, variant_id == 'chr17_45716022_C_T_b38')
ARHGAP27$chr = gsub("chr", "", ARHGAP27$chr) %>% as.integer()
ARHGAP27$bp <- as.numeric(ARHGAP27$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARHGAP27, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.41e-16  7.33e-01  4.36e-17  2.27e-01  4.05e-02 
#[1] "PP abf for shared variant: 4.05%"


mergedcovidgtex17 = merge(subcovid17,ARHGAP27, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.86e-10  9.19e-01  8.58e-12  2.75e-02  5.36e-02 
#[1] "PP abf for shared variant: 5.36%"

#ARHGAP27 in WB
ARHGAP27 <- subset(wb, gene_id == 'ENSG00000159314.11')
ARHGAP27 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARHGAP27
ARHGAP27_exact <- subset(ARHGAP27, variant_id == 'chr17_45716022_C_T_b38')
ARHGAP27$chr = gsub("chr", "", ARHGAP27$chr) %>% as.integer()
ARHGAP27$bp <- as.numeric(ARHGAP27$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARHGAP27, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.52e-16  7.90e-01  3.30e-17  1.72e-01  3.89e-02 
#[1] "PP abf for shared variant: 3.89%"

mergedcovidgtex17 = merge(subcovid17,ARHGAP27, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.90e-10  9.29e-01  7.65e-12  2.44e-02  4.64e-02 
#[1] "PP abf for shared variant: 4.64%"

#ARHGAP27 in FB
ARHGAP27 <- subset(fb, gene_id == 'ENSG00000159314.11')
ARHGAP27 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARHGAP27
ARHGAP27_exact <- subset(ARHGAP27, variant_id == 'chr17_45716022_C_T_b38')
ARHGAP27$chr = gsub("chr", "", ARHGAP27$chr) %>% as.integer()
ARHGAP27$bp <- as.numeric(ARHGAP27$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARHGAP27, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.20e-25  6.26e-10  1.92e-16  1.00e+00  2.92e-10 
#[1] "PP abf for shared variant: 2.92e-08%"

mergedcovidgtex17 = merge(subcovid17,ARHGAP27, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#9.42e-17  2.95e-07  3.20e-10  1.00e+00  1.42e-07 
#[1] "PP abf for shared variant: 1.42e-05%"

#ARHGAP27 in LCL
ARHGAP27 <- subset(lcl, gene_id == 'ENSG00000159314.11')
ARHGAP27 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARHGAP27
ARHGAP27_exact <- subset(ARHGAP27, variant_id == 'chr17_45716022_C_T_b38')
ARHGAP27$chr = gsub("chr", "", ARHGAP27$chr) %>% as.integer()
ARHGAP27$bp <- as.numeric(ARHGAP27$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARHGAP27, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.48e-21  1.29e-05  1.31e-16  6.81e-01  3.19e-01 
#[1] "PP abf for shared variant: 31.9%"

mergedcovidgtex17 = merge(subcovid17,ARHGAP27, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#9.64e-15  3.03e-05  9.79e-11  3.07e-01  6.93e-01 
#[1] "PP abf for shared variant: 69.3%"


#LRRC37A2 in LUNG
LRRC37A2 <- subset(lung, gene_id == 'ENSG00000238083.7')
LRRC37A2 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> LRRC37A2
LRRC37A2_exact <- subset(LRRC37A2, variant_id == 'chr17_45716022_C_T_b38')
LRRC37A2$chr = gsub("chr", "", LRRC37A2$chr) %>% as.integer()
LRRC37A2$bp <- as.numeric(LRRC37A2$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,LRRC37A2, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.56e-80  8.11e-65  1.28e-16  6.68e-01  3.32e-01 
#[1] "PP abf for shared variant: 33.2%"

mergedcovidgtex17 = merge(subcovid17,LRRC37A2, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#4.88e-73  1.57e-63  9.27e-11  2.97e-01  7.03e-01 
#[1] "PP abf for shared variant: 70.3%"

#LRRC37A2 in WB
LRRC37A2 <- subset(wb, gene_id == 'ENSG00000238083.7')
LRRC37A2 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> LRRC37A2
LRRC37A2_exact <- subset(LRRC37A2, variant_id == 'chr17_45716022_C_T_b38')
LRRC37A2$chr = gsub("chr", "", LRRC37A2$chr) %>% as.integer()
LRRC37A2$bp <- as.numeric(LRRC37A2$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,LRRC37A2, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.19e-64  1.66e-48  1.31e-16  6.80e-01  3.20e-01 
#[1] "PP abf for shared variant: 32%"

mergedcovidgtex17 = merge(subcovid17,LRRC37A2, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.92e-57  1.25e-47  9.26e-11  2.96e-01  7.04e-01 
#[1] "PP abf for shared variant: 70.4%"

#LRRC37A2 in FB
LRRC37A2 <- subset(fb, gene_id == 'ENSG00000238083.7')
LRRC37A2 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> LRRC37A2
LRRC37A2_exact <- subset(LRRC37A2, variant_id == 'chr17_45716022_C_T_b38')
LRRC37A2$chr = gsub("chr", "", LRRC37A2$chr) %>% as.integer()
LRRC37A2$bp <- as.numeric(LRRC37A2$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,LRRC37A2, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.80e-51  9.36e-36  1.31e-16  6.83e-01  3.17e-01 
#[1] "PP abf for shared variant: 31.7%"

mergedcovidgtex17 = merge(subcovid17,LRRC37A2, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#5.34e-45  1.67e-35  9.74e-11  3.04e-01  6.96e-01 
#[1] "PP abf for shared variant: 69.6%"

#LRRC37A2 in lcl
LRRC37A2 <- subset(lcl, gene_id == 'ENSG00000238083.7')
LRRC37A2 %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> LRRC37A2
LRRC37A2_exact <- subset(LRRC37A2, variant_id == 'chr17_45716022_C_T_b38')
LRRC37A2$chr = gsub("chr", "", LRRC37A2$chr) %>% as.integer()
LRRC37A2$bp <- as.numeric(LRRC37A2$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,LRRC37A2, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#6.60e-22  3.43e-06  1.29e-16  6.72e-01  3.28e-01 
#[1] "PP abf for shared variant: 32.8%"

mergedcovidgtex17 = merge(subcovid17,LRRC37A2, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.75e-15  8.63e-06  9.77e-11  3.06e-01  6.94e-01 
#[1] "PP abf for shared variant: 69.4%"

# #FOR FAM215B
# FAM <- subset(lung, gene_id == 'ENSG00000232300.1')
# FAM %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> FAM
# FAM_exact <- subset(FAM, variant_id == 'chr17_45716022_C_T_b38')
# FAM$chr = gsub("chr", "", FAM$chr) %>% as.integer()
# FAM$bp <- as.numeric(FAM$bp)
# 
# subipf17 <- subset(ipf, chromosome==17)
# mysnp <- subset(subipf17, rsid=="rs1105569")
# 
# mergedipfgtex17 = merge(subipf17,FAM, by.x = c("chromosome","start"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
# print(Idx)
# mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
# names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'
# 
# tmp1 = subset(mergedipfgtex17, position > start)
# tmp1
# plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedipfgtex17, position < end)
# tmp2
# 
# mergedipfgtex17[1,]
# tmp3 = subset(mergedipfgtex17, position > start & position < end)
# tmp3
# plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# #9.26e-24  4.82e-08  1.28e-16  6.67e-01  3.33e-01 
# #[1] "PP abf for shared variant: 33.3%"
# 
# mergedcovidgtex17 = merge(subcovid17,FAM, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
# mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")
# 
# Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
# print(Idx)
# mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope
# 
# start = 45716022 - 500000
# end = 45716022 + 500000
# 
# tmp1 = subset(mergedcovidgtex17, POS > start)
# tmp1
# plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
# tmp2 = subset(mergedcovidgtex17, POS < end)
# tmp2
# 
# mergedcovidgtex17[1,]
# tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
# tmp3
# plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
# plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))
# 
# mysnp <- subset(tmp3, rsid == 'rs1105569')
# points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')
# 
# tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()
# 
# tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]
# 
# myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
#                    dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
#                                  type="quant", snp=tmp4$rsid),
#                    MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
# 
# #PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# # 3.59e-17  1.15e-07  9.75e-11  3.13e-01  6.87e-01 
# # [1] "PP abf for shared variant: 68.7%"
# 
# #FAM in WB - No eQTLs for FAM215b in Whole Blood, FB, LCL 

#LRRC37A in LUNG
LRRC37A <- subset(lung, gene_id == 'ENSG00000176681.14')
LRRC37A %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> LRRC37A
LRRC37A_exact <- subset(LRRC37A, variant_id == 'chr17_45716022_C_T_b38')
LRRC37A$chr = gsub("chr", "", LRRC37A$chr) %>% as.integer()
LRRC37A$bp <- as.numeric(LRRC37A$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,LRRC37A, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.40e-45  7.30e-30  1.92e-16  1.00e+00  4.53e-07 
#[1] "PP abf for shared variant: 4.53e-05%"

mergedcovidgtex17 = merge(subcovid17,LRRC37A, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#8.56e-23  2.75e-13  9.63e-11  3.09e-01  6.91e-01 
#[1] "PP abf for shared variant: 69.1%"

#LRRC37A in WB
LRRC37A <- subset(wb, gene_id == 'ENSG00000176681.14')
LRRC37A %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> LRRC37A
LRRC37A_exact <- subset(LRRC37A, variant_id == 'chr17_45716022_C_T_b38')
LRRC37A$chr = gsub("chr", "", LRRC37A$chr) %>% as.integer()
LRRC37A$bp <- as.numeric(LRRC37A$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,LRRC37A, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#6.39e-39  3.33e-23  1.92e-16  1.00e+00  4.79e-07 
#[1] "PP abf for shared variant: 4.79e-05%"

mergedcovidgtex17 = merge(subcovid17,LRRC37A, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#6.90e-20  2.21e-10  9.48e-11  3.03e-01  6.97e-01 
#[1] "PP abf for shared variant: 69.7%"

#LRRC37A in FB
LRRC37A <- subset(fb, gene_id == 'ENSG00000176681.14')
LRRC37A %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> LRRC37A
LRRC37A_exact <- subset(LRRC37A, variant_id == 'chr17_45716022_C_T_b38')
LRRC37A$chr = gsub("chr", "", LRRC37A$chr) %>% as.integer()
LRRC37A$bp <- as.numeric(LRRC37A$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,LRRC37A, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.46e-33  7.57e-18  1.92e-16  1.00e+00  5.10e-07 
#[1] "PP abf for shared variant: 5.1e-05%"

mergedcovidgtex17 = merge(subcovid17,LRRC37A, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#5.01e-11  1.57e-01  9.62e-11  3.00e-01  5.43e-01 
#[1] "PP abf for shared variant: 54.3%"

#LRRC37A in lcl
LRRC37A <- subset(lcl, gene_id == 'ENSG00000176681.14')
LRRC37A %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> LRRC37A
LRRC37A_exact <- subset(LRRC37A, variant_id == 'chr17_45716022_C_T_b38')
LRRC37A$chr = gsub("chr", "", LRRC37A$chr) %>% as.integer()
LRRC37A$bp <- as.numeric(LRRC37A$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,LRRC37A, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.35e-16  7.02e-01  4.79e-17  2.49e-01  4.83e-02 
#[1] "PP abf for shared variant: 4.83%"

mergedcovidgtex17 = merge(subcovid17,LRRC37A, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.88e-10  9.05e-01  1.03e-11  3.24e-02  6.26e-02 
#[1] "PP abf for shared variant: 6.26%"

#ARL17A in LUNG
lung = fread("gtex_eqtl_lung_allpairs_chr17_HG38.txt")
ARL17A <- subset(lung, gene_id == 'ENSG00000185829.17')
ARL17A %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARL17A
ARL17A_exact <- subset(ARL17A, variant_id == 'chr17_45716022_C_T_b38')
ARL17A$chr = gsub("chr", "", ARL17A$chr) %>% as.integer()
ARL17A$bp <- as.numeric(ARL17A$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARL17A, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.47e-34  3.89e-18  1.92e-16  1.00e+00  6.23e-09 
#[1] "PP abf for shared variant: 6.23e-07%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,ARL17A, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.07e-19  6.64e-10  9.60e-11  3.08e-01  6.92e-01 
#[1] "PP abf for shared variant: 69.2%"

#ARL17A in WB
wb = fread("gtex_wholeblood_eqtl_allpairs_chr17_HG38.txt")
ARL17A <- subset(wb, gene_id == 'ENSG00000185829.17')
ARL17A %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARL17A
ARL17A_exact <- subset(ARL17A, variant_id == 'chr17_45716022_C_T_b38')
ARL17A$chr = gsub("chr", "", ARL17A$chr) %>% as.integer()
ARL17A$bp <- as.numeric(ARL17A$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARL17A, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.40e-51  1.25e-35  1.92e-16  1.00e+00  2.73e-16 
#[1] "PP abf for shared variant: 2.73e-14%"

mergedcovidgtex17 = merge(subcovid17,ARL17A, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#8.90e-13  2.85e-03  1.00e-10  3.21e-01  6.76e-01 
#[1] "PP abf for shared variant: 67.6%"

#ARL17A in FB
fb <- fread("gtex_fibroblasts_eqtl_allpairs_chr17_HG38.txt")
ARL17A <- subset(fb, gene_id == 'ENSG00000185829.17')
ARL17A %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARL17A
ARL17A_exact <- subset(ARL17A, variant_id == 'chr17_45716022_C_T_b38')
ARL17A$chr = gsub("chr", "", ARL17A$chr) %>% as.integer()
ARL17A$bp <- as.numeric(ARL17A$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARL17A, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#7.59e-41  3.95e-25  1.92e-16  1.00e+00  3.69e-12 
#[1] "PP abf for shared variant: 3.69e-10%"

mergedcovidgtex17 = merge(subcovid17,ARL17A, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.77e-23  1.18e-13  1.03e-10  3.21e-01  6.79e-01 
#[1] "PP abf for shared variant: 67.9%"

#ARL17A in lcl
lcl <- fread("gtex_lcl_eqtl_allpairs_chr17_HG38.txt")
ARL17A <- subset(lcl, gene_id == 'ENSG00000185829.17')
ARL17A %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARL17A
ARL17A_exact <- subset(ARL17A, variant_id == 'chr17_45716022_C_T_b38')
ARL17A$chr = gsub("chr", "", ARL17A$chr) %>% as.integer()
ARL17A$bp <- as.numeric(ARL17A$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARL17A, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#6.05e-18  3.15e-02  1.43e-16  7.45e-01  2.24e-01 
#[1] "PP abf for shared variant: 22.4%"

mergedcovidgtex17 = merge(subcovid17,ARL17A, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#3.10e-11  9.74e-02  8.99e-11  2.82e-01  6.20e-01 
#[1] "PP abf for shared variant: 62%"


#ARL17B
#ARL17B in LUNG
lung = fread("gtex_eqtl_lung_allpairs_chr17_HG38.txt")
ARL17B <- subset(lung, gene_id == 'ENSG00000228696.8')
ARL17B %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARL17B
ARL17B_exact <- subset(ARL17B, variant_id == 'chr17_45716022_C_T_b38')
ARL17B$chr = gsub("chr", "", ARL17B$chr) %>% as.integer()
ARL17B$bp <- as.numeric(ARL17B$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARL17B, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.99e-22  1.03e-06  1.92e-16  1.00e+00  5.07e-08 
#[1] "PP abf for shared variant: 5.07e-06%"

subcovid17 <- subset(covid, `#CHR`==17)
mysnp <- subset(subcovid17, rsid=="rs1105569")

mergedcovidgtex17 = merge(subcovid17,ARL17B, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 515,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.87e-10  9.23e-01  9.39e-12  3.02e-02  4.66e-02 
#[1] "PP abf for shared variant: 4.66%"

#ARL17B in WB
wb = fread("gtex_wholeblood_eqtl_allpairs_chr17_HG38.txt")
ARL17B <- subset(wb, gene_id == 'ENSG00000228696.8')
ARL17B %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARL17B
ARL17B_exact <- subset(ARL17B, variant_id == 'chr17_45716022_C_T_b38')
ARL17B$chr = gsub("chr", "", ARL17B$chr) %>% as.integer()
ARL17B$bp <- as.numeric(ARL17B$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARL17B, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.83e-40  1.47e-24  1.92e-16  1.00e+00  2.73e-16 
#[1] "PP abf for shared variant: 2.73e-14%"

mergedcovidgtex17 = merge(subcovid17,ARL17B, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 670,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#5.36e-12  1.71e-02  9.73e-11  3.11e-01  6.72e-01 
#[1] "PP abf for shared variant: 67.2%"

#ARL17A in FB
fb <- fread("gtex_fibroblasts_eqtl_allpairs_chr17_HG38.txt")
ARL17B <- subset(fb, gene_id == 'ENSG00000228696.8')
ARL17B %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARL17B
ARL17B_exact <- subset(ARL17B, variant_id == 'chr17_45716022_C_T_b38')
ARL17B$chr = gsub("chr", "", ARL17B$chr) %>% as.integer()
ARL17B$bp <- as.numeric(ARL17B$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARL17B, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#4.07e-24  2.12e-08  1.92e-16  1.00e+00  2.79e-09 
#[1] "PP abf for shared variant: 2.79e-07%"

mergedcovidgtex17 = merge(subcovid17,ARL17B, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 483,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#2.88e-10  9.01e-01  1.46e-11  4.56e-02  5.34e-02 
#[1] "PP abf for shared variant: 5.34%"

#ARL17B in lcl
lcl <- fread("gtex_lcl_eqtl_allpairs_chr17_HG38.txt")
ARL17B <- subset(lcl, gene_id == 'ENSG00000228696.8')
ARL17B %>% separate(variant_id, c("chr","bp","ref","alt"), sep="_", remove = FALSE) -> ARL17B
ARL17B_exact <- subset(ARL17B, variant_id == 'chr17_45716022_C_T_b38')
ARL17B$chr = gsub("chr", "", ARL17B$chr) %>% as.integer()
ARL17B$bp <- as.numeric(ARL17B$bp)

subipf17 <- subset(ipf, chromosome==17)
mysnp <- subset(subipf17, rsid=="rs1105569")

mergedipfgtex17 = merge(subipf17,ARL17B, by.x = c("chromosome","start"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedipfgtex17$effect_allele != mergedipfgtex17$alt)
print(Idx)
mergedipfgtex17[Idx,]$slope = - mergedipfgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

names(mergedipfgtex17)[names(mergedipfgtex17) == 'start'] <- 'position'
names(mergedipfgtex17)[names(mergedipfgtex17) == 'end'] <- 'position2'

tmp1 = subset(mergedipfgtex17, position > start)
tmp1
plot(x=-log10(tmp1$p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedipfgtex17, position < end)
tmp2

mergedipfgtex17[1,]
tmp3 = subset(mergedipfgtex17, position > start & position < end)
tmp3
plot(x=-log10(tmp3$p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedipfgtex17$p), y=-log10(mergedipfgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

myres <- coloc.abf(dataset1=list(pvalues = tmp4$p,  N= 4125+20464,type="cc", s=4125/(4125+20464), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf.x, p1=1e-4,p2=1e-4, p12= 1e-5)

#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#6.55e-17  3.41e-01  9.31e-17  4.85e-01  1.74e-01 
#[1] "PP abf for shared variant: 17.4%"

mergedcovidgtex17 = merge(subcovid17,ARL17B, by.x = c("#CHR","POS"), by.y = c("chr","bp"))
mysnp <- subset(mergedipfgtex17, rsid=="rs1105569")

Idx = which(mergedcovidgtex17$alt != mergedipfgtex17$ALT)
print(Idx)
mergedcovidgtex17[Idx,]$slope = - mergedcovidgtex17[Idx,]$slope

start = 45716022 - 500000
end = 45716022 + 500000

tmp1 = subset(mergedcovidgtex17, POS > start)
tmp1
plot(x=-log10(tmp1$all_inv_var_meta_p), y=-log10(tmp1$pval_nominal))
tmp2 = subset(mergedcovidgtex17, POS < end)
tmp2

mergedcovidgtex17[1,]
tmp3 = subset(mergedcovidgtex17, POS > start & POS < end)
tmp3
plot(x=-log10(tmp3$all_inv_var_meta_p), y=-log10(tmp3$pval_nominal))
plot(x=-log10(mergedcovidgtex17$all_inv_var_meta_p), y=-log10(mergedcovidgtex17$pval_nominal))

mysnp <- subset(tmp3, rsid == 'rs1105569')
points(x=-log10(mysnp$all_inv_var_meta_p), y=-log10(mysnp$pval_nominal),col='red')

tmp4 <- tmp3 %>% group_by(rsid) %>% dplyr::slice(which.min(pval_nominal)) %>% ungroup()

tmp4 <- tmp4[grep("^rs", tmp4$rsid), ]

myres <- coloc.abf(dataset1=list(pvalues = tmp4$all_inv_var_meta_p,  N= 18152+1145546, type="cc", s=18152/(18152+1145546), snp=tmp4$rsid),
                   dataset2=list(pvalues = tmp4$pval_nominal, N= 147,
                                 type="quant", snp=tmp4$rsid),
                   MAF=tmp4$maf, p1=1e-4,p2=1e-4, p12= 1e-5)
#PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#1.81e-10  5.69e-01  4.45e-11  1.40e-01  2.91e-01 
#[1] "PP abf for shared variant: 29.1%"

