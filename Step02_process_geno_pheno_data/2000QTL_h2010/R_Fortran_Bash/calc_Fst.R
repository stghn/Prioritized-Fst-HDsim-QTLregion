rm(list=ls())

library(pacman)
pacman::p_load(dplyr)

## Define number of samples in sub-population 1, 2, and total
nS1=750;nS2=750;nS1S2=1500

## Reading there subpopulation
## allele frequency of subset animals with low phenotype 
ldata<-read.table("allele.freq.S1pop",head=F)
## allele frequency of subset animals with high phenotype 
hdata<-read.table("allele.freq.S2pop",head=F)
## allele frequency of total animal population(combined low and high phenotype animals) 
data<-read.table("allele.freq.S1S2pop",head=F)

## p and q allele frequency of low phenotype subset
S1<-ldata %>% dplyr::rename(p1=V2,q1=V3)

## p and q allele frequency of high phenotype subset
S2<-hdata %>% dplyr::rename(p2=V2,q2=V3)

## p and q allele frequency of total population
Total<-data %>% dplyr::rename(p=V2,q=V3)


##Hs from low phenotype subset
Hsl<- S1 %>% mutate(Hsl=2*p1*q1)

##Hs from high phenotype subset
Hsh<-S2 %>% mutate(Hsh=2*p2*q2)

## Hs weighted from Hsl and Hsh 
Hsw <-data.frame(Hsw=((Hsl$Hsl*nS1)+(Hsh$Hsh*nS2))/nS1S2)

## Ht from the total population 
Ht<-Total %>% mutate(Ht=2*p*q)


## Calculate Fst score for each loci (600K SNP)
Fst <- cbind(Hsw,Ht) %>% dplyr::mutate(Fst=(Ht-Hsw)/Ht) %>% 
                         dplyr::select(Fst) %>%
                         tibble::rownames_to_column(.,"ID")


write.table(Fst,"Fst_5%_600K",col.names=T, row.names=F,quote=F)

## Define three different thresold to pick most differential signiture selected SNPs
Quantile<- Fst %>% summarize( q99 = quantile(Fst, probs = .99),
                   q97 = quantile(Fst, probs = .97),
                   q95 = quantile(Fst, probs = .95))

## Prioritize SNPs with Fst score greater that 0.99 quantile 
fst99_ID <- Fst %>% dplyr::filter(Fst>Quantile[[1]]) %>% 
                 dplyr::select(ID) 

write.table(fst99_ID,"fst99.ID",col.names=F, row.names=F,quote=F)


