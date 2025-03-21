rm(list=ls())

library(pacman)
pacman::p_load(dplyr,readr,tibble,gtools,gdata)

## read QTL file (generated from additive_QTL.R script)
QTL<-readr::read_table("QTL_genetic_percentage") %>% dplyr::select(-p, -q) 

## read Markerinfo file (generated from QTLeffect.sh script)
Marker<-readr::read_table("Markerinfo") %>% dplyr::mutate(addQtl.percent=0) 

## rbind QTL and Marker info and # sort by chromosome (chr) and position (pos)
MarkerQTL<- dplyr::bind_rows(Marker,QTL) %>% dplyr::arrange(chr, pos)

## Extract QTLs from MarkerQTL file and keep original index of QTL positions in the rowname
qdata <- MarkerQTL %>% tibble::rownames_to_column('rowname') %>% 
                       dplyr::filter(grepl('Q',Index)) %>%
                       tibble::column_to_rownames('rowname')

##-----------------------------------------------------

## store a vector of index for all QTL positions
rowq <- as.numeric(row.names(qdata))

###########################################################
###   Define the windows of 200 SNPs around 500 QTLs  #####
###########################################################
SNPwindows<-"200win"
SNPmargin<- 100

## function to extract 15 SNPs before and after each QTL index 
extract_rows <- function(df, index_row, margin=SNPmargin) {
  df %>%
    filter(row_number() >= (index_row - margin) & row_number() <= (index_row + margin))
}

## using lapply function to extract windows of SNP around each QTL index 
allQTLs <- lapply(rowq, extract_rows, df = MarkerQTL)

## check the length of vectors within each nested list
vector_lengths <- sapply(allQTLs, length)
is_same_length <- all(vector_lengths == vector_lengths[1])
is_same_length


### rename column Index to Mid across all list dataframe
allQTLs <-lapply(allQTLs,function(x) x %>% dplyr::rename(Mid=Index))

## adding a column to every dataframe in a list class and then merge all list elements into one dataframe
allQTLs <- Map(cbind, allQTLs, Qid=qdata$Index) %>%  dplyr::bind_rows()

##----------------------------------------------------------

## Extract only Marker ids from Mid column
winSNP <- as_tibble(allQTLs %>% dplyr::filter(startsWith(Mid, "M"))) %>% dplyr::select(Mid,chr,pos,Qid)
head(winSNP)

## read Fst file for 600K SNP markers and modify to get desired Fst output 
Fst_600K<-readr::read_table("Fst_5%_600K") %>% dplyr::mutate(Mid=paste0("M", ID)) %>%
          dplyr:: mutate(Fst = ifelse(Fst < 0, 0, Fst)) %>% dplyr::select(Mid, Fst)
head(Fst_600K)


## Set fst threshold from 600K SNP to compare with average fst from each SNP window surrounding QTLs
fst_thr_600K <- quantile(Fst_600K$Fst,probs=0.25)
fst_thr="q25"

## join windows of SNPs around each QTL with Fst scores of all SNP marker panel 
Fst_winSNP <- winSNP %>% dplyr::left_join(Fst_600K, by = "Mid") 

## Evaluate the number of QTLs pass and qualify 
## if average fst from each SNP window surrounding QTLs is greater than defined fst threshold from 600K SNP panel
Num_QTLs_pass_threshold <-  Fst_winSNP %>% dplyr::group_by(Qid) %>%
                            dplyr::summarize(meanfst = mean(Fst)) %>%        
                            dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 6))) %>%
                            dplyr::arrange(order(gtools::mixedorder(Qid))) %>%
                            dplyr::mutate(QTL_pass = case_when(
                                                      meanfst > fst_thr_600K ~ "Pass",
                                                      meanfst <= fst_thr_600K ~ "Fail")) %>%
			    dplyr::filter(QTL_pass == "Pass")
head(Num_QTLs_pass_threshold) 

## merging passed QTLs from "Num_QTLs_pass_threshold" data with percentage of genetic variance explained for each QTL in "QTL" data 
passed_QTLs_summary <- Num_QTLs_pass_threshold  %>% dplyr::rename(Index=Qid)  %>% 
                                                      dplyr::left_join(QTL, by = "Index") %>% 
                                                      dplyr::select(Index,chr,pos,Alpha,addQtl.percent) %>%
						      dplyr::summarize(Num_QTL = n(), Percent_Gvar=sum(addQtl.percent))

## Using the write.fwf from the gdata package if you want to create a Fixed Width File.
## The col names don't seem to be saved to the correct position,
## so to work around this you can format them to the width you want the columns to be using formatC

colnames(passed_QTLs_summary) <- formatC(colnames(passed_QTLs_summary), format = "s",width = 12, flag = "0")

gdata::write.fwf(as.data.frame(passed_QTLs_summary),paste0("Passed_QTLs_", SNPwindows,"_",fst_thr), width = rep(12,ncol(passed_QTLs_summary)),colnames = T,rownames = F,quote = F, sep = " ")


## join windows of SNPs around qualified QTLs where passed the defined Fst threshold for all 600K SNP markers
Passed_QTLs <- Fst_winSNP %>% dplyr::right_join(Num_QTLs_pass_threshold, by = "Qid")

## Prioritize random 12 SNPs surrounding qualified QTLs based on Fst scores
set.seed(123)
Random12SNPsFst <- Passed_QTLs %>% dplyr::group_by(Qid) %>%
                      dplyr::slice_sample(n = 12) %>%
                      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 6))) %>%
                      dplyr::arrange(order(gtools::mixedorder(Mid))) %>%
                      dplyr::arrange(order(gtools::mixedorder(Qid)))
head(Random12SNPsFst)

readr::write_csv(Random12SNPsFst, paste0(SNPwindows,"_","random12snpselection.csv"))

## determine the index of randomly selected SNPs based on Fst scores for qualified QTLs
Random12SNPsFst <- Random12SNPsFst %>% dplyr::mutate(Mid = sub("M", "", Mid)) %>%
                                       ungroup() %>% dplyr::select(Mid) %>%
                                     dplyr::distinct(Mid, .keep_all=TRUE)

write.table(Random12SNPsFst,paste0("Rnd12SNPsFst_", SNPwindows),col.names=F, row.names=F,quote=F)



