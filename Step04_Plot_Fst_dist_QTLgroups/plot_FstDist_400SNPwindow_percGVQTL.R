rm(list=ls())

setwd("C:/Users/Sajjad.Toghiani/Downloads/WQP_h2_rep1/WQP_h2_rep1/h2040/W1Q1P2_h2040_rep1")

library(pacman)
pacman::p_load(dplyr,readr,tibble,gtools,ggplot2,ggpubr,forcats)

## read QTL file (generated from additive_QTL.R script)
QTL<-readr::read_table("QTL_genetic_percentage") %>% dplyr::select(-p, -q)

## read Markerinfo file (generated from QTLeffect.sh script)
Marker<-readr::read_table("Markerinfo") %>% dplyr::mutate(addQtl.percent=0)

## rbind QTL and Marker info and # sort by chromosome (chr) and position (pos)
MarkerQTL<- dplyr::bind_rows(Marker,QTL) %>% dplyr::arrange(chr, pos)

## read prioritize random 12 SNPs from qualified SNP windows (generated from calc_NumQTLsPassed_FstThreshold.R script)
Rnd12_400win <-readr::read_csv("400win_random12snpselection.csv",show_col_types = FALSE) 

## Extract the 2 middle QTLs from %GV explained by QTL distribution
n_rows <-2
M2QTLs <-  QTL %>% dplyr::arrange(-addQtl.percent) %>% 
  dplyr::slice((n() %/% 2 - ceiling(n_rows/2) + 1):(n() %/% 2 + floor(n_rows/2))) %>% # Extract the middle rows
  dplyr::arrange(chr) %>% dplyr::pull(Index) %>% dplyr::nth(2)

## Extract the 2 top QTLs from %GV explained by QTL distribution
T2QTLs <-QTL %>% dplyr::arrange(-addQtl.percent) %>% dplyr::distinct(chr,.keep_all=TRUE) %>% 
  head(2) %>% dplyr::arrange(chr) %>% dplyr::pull(Index) %>% dplyr::nth(2)

## Extract the 2 bottom QTLs from %GV explained by QTL distribution
B2QTLs <-QTL %>% dplyr::arrange(addQtl.percent) %>% dplyr::distinct(chr,.keep_all=TRUE) %>%
  head(2) %>% dplyr::arrange(chr) %>% dplyr::pull(Index) %>% dplyr::nth(1)


## Extract QTLs from MarkerQTL file and keep original index of QTL positions in the rowname
qdata <- MarkerQTL %>% tibble::rownames_to_column('rowname') %>%
  dplyr::filter(grepl('Q',Index)) %>%
  tibble::column_to_rownames('rowname')


## store a vector of index for top, middle and bottom QTLs based on %GV explained by QTL 
rowtopQ <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% T2QTLs)))
rowmiddleQ <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% M2QTLs)))
rowbottomQ <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% B2QTLs)))


#######################################################
###   Define the windows of SNPs around 500 QTLs  #####
#######################################################
SNPwindows<-"400SNPwin"
SNPmargin<- 200

## function to extract 15 SNPs before and after each QTL index
extract_rows <- function(df, index_row, margin=SNPmargin) {
  df %>%
    filter(row_number() >= (index_row - margin) & row_number() <= (index_row + margin))
}

## using lapply function to extract windows of SNP around each QTL index
top2QTLdf <- lapply(rowtopQ, extract_rows, df = MarkerQTL); names(top2QTLdf) <- T2QTLs
middle2QTLdf <- lapply(rowmiddleQ, extract_rows, df = MarkerQTL); names(middle2QTLdf) <- M2QTLs
bottom2QTLdf <- lapply(rowbottomQ, extract_rows, df = MarkerQTL); names(bottom2QTLdf) <- B2QTLs


### rename column Index to Mid across all list dataframe
top2QTLdf <-lapply(top2QTLdf,function(x) x %>% dplyr::rename(Mid=Index))
middle2QTLdf <-lapply(middle2QTLdf,function(x) x %>% dplyr::rename(Mid=Index))
bottom2QTLdf <-lapply(bottom2QTLdf,function(x) x %>% dplyr::rename(Mid=Index))

## adding a column to every dataframe in a list class and then merge all list elements into one dataframe
top2QTLdf <- Map(cbind, top2QTLdf, Qid=names(top2QTLdf)) %>%  dplyr::bind_rows()
middle2QTLdf <- Map(cbind, middle2QTLdf, Qid=names(middle2QTLdf)) %>%  dplyr::bind_rows()
bottom2QTLdf <- Map(cbind, bottom2QTLdf, Qid=names(bottom2QTLdf)) %>%  dplyr::bind_rows()

##----------------------------------------------------------

## Extract only Marker ids from Mid column
SNPwindows_top2QTL <- as_tibble(top2QTLdf %>% dplyr::filter(startsWith(Mid, "M")) %>% dplyr::select(Mid,chr,pos,Qid))
head(SNPwindows_top2QTL)

SNPwindows_middle2QTL <- as_tibble(middle2QTLdf %>% dplyr::filter(startsWith(Mid, "M")) %>% dplyr::select(Mid,chr,pos,Qid))
head(SNPwindows_middle2QTL)

SNPwindows_bottom2QTL <- as_tibble(bottom2QTLdf %>% dplyr::filter(startsWith(Mid, "M")) %>% dplyr::select(Mid,chr,pos,Qid))
head(SNPwindows_bottom2QTL)


## read Fst files for 600K SNP markers
Fst_600K<-readr::read_table("Fst_5%_600K") %>% dplyr::mutate(Mid=paste0("M", ID)) %>%
  dplyr::mutate(Fst = ifelse(Fst < 0, 0, Fst)) %>%  
  dplyr::select(Mid, Fst)
head(Fst_600K)

## calculate mean and sd of Fst score across all 600K SNPs
Fst600K_summary <-  Fst_600K %>%  summarize(Min=min(Fst), Q25=quantile(Fst, 0.25), Median=median(Fst), 
                                            Mean = mean(Fst), Q75=quantile(Fst, 0.75), Max=max(Fst) )
Fst600K_summary

## join a SNP windows around top, middle, and bottom 2 QTLs with Fst scores of all SNP marker panel
Fst_SNPwin_top2QTL <- SNPwindows_top2QTL %>% dplyr::left_join(Fst_600K, by = "Mid") %>%
  dplyr::mutate(Chr = paste("chr", chr, sep = "")) %>%
  dplyr::mutate_at(c('Chr','Qid'), as.factor) %>%
  dplyr::mutate(Chr = factor(Chr, levels = c("chr16")))


Fst_SNPwin_middle2QTL <- SNPwindows_middle2QTL %>% dplyr::left_join(Fst_600K, by = "Mid") %>% 
  dplyr::mutate(Chr = paste("chr", chr, sep = "")) %>%
  dplyr::mutate_at(c('Chr','Qid'), as.factor) %>%
  dplyr::mutate(Chr = factor(Chr, levels = c("chr21")))

Fst_SNPwin_bottom2QTL <- SNPwindows_bottom2QTL %>% dplyr::left_join(Fst_600K, by = "Mid") %>%
  dplyr::mutate(Chr = paste("chr", chr, sep = "")) %>%
  dplyr::mutate_at(c('Chr','Qid'), as.factor) %>%
  dplyr::mutate(Chr = factor(Chr, levels = c("chr4")))


## calculate mean of Fst score for a SNP windows around top, middle and bottom 2 QTLs
win400SNP_Top2QTL_summary <-  Fst_SNPwin_top2QTL %>% dplyr::group_by(Qid,Chr) %>% 
  dplyr::summarize(mean = mean(Fst), .groups = "drop")                    
win400SNP_Top2QTL_summary

win400SNP_Middle2QTL_summary <-  Fst_SNPwin_middle2QTL %>% dplyr::group_by(Qid,Chr) %>%
  dplyr::summarize(mean = mean(Fst), .groups = "drop")
win400SNP_Middle2QTL_summary

win400SNP_Bottom2QTL_summary <-  Fst_SNPwin_bottom2QTL %>% dplyr::group_by(Qid,Chr) %>%
  dplyr::summarize(mean = mean(Fst), .groups = "drop")
win400SNP_Bottom2QTL_summary


##################################################################################################################
## Plot the distribution of fst score for a window of 50 SNPs around top 2 QTLs based on %GV explained by QTL  ###
##################################################################################################################

#par(mar=c(7,4,4,2)+0.1)
#png(file=paste0(SNPwindows,"_","Top2QTLs.png"), units="in", width=20, height=10, res=300)

## define fst mean dataframe for 400 SNP-window around top 2 QTLs
datahline<- Fst_SNPwin_top2QTL  %>% dplyr::group_by(Chr) %>% 
  dplyr::summarize(mean_fst = mean(Fst),sd_fst=sd(Fst)) %>%
  dplyr::mutate(mean_fst=formatC(mean_fst,format = "f",digits = 5)) %>%
  dplyr::mutate(sd_fst=formatC(sd_fst,format = "f",digits = 5)) %>%
  dplyr::mutate_if(is.character, ~ as.character(gtools::mixedsort(.)))

# Select the desired Mids from Rnd12_400win
selected_Mids<- Rnd12_400win %>% dplyr::filter(grepl('Q313',Qid)) %>%
  dplyr::pull(Mid) 

# Add an indicator column to Fst_SNPwin_top2QTL
Fst_SNPwin_top2QTL <- Fst_SNPwin_top2QTL %>%
  dplyr::mutate(indicator = ifelse(Mid %in% selected_Mids, "Selected", "NotSelected"))


mod_data <- Fst_SNPwin_top2QTL %>%
  dplyr::rename(Top_QTL = Qid) %>%
  dplyr::mutate(Top_QTL = factor(Top_QTL, labels = c("QTL313"))) %>%
  dplyr::group_by(Top_QTL) %>%
  dplyr::mutate(Fst = ifelse(Fst < 0, 0, Fst)) 

# Find the middle row of the dataframe
middle_row <- mod_data[(nrow(mod_data) %/% 2) + 1, ]
middle_row

## create barplot
topplot <- ggplot(mod_data, aes(x = Mid, y = Fst, fill = indicator)) +  
  geom_bar(stat = "identity", width = 0.4, position = position_dodge(width = 0.75)) +
  facet_wrap(~ ifelse(Chr == "chr16", "QTL313", Chr), 
             ncol = 2, scales = "free_x",labeller = as_labeller(c("QTL313" = "QTL313; Large QTL variance"))) +
  # Define colors for "Selected" and "NotSelected"
  scale_fill_manual(values = c("Selected" = "red", "NotSelected" = "green4")) +  
  
  # Fix geom_hline layer
  geom_hline(data = datahline, aes(yintercept = as.numeric(mean_fst)), 
             linetype = "dashed", color = "red4",linewidth = 1) +
  
  # Fix geom_label layer
  geom_label(data = datahline, 
             aes(x = Inf, y = Inf, label = paste("Average Fst score = ", mean_fst)), 
             color = "red4", fontface = "bold.italic", hjust = 1, vjust = 1, size = 10, inherit.aes = FALSE) +
  
  labs(y = "Fst score") +
  scale_y_continuous(limits = c(0, 0.004)) +
  geom_vline(xintercept = middle_row[["Mid"]], linetype = "dashed", color = "gray35", linewidth = 1) +  # Add vertical dashed line at median 
  # Ensure the legend has thicker edges
  guides(fill = guide_legend(override.aes = list(size = 5, colour = "black", linewidth = 2))) +
  
  theme_bw() +
  theme_classic() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 28, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 26),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "none",  # <<< This removes the legend
        strip.text = element_text(face = "bold", size = 32, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"))

#print(topplot)

#dev.off()
#graphics.off()



##################################################################################################################
## Plot the distribution of fst score for a window of 50 SNPs around middle 2 QTLs based on %GV explained by QTL ###
##################################################################################################################

#par(mar=c(7,4,4,2)+0.1)
#png(file=paste0(SNPwindows,"_","Middle2QTLs.png"), units="in", width=20, height=10, res=300)

## define fst mean dataframe for 400 SNP-window around middle 2 QTLs
datahline<- Fst_SNPwin_middle2QTL  %>% dplyr::group_by(Chr) %>%
  dplyr::summarize(mean_fst = mean(Fst),sd_fst=sd(Fst)) %>%
  dplyr::mutate(mean_fst=formatC(mean_fst,format = "f",digits = 5)) %>%
  dplyr::mutate(sd_fst=formatC(sd_fst,format = "f",digits = 5)) %>%
  dplyr::mutate_if(is.character, ~ as.character(gtools::mixedsort(.)))

# Select the desired Mids from Rnd12_400win
selected_Mids<- Rnd12_400win %>% dplyr::filter(grepl('Q380',Qid)) %>%
  dplyr::pull(Mid) 

# Add an indicator column to Fst_SNPwin_middle2QTL
Fst_SNPwin_middle2QTL <- Fst_SNPwin_middle2QTL %>%
  dplyr::mutate(indicator = ifelse(Mid %in% selected_Mids, "Selected", "NotSelected"))


mod_data <- Fst_SNPwin_middle2QTL %>%
  dplyr::rename(Top_QTL = Qid) %>%
  dplyr::mutate(Top_QTL = factor(Top_QTL, labels = c("QTL380"))) %>%
  dplyr::group_by(Top_QTL) %>%
  dplyr::mutate(Fst = ifelse(Fst < 0, 0, Fst)) 

# Find the middle row of the dataframe
middle_row <- mod_data[(nrow(mod_data) %/% 2) + 1, ]
middle_row

## create barplot
middleplot <- ggplot(mod_data, aes(x = Mid, y = Fst, fill = indicator)) +  
  geom_bar(stat = "identity", width = 0.4, position = position_dodge(width = 0.75)) +
  facet_wrap(~ ifelse(Chr == "chr21", "QTL380", Chr), 
             ncol = 2, scales = "free_x",labeller = as_labeller(c("QTL380" = "QTL380; Medium QTL variance"))) +
  # Define colors for "Selected" and "NotSelected"
  scale_fill_manual(values = c("Selected" = "red", "NotSelected" = "green4")) +  
  
  # Fix geom_hline layer
  geom_hline(data = datahline, aes(yintercept = as.numeric(mean_fst)), 
             linetype = "dashed", color = "red4",linewidth = 1) +
  
  # Fix geom_label layer
  geom_label(data = datahline, 
             aes(x = Inf, y = Inf, label = paste("Average Fst score = ", mean_fst)), 
             color = "red4", fontface = "bold.italic", hjust = 1, vjust = 1, size = 10, inherit.aes = FALSE) +
  
  labs(y = "Fst score") +
  scale_y_continuous(limits = c(0, 0.004)) +
  geom_vline(xintercept = middle_row[["Mid"]], linetype = "dashed", color = "gray35", linewidth = 1) +  # Add vertical dashed line at median 
  # Ensure the legend has thicker edges
  guides(fill = guide_legend(override.aes = list(size = 5, colour = "black", linewidth = 2))) +
  
  theme_bw() +
  theme_classic() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 28, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 26),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "none",  # <<< This removes the legend
        strip.text = element_text(face = "bold", size = 32, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"))


#print (middleplot)

#dev.off()
#graphics.off()



##################################################################################################################
## Plot the distribution of fst score for a window of 50 SNPs around bottom 2 QTLs based on %GV explained by QTL ###
##################################################################################################################

#par(mar=c(7,4,4,2)+0.1)
#png(file=paste0(SNPwindows,"_","Bottom2QTLs.png"), units="in", width=20, height=10, res=300)

## define fst mean dataframe for 400 SNP-window around bottom 2 QTLs
datahline <- Fst_SNPwin_bottom2QTL %>%
  dplyr::group_by(Chr) %>%
  dplyr::summarize(mean_fst = mean(Fst), sd_fst = sd(Fst)) %>%
  dplyr::mutate(mean_fst = formatC(mean_fst, format = "f", digits = 5)) %>%
  dplyr::mutate(sd_fst = formatC(sd_fst, format = "f", digits = 5)) %>%
  dplyr::mutate_if(is.character, ~ as.character(gtools::mixedsort(.))) 

# Select the desired Mids from Rnd12_400win
selected_Mids<- Rnd12_400win %>% dplyr::filter(grepl('Q85',Qid)) %>%
  dplyr::pull(Mid) 

# Add an indicator column to Fst_SNPwin_bottom2QTL
Fst_SNPwin_bottom2QTL <- Fst_SNPwin_bottom2QTL %>%
  dplyr::mutate(indicator = ifelse(Mid %in% selected_Mids, "Selected", "NotSelected"))

# Rename facet titles in the dataset
mod_data <- Fst_SNPwin_bottom2QTL %>%
  dplyr::rename(Bottom_QTL = Qid) %>%
  dplyr::mutate_at(c('Bottom_QTL'), as.factor) %>%
  dplyr::mutate(Fst = ifelse(Fst < 0, 0, Fst)) 

# Find the middle row of the dataframe
middle_row <- mod_data[(nrow(mod_data) %/% 2) + 1, ]
middle_row

## create barplot
bottomplot <- ggplot(mod_data, aes(x = Mid, y = Fst, fill = indicator)) +  
  geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.75)) +
  facet_wrap(~ ifelse(Chr == "chr4", "QTL85", Chr), 
             ncol = 2, scales = "free_x",labeller = as_labeller(c("QTL85" = "QTL85; Small QTL variance"))) +
  # Define colors for "Selected" and "NotSelected"
  scale_fill_manual(values = c("Selected" = "red", "NotSelected" = "green4")) +  
  
  # Fix geom_hline layer
  geom_hline(data = datahline, aes(yintercept = as.numeric(mean_fst)), 
             linetype = "dashed", color = "red4",linewidth = 1) +
  
  # Fix geom_label layer
  geom_label(data = datahline, 
             aes(x = Inf, y = Inf, label = paste("Average Fst score = ", mean_fst)), 
             color = "red4", fontface = "bold.italic", hjust = 1, vjust = 1, size = 10, inherit.aes = FALSE) +
  
  labs(y = "Fst score") +
  scale_y_continuous(limits = c(0, 0.004)) +
  geom_vline(xintercept = middle_row[["Mid"]], linetype = "dashed", color = "gray35", linewidth = 1) +  # Add vertical dashed line at median 
  # Ensure the legend has thicker edges
  guides(fill = guide_legend(override.aes = list(size = 5, colour = "black", linewidth = 2))) +
  theme_bw() +
  theme_classic() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 28, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 26),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "none",  # <<< This removes the legend
        strip.text = element_text(face = "bold", size = 32, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"))

#print(bottomplot)

#dev.off()
#graphics.off()



#########################################################################################################
####################### To arrange multiple ggplots on one single page ###################################
##########################################################################################################
figure<-ggpubr::ggarrange(topplot, middleplot, bottomplot, ncol = 1, nrow = 3,labels = c("A", "B", "C"), 
                          align = "v", font.label = list(size = 30, color = "black"))

#annotate_figure(figure, top = text_grob("Fst score distribution for a window of 50 SNPs surrounding top (A), middle (B), and bottom (C) 2 QTLs selected base on %GV explained by QTL ", color = "black", size=14,face = "bold"))

ggplot2::ggsave(paste0(SNPwindows,"_","around1QTL_percGVQTL_h2040.png"), unit="in",width = 20, height = 20, dpi = 300)

