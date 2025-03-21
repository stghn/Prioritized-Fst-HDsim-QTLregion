rm(list=ls())

library(pacman)
pacman::p_load(dplyr,readr,tibble,gtools,ggplot2,ggpubr,forcats)
 
## read QTL file (generated from additive_QTL.R script)
QTL<-readr::read_table("QTL_genetic_percentage") %>% dplyr::select(-p, -q) 

## read Markerinfo file (generated from QTLeffect.sh script)
Marker<-readr::read_table("Markerinfo") %>% dplyr::mutate(addQtl.percent=0)

## rbind QTL and Marker info and # sort by chromosome (chr) and position (pos)
MarkerQTL<- dplyr::bind_rows(Marker,QTL) %>% dplyr::arrange(chr, pos)


## Extract the 2 QTLs from the middle of QTL allele effect distribution
n_rows <- 2
M2QTLs <-  QTL %>% dplyr::arrange(-Alpha) %>% 
           dplyr::slice((n() %/% 2 - ceiling(n_rows/2) + 1):(n() %/% 2 + floor(n_rows/2))) %>%  # Extract the middle rows
           dplyr::arrange(chr) %>% dplyr::pull(Index)

## Extract the 2 top QTLs from the QTL allele effect distribution
T2QTLs <-QTL %>% dplyr::arrange(-Alpha) %>% dplyr::distinct(chr,.keep_all=TRUE) %>% 
                 head(2) %>% dplyr::arrange(chr) %>% dplyr::pull(Index)


## Extract the 2 bottom QTLs from the QTL allele effect distribution
B2QTLs <-QTL %>% dplyr::arrange(Alpha) %>% dplyr::distinct(chr,.keep_all=TRUE) %>%
                 head(2) %>% dplyr::arrange(chr) %>% dplyr::pull(Index)


## Extract QTLs from MarkerQTL file and keep original index of QTL positions in the rowname
qdata <- MarkerQTL %>% tibble::rownames_to_column('rowname') %>%
                       dplyr::filter(grepl('Q',Index)) %>%
                       tibble::column_to_rownames('rowname')


## store a vector of index for top, middle and bottom QTLs based on QTL allelic effect
rowtopQ <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% T2QTLs)))
rowmiddleQ <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% M2QTLs)))
rowbottomQ <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% B2QTLs)))


#######################################################
###   Define the windows of SNPs around 500 QTLs  #####
#######################################################
SNPwindows<-"50SNPwin"
SNPmargin<- 25

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
                                             dplyr::mutate_at(c('Chr'), as.factor) %>%
                                             dplyr::mutate(Chr = factor(Chr, levels = c("chr8","chr25")))
                                             
                                                  

Fst_SNPwin_middle2QTL <- SNPwindows_middle2QTL %>% dplyr::left_join(Fst_600K, by = "Mid") %>% 
						   dplyr::mutate(Chr = paste("chr", chr, sep = "")) %>%
                                                   dplyr::mutate_at(c('Chr'), as.factor) %>%
                                                   dplyr::mutate(Chr = factor(Chr, levels = c("chr13","chr19")))


Fst_SNPwin_bottom2QTL <- SNPwindows_bottom2QTL %>% dplyr::left_join(Fst_600K, by = "Mid") %>%
 						   dplyr::mutate(Chr = paste("chr", chr, sep = "")) %>%
                                                   dplyr::mutate_at(c('Chr'), as.factor) %>%
                                                   dplyr::mutate(Chr = factor(Chr, levels = c("chr4","chr24")))


## calculate mean of Fst score for a SNP windows around top, middle and bottom 2 QTLs
win50SNP_Top2QTL_summary <-  Fst_SNPwin_top2QTL %>% dplyr::group_by(Qid,Chr) %>%
                                                    dplyr::summarize(mean = mean(Fst), .groups = "drop")
win50SNP_Top2QTL_summary

win50SNP_Middle2QTL_summary <-  Fst_SNPwin_middle2QTL %>% dplyr::group_by(Qid,Chr) %>%
                                                    dplyr::summarize(mean = mean(Fst), .groups = "drop")
win50SNP_Middle2QTL_summary

win50SNP_Bottom2QTL_summary <-  Fst_SNPwin_bottom2QTL %>% dplyr::group_by(Qid,Chr) %>%
                                                    dplyr::summarize(mean = mean(Fst), .groups = "drop")
win50SNP_Bottom2QTL_summary


###############################################################################################################
## Plot the distribution of fst score for a window of 50 SNPs around top 2 QTLs based on QTL allelic effect ###
###############################################################################################################

#par(mar=c(7,4,4,2)+0.1)
#png(file=paste0(SNPwindows,"_","Top2QTLs.png"), units="in", width=20, height=10, res=300)

## define fst mean dataframe for 50 SNP-window around top 2 QTLs
datahline<- Fst_SNPwin_top2QTL  %>% dplyr::group_by(Chr) %>%
                             dplyr::summarize(mean_fst = mean(Fst),sd_fst=sd(Fst)) %>%
                             dplyr::mutate(mean_fst=formatC(mean_fst,format = "f",digits = 5)) %>%
                             dplyr::mutate(sd_fst=formatC(sd_fst,format = "f",digits = 5)) %>%
                             dplyr::mutate_if(is.character, ~ as.character(gtools::mixedsort(.)))

mod_data <- Fst_SNPwin_top2QTL %>%
              dplyr::rename(Top_QTL = Qid) %>%
              dplyr::mutate(Top_QTL = factor(Top_QTL, labels = c("QTL163", "QTL452"))) %>%
              dplyr::group_by(Top_QTL) %>%
              dplyr::mutate(Fst = ifelse(Fst < 0, 0, Fst))

# create barplot
topplot <- ggplot(mod_data, aes(x = Mid, y = Fst, fill = Top_QTL)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge(width = 0.75)) +
  facet_wrap(~ ifelse(Chr == "chr8", "QTL163", ifelse(Chr == "chr25", "QTL452", Chr)),
             ncol = 2, scales = "free_x",labeller = as_labeller(c("QTL163" = "QTL163; Large QTL effect",
                         "QTL452" = "QTL452; Large QTL effect"))) +
  scale_fill_manual(values = c("green4", "blue4")) +

  # Fix geom_hline layer
  geom_hline(data = datahline, aes(yintercept = as.numeric(mean_fst)),
             linetype = "dashed", color = "red4",linewidth = 0.6) +

  # Fix geom_label layer
  geom_label(data = datahline,
             aes(x = Inf, y = Inf, label = paste("Average Fst score = ", mean_fst)),
             color = "red4", fontface = "bold.italic", hjust = 1, vjust = 1, size = 7, inherit.aes = FALSE) +

  labs(y = "Fst score", fill = "Large QTL effect") +
  scale_y_continuous(limits = c(0, 0.004)) +

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
        strip.text = element_text(face = "bold", size = 24, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"))



#print(topplot)

#dev.off()
#graphics.off()



##################################################################################################################
## Plot the distribution of fst score for a window of 50 SNPs around middle 2 QTLs based on QTL allelic effect ###
##################################################################################################################

#par(mar=c(7,4,4,2)+0.1)
#png(file=paste0(SNPwindows,"_","Middle2QTLs.png"), units="in", width=20, height=10, res=300)

                                
## define fst mean dataframe for 50 SNP-window around middle 2 QTLs
datahline<- Fst_SNPwin_middle2QTL  %>% dplyr::group_by(Chr) %>%
                             dplyr::summarize(mean_fst = mean(Fst),sd_fst=sd(Fst)) %>%
                             dplyr::mutate(mean_fst=formatC(mean_fst,format = "f",digits = 5)) %>%
                             dplyr::mutate(sd_fst=formatC(sd_fst,format = "f",digits = 5)) %>%
                             dplyr::mutate_if(is.character, ~ as.character(gtools::mixedsort(.)))

mod_data <- Fst_SNPwin_middle2QTL %>%
              dplyr::rename(Middle_QTL = Qid) %>%
              dplyr::mutate(Middle_QTL = factor(Middle_QTL, labels = c("QTL259", "QTL356"))) %>%
              dplyr::group_by(Middle_QTL) %>%
              dplyr::mutate(Fst = ifelse(Fst < 0, 0, Fst))
 
# create barplot
middleplot <- ggplot(mod_data, aes(x = Mid, y = Fst, fill = Middle_QTL)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge(width = 0.75)) +
  facet_wrap(~ ifelse(Chr == "chr13", "QTL259", ifelse(Chr == "chr19", "QTL356", Chr)),
             ncol = 2, scales = "free_x",labeller = as_labeller(c("QTL259" = "QTL259; Medium QTL effect",
                         "QTL356" = "QTL356; Medium QTL effect"))) +
  scale_fill_manual(values = c("green4", "blue4")) +

  # Fix geom_hline layer
  geom_hline(data = datahline, aes(yintercept = as.numeric(mean_fst)),
             linetype = "dashed", color = "red4",linewidth = 0.6) +

  # Fix geom_label layer
  geom_label(data = datahline,
             aes(x = Inf, y = Inf, label = paste("Average Fst score = ", mean_fst)),
             color = "red4", fontface = "bold.italic", hjust = 1, vjust = 1, size = 7, inherit.aes = FALSE) +

  labs(y = "Fst score", fill = "Medium QTL effect") +
  scale_y_continuous(limits = c(0, 0.004)) +

  # Ensure the legend has thicker edges
  guides(fill = guide_legend(override.aes = list(size = 5, colour = "black", linewidth = 2))) +

  theme_bw() +
  theme_classic() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 28, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 24),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "none",  # <<< This removes the legend
        strip.text = element_text(face = "bold", size = 24, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"))

#print (middleplot)

#dev.off()
#graphics.off()



##################################################################################################################
## Plot the distribution of fst score for a window of 50 SNPs around bottom 2 QTLs based on QTL allelic effect ###
##################################################################################################################

#par(mar=c(7,4,4,2)+0.1)
#png(file=paste0(SNPwindows,"_","Bottom2QTLs.png"), units="in", width=20, height=10, res=300)

## define fst mean dataframe for 50 SNP-window around bottom 2 QTLs
datahline <- Fst_SNPwin_bottom2QTL %>%
  dplyr::group_by(Chr) %>%
  dplyr::summarize(mean_fst = mean(Fst), sd_fst = sd(Fst)) %>%
  dplyr::mutate(mean_fst = formatC(mean_fst, format = "f", digits = 5)) %>%
  dplyr::mutate(sd_fst = formatC(sd_fst, format = "f", digits = 5)) %>%
  dplyr::mutate_if(is.character, ~ as.character(gtools::mixedsort(.))) %>%
  dplyr::mutate(Chr = recode_factor(Chr, "chr4" = "QTL85", "chr24" = "QTL436"))


# Rename facet titles in the dataset
mod_data <- Fst_SNPwin_bottom2QTL %>%
  dplyr::rename(Bottom_QTL = Qid) %>%
  dplyr::mutate_at(c('Bottom_QTL'), as.factor) %>%
  dplyr::mutate(Fst = ifelse(Fst < 0, 0, Fst)) %>%
  dplyr::mutate(Chr = recode_factor(Chr, "chr4" = "QTL85", "chr24" = "QTL436")) %>%
  dplyr::mutate(Bottom_QTL = fct_recode(Bottom_QTL, "QTL85" = "Q85", "QTL436" = "Q436"))

# Create barplot with fill instead of colour
bottomplot <- ggplot(mod_data, aes(x = Mid, y = Fst, fill = Bottom_QTL)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge(width = 0.75)) +
  facet_wrap(~Chr, ncol = 2, scales = "free_x",labeller = as_labeller(c("QTL85" = "QTL85; Small QTL effect",
                         "QTL436" = "QTL436; Small QTL effect")))  +
  scale_fill_manual(values = c("blue4", "green4")) +  # Use fill instead of colour
  geom_hline(data = datahline, aes(yintercept = as.numeric(mean_fst)), linetype = "dashed", color = "red4",linewidth=0.6) + # Fix geom_hline layer
  geom_label(data = datahline, aes(x = Inf, y = Inf, label = paste("Average Fst score = ", mean_fst)), # Fix geom_label layer
             color = "red4", fontface = "bold.italic", hjust = 1, vjust = 1, size = 7,
             inherit.aes = FALSE) +  # Prevents geom_label() from inheriting fill = Bottom_QTL
  labs(y = "Fst score", fill = "Small QTL effect") +
  scale_y_continuous(limits = c(0, 0.004)) +
  # Ensure the legend has thicker edges
  guides(fill = guide_legend(override.aes = list(size = 5, colour = "black", linewidth = 2),reverse=TRUE)) +
#  guides(fill = guide_legend(override.aes = list(size = 10), reverse = TRUE)) +
  theme_bw() +
  theme_classic() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 28, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 24),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "none",  # <<< This removes the legend
        strip.text = element_text(face = "bold", size = 24, colour = "white"),
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

ggplot2::ggsave(paste0(SNPwindows,"_","around2QTLs_QTLeffect.png"), unit="in",width = 20, height = 20, dpi = 300)


