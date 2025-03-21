rm(list=ls())

library(pacman)
pacman::p_load(dplyr,readr,tibble,gtools,gdata,data.table)

## read QTL file (generated from additive_QTL.R script)
QTL<-readr::read_table("QTL_genetic_percentage") %>% dplyr::select(-p, -q)

## read Markerinfo file (generated from QTLeffect.sh script)
Marker<-readr::read_table("Markerinfo") %>% dplyr::mutate(addQtl.percent=0)

## rbind QTL and Marker info and # sort by chromosome (chr) and position (pos)
MarkerQTL<- dplyr::bind_rows(Marker,QTL) %>% dplyr::arrange(chr, pos)

summary(QTL$Alpha)
length(QTL$Alpha)

# Calculate mean and SD for the top 95th percentile group
top95_stats <- QTL %>%
         dplyr::filter(Alpha >= quantile(Alpha, probs = 0.95)) %>%
         dplyr::summarise(mean_Alpha = mean(Alpha, na.rm = TRUE), 
                          sd_Alpha = sd(Alpha, na.rm = TRUE),
                          n_selected = n() )

# Calculate mean and SD for the middle 25th-75th percentile group
middle2575_stats <- QTL %>%
                  dplyr::filter(Alpha >= quantile(Alpha, probs = 0.25) &  
                          Alpha <= quantile(Alpha, probs = 0.75)) %>%
                  dplyr::summarise(mean_Alpha = mean(Alpha, na.rm = TRUE), 
                                   sd_Alpha = sd(Alpha, na.rm = TRUE),
                                   n_selected = n() )

# Calculate mean and SD for the bottom 5th percentile group
bottom5_stats <- QTL %>%
  		dplyr::filter(Alpha <= quantile(Alpha, probs = 0.05)) %>%
  		dplyr::summarise(mean_Alpha = mean(Alpha, na.rm = TRUE), 
                                 sd_Alpha = sd(Alpha, na.rm = TRUE),
                                 n_selected = n() )

# Set seed for reproducibility (ensures the same 25 QTLs are selected every time)
set.seed(123)

# Randomly select 25 QTLs from within the 25th-75th percentile range of addQtl.percent
sample2575Q_stats <-  QTL %>%
           #Filter within 25th-75th percentile
          dplyr::filter(Alpha >= quantile(Alpha, probs = 0.25) &  
                 Alpha <= quantile(Alpha, probs = 0.75)) %>%  
          dplyr::sample_n(size = 25, replace = FALSE) %>%  # Randomly select 25 QTLs
          # Calculate mean and standard deviation of Alpha for the selected 25 QTLs %>%
          dplyr::summarise(mean_Alpha = mean(Alpha, na.rm = TRUE), 
                           sd_Alpha = sd(Alpha, na.rm = TRUE),
                           n_selected = n() )


# Combine the tibbles into one table
combined_stats <- dplyr::bind_rows(top95_stats %>% mutate(group = "Top 95th percentile"),
			           middle2575_stats %>% mutate(group = "Middle 25th-75th percentile"),
                                   sample2575Q_stats %>% mutate(group = "Random 25 QTLs (25th-75th percentile)"),
				   bottom5_stats %>% mutate(group = "Bottom 5th percentile"))

# Print the combined table
print(combined_stats)

readr::write_csv(combined_stats, "QTLeffect_stats_QTLgroups_W1Q1P2.csv")


# Extract QTLs where Alpha is in the top 5% (>= 95th percentile)
top95qtl <- QTL %>%
         dplyr::filter(Alpha >= quantile(Alpha, probs = 0.95)) %>%  # Filter top 5%
        dplyr::arrange(chr) %>% #Arrange by chromosome
        dplyr::pull(Index) %>%  #Extract the Index column
        gtools::mixedsort()  # Sorting alphanumeric vector

# Get the number of QTLs in the top 95% group
length(top95qtl)

# Set seed for reproducibility (ensures the same 25 QTLs are selected every time)
set.seed(123)

# Randomly select 25 QTLs from within the 25th-75th percentile range of Alpha
sample2575qtl <-  QTL %>%
          dplyr::filter(Alpha >= quantile(Alpha, probs = 0.25) &  
                 Alpha <= quantile(Alpha, probs = 0.75)) %>%  # Filter within  25th-75th percentile
          dplyr::sample_n(size = 25, replace = FALSE) %>%  # Randomly select 25 QTLs
          dplyr::arrange(chr) %>%  #Arrange by chromosome
          dplyr::pull(Index) %>%  #Extract the Index column
          gtools::mixedsort() # Sorting alphanumeric vector

length(sample2575qtl)

# Extract all QTLs within the 25th-75th  percentile range (no random sampling)
middle2575qtl <- QTL %>%
             dplyr::filter(Alpha >= quantile(Alpha, probs = 0.25) &  
                           Alpha <= quantile(Alpha, probs = 0.75)) %>%  # Filter within  25th-75th percentile
             dplyr::arrange(chr) %>% # Arrange by chromosome
             pull(Index) %>%  # Extract the Index column
             gtools::mixedsort() # Sorting alphanumeric vector

# Get the number of QTLs in the middle 25th-75th percentile group
length(middle2575qtl)

# Extract QTLs where Alpha is in the bottom 5% (<= 5th percentile)
bottom5qtl <- QTL %>%
            filter(Alpha <= quantile(Alpha, probs = 0.05)) %>%  #Filter bottom 5%
            dplyr::arrange(chr) %>% # Arrange by chromosome
            dplyr::pull(Index) %>%  # Extract the Index column
            gtools::mixedsort() # Sorting alphanumeric vector

# Get the number of QTLs in the bottom 5% group
length(bottom5qtl)


## Extract QTLs from MarkerQTL file and keep original index of QTL positions in the rowname
qdata <- MarkerQTL %>% tibble::rownames_to_column('rowname') %>%
                       dplyr::filter(grepl('Q',Index)) %>%
                       tibble::column_to_rownames('rowname')

## store a vector of index for top, middle and bottom QTL groups 
rowtop95Q <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% top95qtl)))
#row2575Q <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% middle2575qtl)))
row2575Q <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% sample2575qtl)))
rowbottom5Q <- as.numeric(rownames(qdata %>% dplyr::filter(Index %in% bottom5qtl)))


###################################################################
###   Define the windows of SNPs around each QTL in 3 groups  #####
###################################################################
SNPwindows<-"50SNPwin"
SNPmargin<- 25

## function to extract SNPs before and after each QTL index
extract_rows <- function(df, index_row, margin=SNPmargin) {
  df %>%
    filter(row_number() >= (index_row - margin) & row_number() <= (index_row + margin))
}

## using lapply function to extract windows of SNP around each QTL index
topQTLdf <- lapply(rowtop95Q, extract_rows, df = MarkerQTL); names(topQTLdf) <- top95qtl
#middleQTLdf <- lapply(row2575Q, extract_rows, df = MarkerQTL); names(middleQTLdf) <- middle2575qtl
middleQTLdf <- lapply(row2575Q, extract_rows, df = MarkerQTL); names(middleQTLdf) <- sample2575qtl
bottomQTLdf <- lapply(rowbottom5Q, extract_rows, df = MarkerQTL); names(bottomQTLdf) <- bottom5qtl


### rename column Index to Mid across all list dataframe
topQTLdf <-lapply(topQTLdf,function(x) x %>% dplyr::rename(Mid=Index))
middleQTLdf <-lapply(middleQTLdf,function(x) x %>% dplyr::rename(Mid=Index))
bottomQTLdf <-lapply(bottomQTLdf,function(x) x %>% dplyr::rename(Mid=Index))

## adding a column to every dataframe in a list class and then merge all list elements into one dataframe
topQTLdf <- Map(cbind, topQTLdf, Qid=names(topQTLdf)) %>%  dplyr::bind_rows()
middleQTLdf <- Map(cbind, middleQTLdf, Qid=names(middleQTLdf)) %>%  dplyr::bind_rows()
bottomQTLdf <- Map(cbind, bottomQTLdf, Qid=names(bottomQTLdf)) %>%  dplyr::bind_rows()


## Extract only Marker ids from Mid column
SNPwindows_topQTL <- as_tibble(topQTLdf %>% 
                               dplyr::filter(startsWith(Mid, "M")) %>% 
                               dplyr::select(Mid,chr,pos,Qid))
head(SNPwindows_topQTL)

SNPwindows_middleQTL <- as_tibble(middleQTLdf %>% 
                                  dplyr::filter(startsWith(Mid, "M")) %>% 
                                  dplyr::select(Mid,chr,pos,Qid))
head(SNPwindows_middleQTL)

SNPwindows_bottomQTL <- as_tibble(bottomQTLdf %>% 
                                  dplyr::filter(startsWith(Mid, "M")) %>% 
                                  dplyr::select(Mid,chr,pos,Qid))
head(SNPwindows_bottomQTL)

## read Fst files for 600K SNP markers
Fst_600K<-readr::read_table("Fst_5%_600K") %>% dplyr::mutate(Mid=paste0("M", ID)) %>%
                                               dplyr::mutate(Fst = ifelse(Fst < 0, 0, Fst)) %>%
                                               dplyr::select(Mid, Fst)
head(Fst_600K)

## calculate mean and sd of Fst score across all 600K SNPs
Fst600K_summary <-  Fst_600K %>%   dplyr::summarise(min_fst=min(Fst), Q25_fst=quantile(Fst, 0.25), median_fst=median(Fst),
                   				    mean_fst = mean(Fst,na.rm = TRUE), Q75_fst=quantile(Fst, 0.75), max_fst=max(Fst),
                                                    sd_fst = sd(Fst, na.rm = TRUE), N_marker = n() )
Fst600K_summary

readr::write_csv(Fst600K_summary, "Fst600K_summary_h2040.csv")

## join a SNP windows around top, middle, and bottom  QTLs with Fst scores of all SNP marker panel
Fst_SNPwin_topQTL <- SNPwindows_topQTL %>% dplyr::left_join(Fst_600K, by = "Mid") 
Fst_SNPwin_middleQTL <- SNPwindows_middleQTL %>% dplyr::left_join(Fst_600K, by = "Mid")
Fst_SNPwin_bottomQTL <- SNPwindows_bottomQTL %>% dplyr::left_join(Fst_600K, by = "Mid") 


## calculate mean and sd of Fst score for a SNP windows around top, middle and bottom 25 QTLs
win50SNP_TopQTL_stats <-  Fst_SNPwin_topQTL %>% dplyr::group_by(Qid) %>%
                                                  dplyr::summarize(fst_win_mean = mean(Fst,na.rm = TRUE),fst_win_sd = sd(Fst, na.rm = TRUE),N_win = n()) %>%
                                                  dplyr::summarise(average_fst_win_mean = mean(fst_win_mean,na.rm = TRUE),
                                                                   average_fst_win_sd = mean(fst_win_sd, na.rm = TRUE), average_N_win = round(mean(N_win,na.rm = TRUE)))

win50SNP_MiddleQTL_stats <- Fst_SNPwin_middleQTL %>% dplyr::group_by(Qid) %>%
                                                  dplyr::summarize(fst_win_mean = mean(Fst,na.rm = TRUE),fst_win_sd = sd(Fst, na.rm = TRUE),N_win = n()) %>%
                                                  dplyr::summarise(average_fst_win_mean = mean(fst_win_mean,na.rm = TRUE),
                                                                   average_fst_win_sd = mean(fst_win_sd, na.rm = TRUE), average_N_win = round(mean(N_win,na.rm = TRUE)))

win50SNP_BottomQTL_stats <-  Fst_SNPwin_bottomQTL %>% dplyr::group_by(Qid) %>%
                                                  dplyr::summarize(fst_win_mean = mean(Fst,na.rm = TRUE),fst_win_sd = sd(Fst, na.rm = TRUE),N_win = n()) %>%
                                                  dplyr::summarise(average_fst_win_mean = mean(fst_win_mean,na.rm = TRUE),
                                                                   average_fst_win_sd = mean(fst_win_sd, na.rm = TRUE), average_N_win = round(mean(N_win,na.rm = TRUE)))

# Combine the tibbles into one table
combined_SNPwin_stats <- dplyr::bind_rows(win50SNP_TopQTL_stats %>% mutate(group = "Fst_stats_SNPwin_top95Q"),
                                   win50SNP_MiddleQTL_stats %>% mutate(group = "Fst_stats_SNPwin_middile2575Q"),
                                   win50SNP_BottomQTL_stats %>% mutate(group = "Fst_stats_SNPwin_bottom5Q"))

# Print the combined table
print(combined_SNPwin_stats)

readr::write_csv(combined_SNPwin_stats, "Fst_stats_50SNPwin_QTLeffect_QTLgroups_W1Q1P2.csv")


