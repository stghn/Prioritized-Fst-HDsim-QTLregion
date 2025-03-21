###########################################################
## Additive genetic variance of QTLs (Aqtl=2pq(alpha)^2)###
###########################################################
rm(list=ls())

library(pacman)
pacman::p_load(dplyr)

# Read the data from the file "QTLinfo" into a data frame and consider the first row as headers
QTLinfo <- read.table("QTLinfo", header = TRUE)

# Set the option to prevent scientific notation
options(scipen = 999)

# Calculate the minimum of Freq1 and Freq2 and store it in column 'q'
QTLinfo <- QTLinfo %>% dplyr::mutate(q = pmin(Freq1, Freq2))

# Calculate the maximum of Freq1 and Freq2 and store it in column 'p'
QTLinfo <- QTLinfo %>% dplyr::mutate(p = pmax(Freq1, Freq2))

# Calculate the additive genetic variance of QTLs and store it in 'Aqtl'
QTLinfo <- QTLinfo %>% dplyr::mutate(Aqtl = 2 * p * q * (Alpha)^2)

# Set the constant percentage of Additive genetic variance of QTLs
Vg <- 0.10

# Calculate the percentage of Additive genetic variance of QTLs and store it in 'Aqtl.percent'
QTLinfo <- QTLinfo %>% dplyr::mutate(Aqtl.percent = (Aqtl / Vg) * 100)

# Exclude columns Freq1, Freq2, and Aqtl from the selection and store the result in QTL1
QTL1 <- QTLinfo %>% 
          dplyr::select(-Freq1, -Freq2,-Aqtl)

write.table(QTL1,"QTL_genetic_percentage",col.names=T, row.names=F,quote=F)

