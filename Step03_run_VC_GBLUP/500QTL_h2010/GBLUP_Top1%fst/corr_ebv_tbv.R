rm(list = ls())

## specify the package names and loading them using pacman package
library(pacman)
pacman::p_load(dplyr,tidyr,tidyverse,data.table)

###############################################################################################################
###### Reading multiple *.Rdata files from a directory and load all of them to environment workspace in R #####
###############################################################################################################
# Get the present working directory
pwd <- getwd()
# Add trailing slash to the directory path
pwd <- paste0(pwd, "/")      

file_list <- list.files(path=pwd, pattern="sol.tbv")           # create list of files

for (i in 1:length(file_list)){
  assign(file_list[i],as.data.frame(fread(paste(pwd, file_list[i], sep=''),head=F,stringsAsFactors = F, na.strings = c("."))))
}

rm(i);rm(file_list);rm(pwd)


######################################################################################################
##read dataframes contained "sol.tbv" from the environment in R and make a list file of those dataframes
GBLUP_pred<-grep("sol.tbv",names(.GlobalEnv),value=TRUE)
GBLUP_pred_list<-do.call("list",mget(GBLUP_pred))

lapply(GBLUP_pred_list, function(x) dim(x))  ### show the dimension for each dataframe in the list file

##calculate correlation between estimated and true breeding value predicted from GBLUP approach
corr_accu<-unlist(lapply(GBLUP_pred_list, function(x) cor(x[[3]],x[[4]],use="complete.obs",method = "pearson")))
corr_accu


