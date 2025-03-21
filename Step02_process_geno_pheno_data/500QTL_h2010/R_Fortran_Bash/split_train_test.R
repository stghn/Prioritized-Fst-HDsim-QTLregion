#########################################################

#' Splitdf splits a dataframe into a training sample and test sample with a given proportion
#'
#' This function takes a data frame and according to predefined proportion "prop" it will return a training and a test sample
#'
#' @param input a n x p dataframe of n observations and p variables or a vector
#' @param seed the seed to be set in order to ensure reproductability of the split
#' @param prop the proportion of the training sample [0-1]
#' @return a list with two slots: trainset and testset
#' @author BlackGuru
#' @details
#' This function takes a data frame or a vector and according to predefined proportion "prop" it will return a training and a test sample. "prop" corresponds to the proportion of the training sample.
#' @export

splitdf <- function(input, prop=0.5, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.data.frame(input)){
    index <- 1:nrow(input)
    trainindex <- sample(index, trunc(length(index)*prop))
    trainset <- input[trainindex, ]
    testset <- input[-trainindex, ]
  }else if (is.vector(input)){
    index<-1:length(input)
    trainindex <- sample(index, trunc(length(index)*prop))
    trainset <- input[trainindex]
    testset <- input[-trainindex]
  }else{
    print("Input must be a dataframe or a vector")
  }
  list(trainset=trainset,testset=testset)
}


library(dplyr)
data <- read.table("output_rep1/p1_data_001.txt",head=T)
ID.gen10 <-data %>% dplyr::filter(G==10) %>% dplyr::select(Progeny)

splits <- splitdf(ID.gen10,prop=2/3, seed=808)

training <- splits$trainset
testing <- splits$testset

write.table(testing,"test.gen10.ID",col.names=F, row.names=F)


