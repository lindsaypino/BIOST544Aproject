options(digits = 3) ## Formats output to 3 digits
library(ggplot2)
library(dplyr)
library(data.table)



#
# CREATION OF DATAFRAME AND PREPROCESSING
#

dir <- getwd()
setwd(dir)
all.df <- read.csv("./data/prp_tendon_cleaned2_reshaped_delim.csv")

## Subset for the variables I know
cols <- colnames(all.df)
keepcols <- c(1,2,3,4,5,6,28,29,30,31,33,34,36,37,38,40,41,42,43,44,45,46,49,50,52,53,55,56,57,58)
subset.df <- all.df[, keepcols]

## Data cleaning
# TODO: prptype -- what is biomet, what is cascade?
# missing prptype for pt AGFT199 -- discarded


## Missingness by feature, ie how many N/As per column
apply(subset.df, 2, function(x){sum(is.na(x))}) 
complete <- subset.df %>% complete.cases()
mean(complete)  ## proportion of cases complete
length(unique(subset.df$patientid))  ## total unique patients

## subset for metadata features
meta.df <- subset.df %>% select(patientid, treatment, age, race, gender, duration_months)  # subset for features we'll be using 
apply(meta.df, 2, function(x){sum(is.na(x))}) 
complete <- meta.df %>% complete.cases()


## Summary statistics
ages <- subset.df %>% 
    group_by(as.factor(gender)) %>% 
    summarise(mean(duration_months))

hist(subset.df$age)
hist(subset.df$duration_months)





##
## LASSO REGRESSION
##


















##
## TIMECOURSE ANALYSIS
##
## For now I'm not working much on this, 
## I'll have to think about it a little more. Namely, the package was 
## designed to be used on data that has many genes, and finding the genes
## that trend significantly over time. We don't have genes here, just
## the two groups. Maybe we don't need to do the significance analysis then?

library(reshape2)

time.prptype.df <- subset.df[, c(2,3,4,5,6,14)]

matrixSubset <- as.matrix(as.data.frame(lapply(time.df, as.numeric)))
colnames(matrixRank) <- sub.samp


source("http://bioconductor.org/biocLite.R")
biocLite("edge")
library(edge)


sub.samp.unique <- as.numeric(c(23.1, 12.3, 5.5, 0.7, 0.6, 0.2, 23.11, 12.31, 5.51, 0.71, 0.61, 0.21, 23.09, 12.29, 5.49, 0.69, 0.59))

colnames(matrixRank) <- sub.samp.unique
sub.samp.sort <- as.character(sort(sub.samp.unique))
matrixRank2 <- matrixRank[,sub.samp.sort]  

meltedRank <- melt(as.data.frame(matrixRank2))
meltedRank$peptide <- 1:nrow(matrixRank2)
meltedRank$variable <- as.numeric(meltedRank$variable)

plot_trend <- function(slice){
  ggplot(subset(meltedRank, peptide %in% slice), 
         aes(x=variable,y=value)) +
    geom_point(aes(color=peptide), size=3) +
    stat_summary(aes(group=peptide, color=peptide), fun.y=mean, geom="line")
}

plot_trend(c(1,2,3,4,5,6,7,8,9,10))
plot_trend(c(11,12,13,14,15,16,17,18,19,20))
plot_trend(c(9000,9001,9002,9003,9004,9005,9006,9007,9008,9009))


## creating full and null models

library(splines)
de_obj <- build_study(data = matrixRank2, tme = sub.samp.unique, sampling = "timecourse")
full_matrix <- fullMatrix(de_obj)
null_matrix <- nullMatrix(de_obj)
cr.expr <- exprs(de_obj)  # seems to be the same as the expression matrix


## fitting the data

ef_obj <- fit_models(de_obj, stat.type = "odp")
alt_fitted <- fitFull(ef_obj)
null_fitted <- fitNull(ef_obj)


## significance analysis


de_odp <- odp(ef_obj, bs.its = 50, verbose = FALSE, n.mods = 50)  # optimal discovery procedure (odp) doesn't work
de_lrt <- lrt(de_obj, nullDistn = "normal")  # likelihood ratio test does work
summary(de_lrt)

sig_results <- qvalueObj(de_lrt)
qvals <- sig_results$qvalues
pvals <- sig_results$pvalues
lfdr <- sig_results$lfdr
pi0 <- sig_results$pi0


# list of all significant genes at an FDr of 0.1
fdr.level <- 0.1
sigGenes <- qvals < fdr.level


##
##
##
