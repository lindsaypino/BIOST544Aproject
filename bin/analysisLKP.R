options(digits = 3) ## Formats output to 3 digits
library(gtools)
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

# note: missing prptype for pt AGFT199 -- discarded in excel


#
# DESCRIPTIVE STATISTICS
#

## Missingness by feature, ie how many N/As per column
apply(subset.df, 2, function(x){sum(is.na(x))}) 
complete <- subset.df %>% complete.cases()
mean(complete)  ## proportion of cases complete
length(unique(subset.df$patientid))  ## total unique patients

## subset for metadata features
meta.df <- subset.df %>% select(patientid, treatment, age, race, gender_num, duration_months)  # subset for features we'll be using 
apply(meta.df, 2, function(x){sum(is.na(x))}) 
complete <- meta.df %>% complete.cases()

## Summary statistics
ages <- subset.df %>% 
  group_by(as.factor(gender_num)) %>% 
  summarise(mean(duration_months))

hist(subset.df$age)
hist(subset.df$duration_months)

## Creating binary variable for improvement at 52 weeks- 
##will need write a function if no 52 week time piont available- to searh previos time 26, 12 or 6 for a value 

improved.last <- function(x) {
  if(x <= 4) {0} else {1}  } ##does R include NA in this data count?

improved.52 <- sapply(subset.df$time52lpt,improved.last)

##Summary mean change score scores first


all.df$pda0.transf <- 100 - all.df$time0pda
all.df$pda6.transf <- 100 - all.df$time6pda
all.df$pda12.transf <- 100 - all.df$time12pda
all.df$pda26.transf <- 100 - all.df$time26pda
all.df$pda52.transf <- 100 - all.df$time52pda

all.df$pse0.transf <- 100 - all.df$time0pse
all.df$pse6.transf <- 100 - all.df$time6pse
all.df$pse12.transf <- 100 - all.df$time12pse
all.df$pse26.transf <- 100 - all.df$time26pse
all.df$pse52.transf <- 100 - all.df$time52pse

##Standardize the outcomes scores


all.df %>%
  mutate(score = (x - mean(x)) / sd(x))

##Developing the propensity scores

all.ps <- glm(prponly_num ~ duration_m + gender_num + age + time0pda + time0pse + time0,
                 family = binomial(), data = all.df)
summary(subset.ps)

prs_df <- data.frame(pr_score = predict(all.ps, type = "response"),
                     prponly_num = all.df$model$prponly_num)
head(prs_df)

labs <- paste("Actual PRP Type:", c("LR-PRP", "LP-PRP"))
prs_df %>%
  mutate(prponly_num = ifelse(prponly == 1, labs[1], labs[2])) %>%
  ggplot(aes(x = pr_score)) +
  geom_histogram(color = "white") +
  facet_wrap(~prponly_num) +
  xlab("Probability of getting LR-PRP injection") +
  theme_bw()

##don't think we can really match by ps for this analysis bc of the imbalance in groups

## code prptype as a factor
subset.df$prptype <- as.character(subset.df$prptype)
subset.df$prpcode[subset.df$prptype == "cascade"] <- 1
subset.df$prpcode[subset.df$prptype == "biomet"] <- 2
subset.df$prpcode[subset.df$prptype == "abi"] <- 0
subset.df$prpcode <- as.factor(subset.df$prpcode)

## subset the data for the 'easy' variables of patientid, timepoints, and the coded prp type

timecols <- c('time0','time6','time12', 'time26', 'time52')
time.prptype.df <- subset.df[, c('patientid',timecols,'prpcode')]

## plot the trends over time for each patient (the other scale, not the likert scale)

meltedSubset <- melt(as.data.frame(time.prptype.df), id=c('patientid','prpcode'))

ggplot(meltedSubset, aes(x=variable, y=value, factor(prpcode))) +
  stat_summary(aes(group=patientid, color = as.factor(prpcode)), alpha=0.5, fun.y=mean, geom="line")

## plot subset of the patients trend over time (whatever subset you want)

plot_trend <- function(slice){
  ggplot(subset(meltedSubset, patientid %in% slice), 
         aes(x=variable,y=value)) +
    geom_point(aes(color=patientid), size=3) +
    stat_summary(aes(group=patientid, color=patientid), fun.y=mean, geom="line")
}

plot_trend(time.prptype.df$patientid[1:10])  # note: you can change this "[1:10]" to however many or few you want! 


## boxplot of trends
ggplot(meltedSubset, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill = factor(prpcode)))


## above with likert
likertcols <- c('time6lpt', 'time12lpt', 'time26lpt', 'time52lpt')
time.likert.df <- subset.df[, c('patientid',likertcols,'prpcode')]
meltedLikert <- melt(as.data.frame(time.likert.df), id=c('patientid','prpcode'))


ggplot(meltedLikert, aes(x=variable, y=value, factor(prpcode))) +
  stat_summary(aes(group=patientid, color = as.factor(prpcode)), alpha=0.5, fun.y=mean, geom="line")

## [LIKERT] plot subset of the patients trend over time (whatever subset you want)

plot_trend <- function(slice){
  ggplot(subset(meltedLikert, patientid %in% slice), 
         aes(x=variable,y=value)) +
    geom_point(aes(color=patientid), size=3) +
    stat_summary(aes(group=patientid, color=patientid), fun.y=mean, geom="line")
}

plot_trend(time.likert.df$patientid[1:10])  # note: you can change this "[1:10]" to however many or few you want! 


## [LIKERT] boxplot of trends
ggplot(meltedLikert, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill = factor(prpcode)))




#
# DATA CLEANING
#

# find last time point
find.last.time.point <- function(row, cols){
  # sort column names alphanumeric so that they're 'in order'
  cols <- mixedsort(cols)
  timepoints <- row[,cols]
  # if all NA then skip
  
  # else
  # go backwards until not NA
  
}
  
  










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
source("http://bioconductor.org/biocLite.R")
biocLite("edge")
library(edge)


## prepare matrix for modeling
matrixSubset <- as.matrix(as.data.frame(lapply(time.prptype.df[,timecols], as.numeric)))

## creating full and null models
library(splines)
timevariable <- c(0,6,12,26,52)
colnames(matrixSubset) <- timevariable
de_obj <- build_study(data = matrixSubset, tme=timevariable, sampling = "timecourse")
full_matrix <- fullMatrix(de_obj)
null_matrix <- nullMatrix(de_obj)
prpcode.expr <- exprs(de_obj)  # seems to be the same as the expression matrix


## fitting the data

ef_obj <- fit_models(de_obj, stat.type = "odp")
alt_fitted <- fitFull(ef_obj)
null_fitted <- fitNull(ef_obj)


## significance analysis


#de_odp <- odp(ef_obj, bs.its = 50, verbose = FALSE, n.mods = 50)  # optimal discovery procedure (odp) doesn't work
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
