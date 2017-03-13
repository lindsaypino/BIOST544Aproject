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
all.df <- read.csv("../data/prp_tendon_cleaned2_reshaped_delim_20170310.csv")

## Subset for the variables I know
cols <- colnames(all.df)
keepcols <- c(1,2,3,4,5,6,7,8,9,10,11,12,18,19,20,21,22,28,29,30,31,32,33,34,36,37,38,40,41,42,43,44,45,46,49,50,52,53,55,56,57,58,59,71,72,73,74,75,80,81,82,83,84)
submz.df <- all.df[, keepcols]
subset(submz.df, treatment!="ABI")


# note: missing prptype for pt AGFT199 -- discarded in excel

##make prponly_num a binary variable
submz.df$prponly_num <- submz.df$prponly_num-1 

##Mean Difference before confounding

(improvement.by.prptype <- submz.df %>% 
    group_by(as.factor(prptype)) %>% 
    summarise(prop.inf = mean(time52lpt)))

(prop.diff <- improvement.by.prptype$prop.inf[1] - improvement.by.prptype$prop.inf[0])

do.one <- function(outcome, label){
  perm.label <- sample(label)
  return(mean(outcome[perm.label == 1]) - mean(outcome[perm.label == 0]))
}

sampling.dist <- with(submz.df,
                      replicate(1e4, do.one(pct.likert52, prptype)))

ggplot(data.frame(perm.prop.diff = sampling.dist), aes(x = perm.prop.diff, y=..density..)) +
  geom_density() + 
  geom_vline(xintercept = prop.diff, color = "red")

mean(sampling.dist > prop.diff)

## Regression before PS Score weighting
prptrtraw <- lm(pct.likert52 ~ prponly_num + pct.base52 + age + gender + duration_months + dxcode + time52pda + time52pse, data = submz.df)
                summary(prptrtraw)
                
time0pda.ATE <- lm(time0pda ~ prponly_num, data=submz.df)
time0pse.ATE <-lm(time0pse ~ prponly_num, data=submz.df)

summary(time0pda.ATE)
summary(time0pse.ATE)

##Calculate the propensity score

apply(submz.df, 2, function(x){sum(is.na(x))}) 
complete <- submz.df %>% complete.cases()
mean(complete)  ## proportion of cases complete
length(unique(submz.df$patientid))  ## total unique patients


all.ps <- glm(prponly_num ~ duration_months + gender + age + time0pda + time0pse + time0 + dxcode,
              data = submz.df, family = binomial())
summary(all.ps)

##Attach PS to the data
psvalue <- predict(all.ps, data=submz.df, type = "response")

ggplot(data.frame(psvalue = psvalue), aes(x = psvalue, y = ..density..)) + geom_histogram()

range(psvalue)

trunc.propen <- psvalue %>% pmin(0.95) %>% pmax(0.05)


##Look at ratio of enrollment probability
npat <- nrow(submz.df)
weights <- rep(0, npat)

## for patients who were assigned to Cascade:
representative.propen <- sum(submz.df$prponly_num) / npat
actual.propen <- trunc.propen

prptype.ind <- which(submz.df$prponly_num == 1)
weights[prptype.ind] <- representative.propen/actual.propen[prptype.ind]
weights[-prptype.ind]<- (1 - representative.propen)/(1- actual.propen[-prptype.ind])

ggplot(data.frame(weights = weights), aes(x = weights, y = ..density..)) + geom_histogram()


(cascade.prob.est <- with(submz.df,
                          mean((weights*pct.likert52)[prptype.ind])))

(biomet.prob.est <- with(submz.df,
                          mean((weights*pct.likert52)[-prptype.ind])))

(diff.est <- cascade.prob.est - biomet.prob.est)

##Inverse weight and use propensity estimates in determining treatment labels

do.one.propen <- function(outcome, propen){
  n <- length(outcome)
  label <- rbinom(n,1,propen)
  
  weights <- rep(0,n)  
  representative <- mean(label)
  actual <- propen
  ind.t <- which(label == 1)
  weights[ind.t] <- (representative/actual)[ind.t]
  weights[-ind.t] <- ((1-representative)/(1-actual))[-ind.t]
  
  return(mean((weights*outcome)[ind.t]) - mean((weights*outcome)[-ind.t]))
}

rerandomized.diffs <- replicate(1e3, do.one.propen(submz.df$pct.likert52, trunc.propen))

ggplot(data.frame(diffs = rerandomized.diffs), aes(x = diffs, y = ..density..)) +
  geom_density() + 
  geom_vline(xintercept = diff.est, color = "red")

mean(rerandomized.diffs > diff.est)


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
  if(x <= 4) {0} else {1}  } ##does R include NA in this data count? How do you tell R to lioke at one precvious?

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

##prp_dich <- funcion(all.df$prponly_num){
  ##if(all.df$prponly_num==2) (prp_dich == 0)
  ##if (all.df$prponly_num==1) (prp_dich == 1)
  ##else NA
##}
  
#all.df$prponly_num <- is.numeric(all.df$prponly_num)
#all.df$prpcode[all.df$prponly_num == is.NA] <- is.NA
#all.df$prpcode[all.df$prponly_num == 1] <- 1
#all.df$prpcode[all.df$prponly_num == 2] <- 0
#all.df$prpcode <- as.factor(all.df$prpcode) ##all listed as NA right now!



all.df %>%
  mutate(score = (x - mean(x)) / sd(x))

##Developing the propensity scores

all.ps <- glm(prponly_num ~ duration_m + gender_num + age + time0pda + time0pse + time0,
                 family = binomial(), data = all.df)
summary(all.ps)

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
