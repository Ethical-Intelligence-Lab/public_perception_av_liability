## ================================================================================================================
##                                 Harvard Business School, Ethical Intelligence Lab
## ================================================================================================================
##                                DATA ANALYSIS | AV SCENARIOS | WITHIN SUBJECTS               
## ================================================================================================================

## clear workspace
rm(list = ls()) 

options(download.file.method="libcurl")

## install packages
library(ggpubr)
library(rstatix)
if (!require(pacman)) {install.packages("pacman")}
pacman::p_load('ggplot2',         # plotting
               'ggsignif',        # plotting significance bars
               'lme4',            # functions for fitting linear regression models
               'ggforce',         # make ggplot even fancier
               'ggpubr',          # arrange plots in a grid, if needed
               'ltm',             # probably not using..
               'tidyr',           # tools for cleaning messy data
               'stringr',         # perform string substitutions easily
               'assertthat',      # allows me to check whether a variable is a string, with is.string
               'lsmeans',         # contrast analysis for regression models
               'stats',           # use function to adjust for multiple comparisons
               'filesstrings',    # create and move files
               'simr',            # power analysis for mixed models
               'compute.es',      # effect size package
               'effsize',         # another effect size package
               'pwr',             # package for power calculation
               'nlme',            # get p values for mixed effect model
               'DescTools',       # get Cramer's V
               'Hmisc',
               'effsize'          # effect size
)

library("lmerTest")

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
# if importing from Qualtrics: (i) export data as numeric values, and (ii) delete rows 2 and 3 of the .csv file.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory
d <- read.csv('e4_ws.csv') 

## perform attention exclusions: 
# this will remove responses from the dataframe that failed attention checks (i.e., "1" or "2")
d <- subset(d, (d$att1 == 2 & d$att2 == 2))
dim(d)

## ================================================================================================================
##                                              PERFORM EXCLUSIONS                
## ================================================================================================================

## get number of participants BEFORE exclusions: 
n_original <- dim(d)[1] # extracting number of rows only, not columns
n_original 

## perform comprehension exclusions: 
# this will remove responses from the dataframe that failed comprehension checks (i.e., "2")
d <- subset(d, ( d$Q132 == 1 & d$Q133 == 1 & d$comp1 == 3 & d$comp_accident == 1))
dim(d) # number of participants should decrease after comprehension exclusions

## get number of participants AFTER exclusions: 
n_final <- dim(d)[1] # extracting number of rows only, not columns
n_final 
percent_excluded <- (n_original - n_final)/n_original 
percent_excluded

## ================================================================================================================
##                                                    SUBSETTING                 
## ================================================================================================================

## define new data frame to extract pre-processed data into:
d_subset <- array(dim=c(dim(d)[1], 20))
colnames(d_subset) <- c('order', 
                        'vA_sue_AV', 'vB_sue_AV', 'defective_AV', 
                        'negligence_AV', 'counterfactual_AV', 'capability_AV', 'fault_AV', 
                        'vA_sue_HDV', 'vB_sue_HDV', 'defective_HDV', 
                        'negligence_HDV', 'counterfactual_HDV', 'capability_HDV', 'fault_HDV',
                        'superhuman', 'ai_familiarity', 'insurance_fmlrty',
                        'mot_insurance_fmlrty', 'prod_liab_fmlrty')
d_subset <- as.data.frame(d_subset, stringsAsFactors=FALSE) 

## extract good data from the middle part of raw data:
for(i in 1:dim(d)[1]) {
    curr_AV <- d[i,c(29:35)][!is.na(d[i,c(29:35)])] # for a given row, get only the non-NA values
    curr_HDV <- d[i,c(46:52)][!is.na(d[i,c(46:52)])]
    curr_other <- d[i,c(55:59)][!is.na(d[i,c(55:59)])]
    d_subset[i,2:8] <- as.numeric(curr_AV[curr_AV!= ""])# and only the non-empty values
    d_subset[i,9:15] <- as.numeric(curr_HDV[curr_HDV!= ""])
    d_subset[i,16:20] <- as.numeric(curr_other[curr_other!= ""])
    d_subset[i,1] <- d[i,'order']
}

d_merged <- array(dim=c(dim(d)[1]*2, 10))
colnames(d_merged) <- c('order', 'cond',
                        'vA_sue', 'vB_sue', 'defective', 
                        'negligence', 'counterfactual', 'capability', 'fault')
d_merged <- as.data.frame(d_merged, stringsAsFactors=FALSE) 
d_merged$order <- c(d_subset$order, d_subset$order)
d_merged$cond <- c(rep(c(1),each=dim(d)[1]), rep(c(2),each=dim(d)[1])) #1= AI, 2 =Human

d_merged$vA_sue <- c(d_subset$vA_sue_AV, d_subset$vA_sue_HDV)
d_merged$vB_sue <- c(d_subset$vB_sue_AV, d_subset$vB_sue_HDV)
d_merged$defective <- c(d_subset$defective_AV, d_subset$defective_HDV)
d_merged$negligence <- c(d_subset$negligence_AV, d_subset$negligence_HDV)
d_merged$counterfactual <- c(d_subset$counterfactual_AV, d_subset$counterfactual_HDV)
d_merged$capability <- c(d_subset$capability_AV, d_subset$capability_HDV)
d_merged$fault <- c(d_subset$fault_AV, d_subset$fault_HDV)

# creates a new variable for the string version of condition
d_merged$cond_name <- ifelse(d_merged$cond==1, "av", "human")

## ================================================================================================================
##                                            PARTICIPANT CHARACTERISTICS                 
## ================================================================================================================

# ## % prior experience apps
# table(d$ai_familiarity)[1]/sum(table(d$ai_familiarity))
# d$ai_familiarity
# #more like 11
# 11/dim(d)[1]

## age
mean(as.numeric(d$age), trim = 0, na.rm = TRUE) ## mean age 

## gender
table(d$gender)[2]/sum(table(d$gender)) ## percentage of females

## ai familiarity
mean(as.numeric(d$ai_familiarity), trim = 0, na.rm = TRUE) ## mean ai familiarity 

## insurance familiarity
mean(as.numeric(d$insurance_fmlrty), trim = 0, na.rm = TRUE) ## mean insurance familiarity 

## motor insurance familiarity
mean(as.numeric(d$mot_insurance_fmlrty), trim = 0, na.rm = TRUE) ## mean motor insurance familiarity

## motor insurance familiarity
mean(as.numeric(d$prod_liab_fmlrty), trim = 0, na.rm = TRUE) ## mean product liability familiarity


## ================================================================================================================
##                              DATA ANALYSIS - SUMMARY, T-TESTS, AND LINEAR REGRESSION               
## ================================================================================================================

## (1) VA SUE
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(vA_sue, type = "mean_sd")
mean(as.numeric(d_subset$vA_sue_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vA_sue_AV), na.rm = TRUE)
mean(as.numeric(d_subset$vA_sue_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vA_sue_HDV), na.rm = TRUE)

vA_sue_mod <- lmer(vA_sue ~ cond + (1 | order), data=d_merged)
summary(vA_sue_mod)

## (2) VB SUE
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(vB_sue, type = "mean_sd")
mean(as.numeric(d_subset$vB_sue_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vB_sue_AV), na.rm = TRUE)
mean(as.numeric(d_subset$vB_sue_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vB_sue_HDV), na.rm = TRUE)

vB_sue_mod <- lmer(vB_sue ~ cond + (1 | order), data=d_merged)
summary(vB_sue_mod)

## (3) DEFECTIVE
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(defective, type = "mean_sd")
mean(as.numeric(d_subset$defective_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$defective_AV), na.rm = TRUE)
mean(as.numeric(d_subset$defective_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$defective_HDV), na.rm = TRUE)

defective_mod <- lmer(defective ~ cond + (1 | order), data=d_merged)
summary(defective_mod)

## (4) COUNTERFACTUAL
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(counterfactual, type = "mean_sd")
mean(as.numeric(d_subset$counterfactual_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$counterfactual_AV), na.rm = TRUE)
mean(as.numeric(d_subset$counterfactual_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$counterfactual_HDV), na.rm = TRUE)

counterfactual_mod <- lmer(counterfactual ~ cond + (1 | order), data=d_merged)
summary(counterfactual_mod)

## Other variables:

## NEGLIGENCE
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(negligence, type = "mean_sd")

negligence_mod <- lmer(negligence ~ cond + (1 | order), data=d_merged)
summary(negligence_mod)

## CAPABILITY
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(capability, type = "mean_sd")

capability_mod <- lmer(capability ~ cond + (1 | order), data=d_merged)
summary(capability_mod)

## FAULT
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(fault, type = "mean_sd")

fault_mod <- lmer(fault ~ cond + (1 | order), data=d_merged)
summary(fault_mod)

## SUPERHUMAN
## get summary statistics
mean(as.numeric(d_subset$superhuman), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$superhuman), na.rm = TRUE)

## ================================================================================================================
##                                              MEDIATION ANALYSIS                
## ================================================================================================================

source("../process.R")

fit <- process(data = d_merged, y = "vB_sue", x = "cond", 
        m =c("counterfactual", "defective"), model = 6, effsize =1, total =1, stand =1, 
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)

## ================================================================================================================
##                                              PLOTTING MAIN FIGURES                 
## ================================================================================================================

## plotting all measures
t_names <- c("AV", "HDV")
title_size <- 12
axis_size <- 16

## (1) VA SUE
p1 <- ggplot(d_merged,aes(x=factor(cond_name),y=vA_sue)) +  
  theme_bw() +coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+ 
  geom_signif(comparisons = list(c("av", "human")), annotation="*", textsize = 5.5)

p1 <- p1 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Sue Veh. A Driver") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=title_size, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1 


## (2) VB SUE
p2 <- ggplot(d_merged,aes(x=factor(cond_name),y=vB_sue)) +  
  theme_bw() +coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="***", textsize = 5.5)

p2 <- p2 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Sue Veh. B Manufacturer") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=title_size, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p2 


## (3) DEFECTIVE
p3 <- ggplot(d_merged,aes(x=factor(cond_name),y=defective)) +  
  theme_bw() +coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="***", textsize = 5.5)

p3 <- p3 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Defective") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=title_size, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75)+  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p3

## (4) NEGLIGENCE
p4 <- ggplot(d_merged,aes(x=factor(cond_name),y=negligence)) +  
  theme_bw() +coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="*", textsize = 5.5)

p4 <- p4 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Negligence") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=title_size, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75)+  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p4

## (5) COUNTERFACTUAL
p5 <- ggplot(d_merged,aes(x=factor(cond_name),y=counterfactual)) +  
  theme_bw() +coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="***", textsize = 5.5)

p5 <- p5 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Counterfactual") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=title_size, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75)+  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p5

## (6) CAPABILITY
p6 <- ggplot(d_merged,aes(x=factor(cond_name),y=capability)) +  
  theme_bw() +coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="***", textsize = 5.5)

p6 <- p6 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Capability to Avoid") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=title_size, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75)+  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p6

## (7) FAULT
p7 <- ggplot(d_merged,aes(x=factor(cond_name),y=fault)) +  
  theme_bw() +coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="NS", textsize = 3)

p7 <- p7 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Avoid when not at fault") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=title_size, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75)+  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p7

## (4) ALL FIGURES
dev.new(width=12,height=4,noRStudioGD = TRUE)

figure <- ggarrange(p1, p2, p5, p3, nrow=1,ncol=4,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure <- annotate_figure(figure,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                bottom = text_grob("Vehicle Type", color="black", face ="plain",size=16)) 
plot(figure)

## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================