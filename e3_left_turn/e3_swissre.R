## clear workspace
rm(list = ls()) 

source("../process.R")

options(download.file.method="libcurl")

## install packages
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

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
# if importing from Qualtrics: (i) export data as numeric values, and (ii) delete rows 2 and 3 of the .csv file.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory
d <- read.csv('e3_data.csv')

## explore dataframe: 
dim(d) # will provide dimensions of the dataframe by row [1] and column [2]
colnames(d) # will provide all column names
summary(d)

## perform attention exclusions: 
# this will remove responses from the dataframe that failed attention checks (i.e., "1" or "2")
d <- subset(d, (d$att1 == 2 & d$att2 == 2))
dim(d) # number of participants should decrease after attention exclusions

## split up dataframes between AV and HDV conditions
## this is necessary before comprehension exclusions
d_AV <- subset(d, (d$FL_4_DO == "FL_39"))
d_HDV <- subset(d, (d$FL_4_DO == "FL_40"))

## get number of participants BEFORE exclusions: 
n_original <- dim(d)[1] # extracting number of rows only, not columns
n_original_AV <- dim(d_AV)[1]
n_original_HDV <- dim(d_HDV)[1]

## perform comprehension exclusions separately for AV and HDV: 
# this will remove responses from the dataframe that failed comprehension checks (i.e., "2")
d_AV <- subset(d_AV, (d_AV$comp1 == 1 & d_AV$comp_accident == 1))
d_HDV <- subset(d_HDV, (d_HDV$comp1 == "" & d_HDV$comp_accident == 1))
dim(d_AV) # number of participants should decrease after comprehension exclusions
dim(d_HDV)

## get number of participants AFTER exclusions: 
n_final_AV <- dim(d_AV)[1] # extracting number of rows only, not columns
n_final_HDV <- dim(d_HDV)[1]
n_final <- n_final_AV + n_final_HDV
percent_excluded <- (n_original - n_final)/n_original
percent_excluded_AV <- (n_original_AV - n_final_AV)/n_original_AV 
percent_excluded_HDV <- (n_original_HDV - n_final_HDV)/n_original_HDV

## remove unused columns according to condition
d_AV <- d_AV[-c(21:28,36:50)]
d_HDV <- d_HDV[-c(21:35,36:43)]

## get mean age and gender:
mean_age = mean(as.numeric(d$age), na.rm = TRUE); mean_age # removing NAs from the dataframe before computing mean 
gender = table(d$gender)[2]/sum(table(d$gender)); gender

## ================================================================================================================
##                                                    SUBSETTING                 
## ================================================================================================================

## define new data frame to extract pre-processed data into:
d_subset <- array(dim=c(n_final, 13))
colnames(d_subset) <- c('cond', 'vA_sue', 'vB_sue', 'defective', 'negligence', 'counterfactual', 
                        'capability', 'fault', 'superhuman', 'comp1', 'comp2', 'familiarity', 'mod')
d_subset <- as.data.frame(d_subset, stringsAsFactors=FALSE) 

## assess moderator data from both AV and HDV
moderator_mat = rbind(d_AV[31:45], d_HDV[31:45])
moderator_mat <- data.frame(sapply(moderator_mat, as.numeric))
cb_alpha = cronbach.alpha(moderator_mat)

#d_AV$moderator <- rowMeans(sapply(d_AV[31:35], as.numeric))
#d_HDV$moderator <- rowMeans(d_HDV[31:35])

moderator_mat$moderator <- rowMeans(moderator_mat)

## extract good data from the middle part of raw data in AV:
for(i in 1:n_final_AV) {
    curr <- d_AV[i,21:31][!is.na(d_AV[i,21:31])] # for a given row, get only the non-NA values
    d_subset[i,2:12] <- as.numeric(curr[curr!= ""]) # and only the non-empty values
    d_subset[i,13] <- moderator_mat$moderator[i]
    d_subset[i,1] <- d_AV[i,53][!is.na(d_AV[i,53])]
}

## extract good data from the middle part of raw data in HDV:
for(i in 1:n_final_HDV) {
    j = i+n_final_AV
    curr <- d_HDV[i,21:31][!is.na(d_HDV[i,21:31])] # for a given row, get only the non-NA values
    d_subset[j,2:12] <- as.numeric(curr) # and only the non-empty values
    d_subset[j,13] <- moderator_mat$moderator[j]
    d_subset[j,1] <- d_HDV[i,53][!is.na(d_HDV[i,53])]
}

## just to keep the df names straight for next section
d_merged <- d_subset

names(d_merged)[names(d_merged) == 'defective'] <- 'defec'
names(d_merged)[names(d_merged) == 'negligence'] <- 'negl'
names(d_merged)[names(d_merged) == 'counterfactual'] <- 'countf'
names(d_merged)[names(d_merged) == 'capability'] <- 'capab'
names(d_merged)[names(d_merged) == 'superhuman'] <- 'superh'

d_merged$countf2 <- d_merged$countf

d_merged$cond_name <- ifelse(d_merged$cond=="FL_39", "av", "human")
d_merged$cond_n <- ifelse(d_merged$cond_name=="av", 1, 2)

## ================================================================================================================
##                                             DATA ANALYSIS - T-TESTS               
## ================================================================================================================

table(d_merged$con) #give us table of number of people in each condition - want to have equal number of people in each condition

## (1) SUE VEHICLE A DRIVER
vA_sue_T <- t.test(vA_sue ~ cond_name, data = d_merged, paired = FALSE); vA_sue_T
cohen.d(vA_sue ~ cond_name, data = d_merged, paired = FALSE)
sd(d_merged[d_merged$cond_name == "av",]$vA_sue)
sd(d_merged[d_merged$cond_name == "human",]$vA_sue)

## (2) SUE VEHICLE B MANUFACTURER
vB_sue_T <- t.test(vB_sue ~ cond_name, data = d_merged, paired = FALSE); vB_sue_T
cohen.d(vB_sue ~ cond_name, data = d_merged, paired = FALSE)
sd(d_merged[d_merged$cond_name == "av",]$vB_sue)
sd(d_merged[d_merged$cond_name == "human",]$vB_sue)

#one-sided test
vB_sue_T <- t.test(vB_sue ~ cond_name, data = d_merged, paired = FALSE, alternative='greater'); vB_sue_T
cohen.d(vB_sue ~ cond_name, data = d_merged, paired = FALSE, alternative='greater')

## (3) VEHICLE B DEFECTIVE
defective_T <- t.test(defec ~ cond_name, data = d_merged, paired = FALSE); defective_T
cohen.d(defec ~ cond_name, data = d_merged, paired = FALSE)
sd(d_merged[d_merged$cond_name == "av",]$defec)
sd(d_merged[d_merged$cond_name == "human",]$defec)

## (4) COUNTERFACTUAL
counterfactual_T <- t.test(countf ~ cond_name, data = d_merged, paired = FALSE); counterfactual_T 
cohen.d(countf ~ cond_name, data = d_merged, paired = FALSE)
sd(d_merged[d_merged$cond_name == "av",]$countf)
sd(d_merged[d_merged$cond_name == "human",]$countf)


## Other variables:

## VEHICLE B NEGLIGENT
negligent_T <- t.test(negl ~ cond_name, data = d_merged, paired = FALSE); negligent_T 
sd(d_merged[d_merged$cond_name == "av",]$negl)
sd(d_merged[d_merged$cond_name == "human",]$negl)

## CAPABILITY
capability_T <- t.test(capab ~ cond_name, data = d_merged, paired = FALSE)
sd(d_merged[d_merged$cond_name == "av",]$capab)
sd(d_merged[d_merged$cond_name == "human",]$capab)

## FAULT
fault_T <- t.test(fault ~ cond_name, data = d_merged, paired = FALSE) 
sd(d_merged[d_merged$cond_name == "av",]$fault)
sd(d_merged[d_merged$cond_name == "human",]$fault)

## SUPERHUMAN
superhuman_T <- t.test(superh ~ cond_name, data = d_merged, paired = FALSE) 
sd(d_merged[d_merged$cond_name == "av",]$superh)
sd(d_merged[d_merged$cond_name == "human",]$superh)

cor(d_merged[,2:9])

## ================================================================================================================
##                                             MEDIATION ANALYSIS              
## ================================================================================================================

d_merged$cond_n <- ifelse(d_merged$cond=="FL_39", 1, 2)

# MODERATED SERIAL MEDIATION
# 87 = B path, 83 = A path, 91 = center path
process(data = d_merged, y = "vB_sue", x = "cond_n", 
        m =c("countf", "defec"), w = "mod", model = 87, effsize =1, total =1, stand =1, 
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)

# SERIAL MEDIATION
process(data = d_merged, y = "vB_sue", x = "cond_n", 
        m =c("countf", "defec"), model = 6, effsize =1, total =1, stand =1, 
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)

## ================================================================================================================
##                                              PLOTTING MAIN FIGURES                 
## ================================================================================================================

## plotting all measures
## FL39 --> AV condition; FL40 --> HDV condition
t_names <- c("AV", "HDV")
title_size <- 20

## (1) Sue VA driver
p1_1 <- ggplot(d_merged,aes(x=factor(cond_name),y=vA_sue)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="**", textsize = 5.5)

p1_1 <- p1_1 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Sue Veh. A Driver") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_1

## (2) Sue VB manufacturer
p1_2 <- ggplot(d_merged,aes(x=factor(cond_name),y=vB_sue)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="NS", textsize = 3)

p1_2 <- p1_2 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Sue Veh. B Manufacturer") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_2

## (3) VB Defective
p1_3 <- ggplot(d_merged,aes(x=factor(cond_name),y=defec)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="***", textsize = 5.5)

p1_3 <- p1_3 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Defective") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_3

## (4) VB Negligence
p1_4 <- ggplot(d_merged,aes(x=factor(cond_name),y=negl)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="***", textsize = 5.5)

p1_4 <- p1_4 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Negligence") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_4

## (5) Counterfactual
p1_5 <- ggplot(d_merged,aes(x=factor(cond_name),y=countf)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="***", textsize = 5.5)

p1_5 <- p1_5 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Counterfactual") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_5

## (6) Capability to Avoid
p1_6 <- ggplot(d_merged,aes(x=factor(cond_name),y=capab)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="***", textsize = 5.5)

p1_6 <- p1_6 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Capability to Avoid") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_6

## (7) Avoid when not at fault
p1_7 <- ggplot(d_merged,aes(x=factor(cond_name),y=fault)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="*", textsize = 5.5)

p1_7 <- p1_7 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Avoid when not at fault") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_7

## (8) Superhuman
p1_8 <- ggplot(d_merged,aes(x=factor(cond_name),y=superh)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation="*", textsize = 5.5)

p1_8 <- p1_8 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Superhuman") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_8

## PLOT SERIES 1
#dev.new(width=13,height=6,noRStudioGD = TRUE)
#figure1 <- ggarrange(p1_1, p1_2, p1_3, p1_4, p1_5, p1_6, p1_7, p1_8, nrow=2,ncol=4,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
#annotate_figure(figure1,left = text_grob("Mean Rating", color="black", face ="plain",size=16, rot=90),
#                bottom = text_grob("Scenario Condition", color="black", face ="plain",size=16)) 

dev.new(width=12,height=4,noRStudioGD = TRUE)
figure1 <- ggarrange(p1_1, p1_2, p1_5, p1_3, nrow=1,ncol=4,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure1 <- annotate_figure(figure1,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                bottom = text_grob("Vehicle Type", color="black", face ="plain",size=16))

plot(figure1)

#export data for analysis in E5
write.csv(d_merged, 'e3_processed.csv')

## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================