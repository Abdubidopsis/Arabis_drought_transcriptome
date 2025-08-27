#### Dry_down Phenotype data analyses, updated Feb 2025  --- Abdul Saboor Khan, PhD University of Cologne ####
while (!is.null(dev.list()))  dev.off()
dev.off()



# Load BiocManager or install if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(BiocManager)

# Complete list of packages
packages <- c(
  # CRAN packages
  "reshape2", "stringi", "stringr", "ggplot2", "readr", "dplyr", 
  "lme4", "lmerTest", "tidyverse", "tidyr", "multcomp", "pegas", 
  "hierfstat", "adegenet", "vcfR", "ape", "heatmap3", "broom", 
  "ggmap", "ggsn", "ggpubr", "agricolae", "GGally", "mlbench", 
  "caret", "RColorBrewer", "corrplot", "PerformanceAnalytics", 
  "qgraph", "grid", "data.table", "gridExtra", "StAMPP",
  
  # Additional packages from recent lists
  "pasilla", "ashr", "apeglm", "DESeq2", "BiocParallel", "IHW", 
  "iSEE", "Glimma", "vsn", "pheatmap", "EnhancedVolcano", "edgeR",
  "poppr", "SNPRelate"
)

# Separate CRAN and Bioconductor packages
cran_packages <- c(
  "reshape2", "stringi", "stringr", "ggplot2", "readr", "dplyr", 
  "lme4", "lmerTest", "tidyverse", "tidyr", "multcomp", "pegas", 
  "hierfstat", "adegenet", "vcfR", "ape", "heatmap3", "broom", 
  "ggmap", "ggsn", "ggpubr", "agricolae", "GGally", "mlbench", 
  "caret", "RColorBrewer", "corrplot", "PerformanceAnalytics", 
  "qgraph", "grid", "data.table", "gridExtra", "StAMPP", "poppr"
)

bioconductor_packages <- c(
  "pasilla", "ashr", "apeglm", "DESeq2", "BiocParallel", "IHW", 
  "iSEE", "Glimma", "vsn", "pheatmap", "EnhancedVolcano", "edgeR", 
  "SNPRelate"
)

# Install missing CRAN packages
cran_to_install <- setdiff(cran_packages, rownames(installed.packages()))
if (length(cran_to_install) > 0) install.packages(cran_to_install, dependencies = TRUE)

# Install missing Bioconductor packages
bioconductor_to_install <- setdiff(bioconductor_packages, rownames(installed.packages()))
if (length(bioconductor_to_install) > 0) BiocManager::install(bioconductor_to_install, dependencies = TRUE)

# Load all packages
lapply(packages, library, character.only = TRUE)













library(reshape2)
library(stringi)
library(stringr)
library(ggplot2)
library(readr)
library(reshape2)
library(stringi)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(tidyverse)
library(dplyr)
library(tidyr)
library(multcomp)
library(pegas)
library(hierfstat)
library(adegenet)
library(vcfR)
library(ape)
library(heatmap3)
library(broom)
library(ggmap)
library(ggsn)
library(ggpubr)
library(agricolae)
library(GGally)
library(dplyr)
library(mlbench)
library(caret)
library(RColorBrewer)
library(corrplot)
library(PerformanceAnalytics)
library(qgraph)
#library(PopGenome)
##########################################################################
##################### Import data file ###################################
##########################################################################
setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/")
data<-read.csv("Dataphenotype.csv", header=TRUE, sep=",", na = "NA")
summary(data)
##############################################################################################
################### statistical analysis (models and plots) ##################################
##############################################################################################
#########################################################################
##### melt the days data during dry down and plot #######################
#########################################################################
pheno=melt(data,id.vars = c("Pot_no","tray","genotype","surv","Smw","LT","LTatwilt","DaystoWilt", "DoD","RA","DaystoRecov", "LfW","LdW","LA","LL","SD","SS","wounding"),measure.vars = c("Day_00","Day_01","Day_02","Day_03","Day_04","Day_05","Day_06","Day_07","Day_08","Day_09","Day_10","Day_11","Day_12","Day_13","Day_14","Day_15"))
head(pheno)
glimpse(pheno)
summary(pheno)
## Make separate column for "days" variable
pheno$days=as.numeric(stri_sub(pheno$variable,-2, -1))
#mositure content of the soil from start until wilting
#ggplot(data=pheno,aes(x=days,y=value,group=Pot_no,col=genotype))+geom_line()+xlab("Days")+ylab("Value")
pheno_m_mean=group_by(pheno,genotype,days) %>% summarise(mean_bz=mean(value,na.rm=T),sd_bz=sd(value,na.rm=T))
a <- ggplot(data=pheno_m_mean,aes(x=days,y=mean_bz, color=genotype, fill=genotype))+geom_line()+ geom_ribbon(aes(ymin=mean_bz-sd_bz,ymax=mean_bz+sd_bz),alpha=0.2)+xlab("Days")+ylab("Weight(g)") +ylab("Soil water content")+ theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24), axis.title = element_text(size = 24), legend.text = element_text(size = 14, face = "italic"), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"))
a
#ggsave("moisuture_loss_per_day.png", plot = a, width = 10, height = 7, dpi = 300)

#########################
##### only nonwound #####
#########################
nonwo <- read_csv("nonwound.csv", na = "NA")
mod1=glm(DaystoRecov~genotype + tray, data=nonwo, family="quasipoisson")
summary(mod1)

##get Fx, DF, and P-value.
anova(mod1, test = "F")
exp(0.02)

summary_mod1 <- summary(mod1)
#capture.output(summary_mod1, file = "nonwound_daystorec_genotype_glm.txt")
boxplot(nonwo$DaystoRecov ~ nonwo$un_touched * nonwo$genotype, fill = genotype)
b <- ggplot(data=nonwo,aes(x=genotype,y=DaystoRecov, fill = genotype))+geom_violin(trim=F, fill = "lightgrey")+ geom_jitter(width = 0.2, alpha = 0.5)+xlab("")+geom_boxplot(width = 0.1, alpha = 0.2)+xlab("")+ylab("Days to recovery") + theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), axis.title = element_text(size = 24),legend.text = element_text(size = 24, face = "italic"), legend.title = element_text(size = 24), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 24, face = "bold")) + theme(legend.position = "top") + scale_fill_discrete(name = "genotype") + scale_x_discrete(labels = c("A. nemorensis", "A. sagittata")) + ylim(0,18)
b
#ggsave("nonwounded_plants_border.png", plot = b, width = 10, height = 7, dpi = 300)
plant<-aov(DaystoRecov~genotype, data=nonwo)
summary(plant)
##plot(TukeyHSD(aov(mod1)))
hsd_res3=HSD.test(mod1, c("genotype"), group = TRUE, console = TRUE)
hsd_res3
sigtab3=hsd_res3$groups
sigtab3
sigtab3$factors=sub(":",".",row.names(sigtab3))
sigtab3
ggplot(data=data,aes(x=genotype,y=DaystoRecov, fill = genotype))+geom_violin(trim=F, fill = "lightgrey")+geom_boxplot(width = 0.1, alpha = 0.2)+xlab("Species")+ylab("Days to Recovery") + theme(text = element_text(size = 20))+ geom_text(data=sigtab3,aes(x=factors,y=16, label = groups), size=8, inherit.aes = F) + theme(legend.position = "top")


##############################
###### Rosette area ##########
##############################
mod2=glm(RA~genotype + tray, data=data, family="quasipoisson")
summary(mod2)

##get Fx, DF, and P-value.
anova(mod2, test = "F")
#exp(-0.24)

###############################
###### RA & LTatWilt ##########
###############################
mod3=glm(RA~genotype*LTatwilt + tray, data=data, family="quasipoisson")
summary(mod3)

##get Fx, DF, and P-value.
anova(mod3, test = "F")
exp(-3.77)

###############################
###### LT & LTatWilt ##########
###############################
LTtpt <- read_csv("LTtimepoint.csv", na = "NA")
mod4<-glm(LT~genotype*timepoint, data = LTtpt, family = "quasipoisson")
summary(mod4)

##get Fx, DF, and P-value.
anova(mod4, test = "F")
exp(0.09)

anova_table <- summary(mod4)[[1]]
#summary_df <- as.data.frame(anova_table)
#write.csv(summary_df, "LTtimepoints_aov.csv")
#tukeyplant <- TukeyHSD(mod4)[[3]]
#write.csv(tukeyplant, "LTtimepoints_tukey.csv")
boxplot(LTtpt$LT~LTtpt$timepoint*LTtpt$genotype, col = c("brown", "blue", "brown", "blue"), main ="Interaction of genotypes at different timepoints")


###############################
###### Only Wounded ###########
###############################
#wo <- read_csv("wound.csv", na = "NA")
#mod5=glm(DaystoRecov~genotype + tray, data=wo, family="quasipoisson")
#summary(mod5)

##get Fx, DF, and P-value.
#anova(mod5, test = "F")
#exp(-0.35)

summary_mod5 <- summary(mod5)
capture.output(summary_mod5, file = "wound_daystorec_genotype_glm.txt")
means <- tapply(wo$DaystoRecov, wo$surv, mean)
#ggplot(data=wo,aes(x=genotype,y=DaystoRecov, fill = genotype))+geom_violin(trim=F, fill = "lightgrey")+geom_boxplot(width = 0.1, alpha = 0.2)+geom_jitter(width = 0.2, alpha = 0.5)+xlab("")+ylab("Days to recovery") + theme(text = element_text(size = 20)) + theme(legend.position = "top") + scale_x_discrete(labels = c("A. nemorensis", "A. sagittata")) + ylim(0,10)
#boxplot(DaystoRecov~genotype*surv, data=wo, fill= genotype)
c <- ggplot(data=wo,aes(x=genotype,y=DaystoRecov, fill = genotype))+geom_violin(trim=F, fill = "lightgrey")+ geom_jitter(width = 0.2, alpha = 0.5)+xlab("")+geom_boxplot(width = 0.1, alpha = 0.2)+xlab("")+ylab("Days to recovery") + theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), axis.title = element_text(size = 24), legend.text = element_text(size = 24, face = "italic"), legend.title = element_text(size = 24), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 24, face = "bold"), legend.position = "top") + scale_x_discrete(labels = c("wounded A. nemorensis", "wounded A. sagittata")) + scale_y_continuous(limits = c(0, 12))
c
#ggsave("wounded_plants_broder.png", plot = c, width = 10, height = 7, dpi = 300)


###################################
###### wound & nonwound ###########
###################################
##### Interaction b/w survival and touched un_touched , just model; (no need to plot)
mod6=glm(surv~un_touched * genotype + tray, data=data, family="binomial")
summary(mod6)

##get Fx, DF, and P-value.
anova(mod6, test = "F")
exp(-1.16)

summary_mod6 <- summary(mod6)
capture.output(summary_mod6, file = "surv_un_touched_glm.txt")
coef_table <- summary_mod6$coefficients
summary_df <- as.data.frame(coef_table)
#write.csv(summary_df, "surv_un_touched_glm.csv")


###################################
###### RA & DaystoWilt ############
###################################

#### Interaction b/w genotypes for days to wilting
mod7=glm(DaystoWilt~genotype*RA + tray, data=data, family="quasipoisson")
summary(mod7)

##get Fx, DF, and P-value.
anova(mod7, test = "F")
exp(0.22)


summary_mod7 <- summary(mod7)
capture.output(summary_mod7, file = "RA_glm.txt")
coef_table <- summary_mod7$coefficients
summary_df <- as.data.frame(coef_table)
#write.csv(summary_df, "RA_glm.csv")


###############################################################
###### Frequencing distribution of days to wilting ############
###############################################################
asd = count(data, 'DaystoWilt[101:200]')
zzy = count(data, 'DaystoWilt[1:100]')
freq <- read_csv("freq.csv")
ggbarplot(freq, "DaystoWilt", "freq", fill = "geno", position = position_dodge(0.9))

#############################################
###### survival of the genotypes ############
#############################################
mod8=glm(surv~genotype + tray, data=data, family="binomial")
summary(mod8)

##get Fx, DF, and P-value.
anova(mod8, test = "F")
exp(2.03)

summary_mod8 <- summary(mod8)
capture.output(summary_mod4, file = "surv_genotype_glm.txt")

#############################################
###### days to wilting & RA #################
#############################################

mod9=glm(DaystoWilt ~ genotype * RA + tray , data=data, family="quasipoisson")
summary(mod9)

##get Fx, DF, and P-value.
anova(mod9, test = "F")
exp(0.22)

summary_mod9 <- summary(mod9)
capture.output(summary_mod9, file = "DaystoWilt_RA_lm_glm.txt")
d <- ggplot(data=data, aes(x=RA, y=DaystoWilt, col=genotype)) + geom_point(size = 5, shape = 20, alpha = 0.4) + 
  geom_smooth(method= "lm", size = 3, alpha = 0.5) + xlab("Rosette areas (mm2)") +ylab("Days to wilting") + 
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"),
        axis.title = element_text(size = 24), legend.text = element_text(size = 24, face = "italic"), 
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"), 
        legend.position = c(0.85, 0.85), legend.background = element_blank(), legend.title = element_text(size = 24)) 
d

ggsave("Daystowilt_RA.png", plot = d, width = 10, height = 7, dpi = 300)
ggsave("Daystowilt_RA.pdf", plot = d, width = 10, height = 7, dpi = 300)
ggsave("Daystowilt_RA.svg", plot = d, width = 10, height = 7, dpi = 300)

########################################
###### days to wilting #################
########################################
mod10=glm(DaystoWilt ~ genotype + tray , data=data, family="quasipoisson")
summary(mod10)

##get Fx, DF, and P-value.
anova(mod10, test = "F")
exp(0.0319)

f <- ggplot(data=data,aes(x=genotype,y=DaystoWilt, fill = genotype))+geom_violin(trim=T, fill = "lightgrey")+ geom_jitter(width = 0.2, alpha = 0.5)+xlab("")+geom_boxplot(width = 0.1, alpha = 0.2)+xlab("")+ylab("Days to wilting") + theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), axis.title = element_text(size = 24), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 24, face = "bold"), legend.position = "top") + scale_x_discrete(labels = c("A. nemorensis", "A. sagittata")) + scale_y_continuous(limits = c(0, 20))


f <- ggplot(data=data,aes(x=genotype,y=DaystoWilt, fill = genotype)) + 
  geom_violin(trim=T, alpha = 0.3)+geom_boxplot(width = 0.1, alpha = 0.5)+
  xlab("Genotype")+ylab("Days to wilting") + scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) +
  geom_jitter(width = 0.2) + theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), 
                                   axis.title = element_text(size = 24), legend.text = element_text("none"), 
                                   panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), 
                                   panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"), 
                                   legend.position = "none") +  ylim(0, 20) # + theme(legend.position = "top")
f




ggsave("Daystowilt_both_species.png", plot = f, width = 10, height = 7, dpi = 300)
ggsave("Daystowilt_both_species.pdf", plot = f, width = 10, height = 7, dpi = 300)
ggsave("Daystowilt_both_species.svg", plot = f, width = 10, height = 7, dpi = 300)
########################################
###### moisture loss per day ###########
########################################
mod11=glm.nb(value ~ genotype * days + tray, data=pheno)
summary(mod11)

##get Fx, DF, and P-value.
anova(mod11, test = "F")
exp(0.018)

#mod11=glm(value ~ genotype * days + tray, data=pheno, family="quasipoisson")
#summary(mod11)

summary_mod11 <- summary(mod11)
capture.output(summary_mod11, file = "moisture_loss_per_day_genotype_glm.txt")
boxplot(pheno$value ~ pheno$genotype * pheno$days)
png("moisture_loss_per_day_genotype_glm.png", width = 800, height = 600)  # Adjust width and height as needed
boxplot(pheno$value ~ pheno$genotype * pheno$days)
dev.off()
########################################
###### RA and days to wilting ##########
########################################
mod12=glm.nb(RA ~ DaystoWilt*genotype + tray, data=data)
summary(mod12)

##get Fx, DF, and P-value.
anova(mod12, test = "F")
exp(-0.021)

#mod12=glm(RA ~ DaystoWilt*genotype + tray, data=data, family="quasipoisson")
#summary(mod12)

summary_mod12 <- summary(mod12)
capture.output(summary_mod12, file = "RA_btw_genotype_glm.txt")
e <- ggplot(data=data,aes(x=genotype,y=RA, fill = genotype))+geom_jitter(width=0.1, height = 0.02)+geom_boxplot(width = 0.1, alpha = 0.2)+ylab("rosette area (mm2)") + theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), axis.title = element_text(size = 24),legend.title = element_text(size = 24), legend.text = element_text(size = 24, face = "italic"), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 24, face = "bold"), legend.position = "top", panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) + scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) + scale_fill_discrete(name = "genotype")
e
#ggsave("initial_rosette_area1.png", plot = e, width = 10, height = 7, dpi = 300)

mod13=glm(DaystoWilt ~ RA * genotype + tray, data=data , family="quasipoisson")
summary(mod13)

##get Fx, DF, and P-value.
anova(mod13, test = "F")
exp(0.22)

summary_mod13 <- summary(mod13)
capture.output(summary_mod13, file = "RA_days_to_wilting_glm.txt")
asd <- ggplot(data=data, aes(x=RA, y=DaystoWilt, col=genotype)) + geom_point() + geom_smooth(method= "lm") + theme_bw()+xlab("Rosette areas (mm2)") +ylab("Days to wilting") + theme(text = element_text(size = 20)) + theme(legend.position = "top")
asd
#ggsave("Daystowilt_RA1.png", plot = asd, width = 10, height = 7, dpi = 300)

########################################
###### correlation networks ############
########################################
## Wheth Days to wilting influenced by RA
network_sag <- read_csv("phenetwork_sag.csv", na = "NA")
network_nem <- read_csv("phenetwork_nem.csv", na = "NA")
##############################################################
###### Interaction b/w genotypes for suvival only ############
##############################################################
###### Interaction b/w genotypes for suvival only
mod14=glm(surv~genotype + tray, data=data, family="quasibinomial")
summary(mod14)

##get Fx, DF, and P-value.
anova(mod14, test = "F")
exp(2.03)
###################################################
###### Interaction b/w genotype for LT ############
###################################################
LTtpt <- read_csv("LTtimepoint.csv")
mod15=glm(LT~genotype*timepoint + tray, data=LTtpt)
#mod9=glm.nb(lLT)~genotype/timepoint + tray, data=LTtpt)
summary(mod15)
##get Fx, DF, and P-value.
anova(mod15, test = "F")
exp(0.1026)


####### check with negative binomial ####

mod4_nb <- MASS::glm.nb(LT ~ genotype*timepoint + tray, data=LTtpt)

mod4_nb2 <- glm(LT~ genotype*timepoint + tray, data = LTtpt,family= negative.binomial(theta=mod4_nb[["theta"]]))

summary(mod4_nb2)

anova(mod4_nb2, test = "F")


#######################################
###### ggplot boxplot initialLT #######
#######################################
ggplot(data=data,aes(x=genotype,y=LT, fill = genotype))+geom_jitter(width=0.1, height = 0.02)+geom_boxplot(width = 0.1, alpha = 0.2)+xlab("Species")+ylab("Initial leaf thickness") + theme(text = element_text(size = 20)) + theme(legend.position = "top")
boxplot(data$LT ~ data$genotype)
boxplot(LTtpt$LT~LTtpt$timepoint*LTtpt$genotype, col = c("brown", "brown", "darkcyan", "darkcyan"), main ="Interaction of genotypes at different timepoints")

# Save as PDF
pdf("leaf_thick_initial_and_wilting.pdf")
boxplot(LTtpt$LT ~ LTtpt$timepoint * LTtpt$genotype,
        col = c("brown", "brown", "darkcyan", "darkcyan"),
        main = "Interaction of genotypes at different timepoints")
dev.off()

# Save as PNG
png("leaf_thick_initial_and_wilting.png", width = 800, height = 600)
boxplot(LTtpt$LT ~ LTtpt$timepoint * LTtpt$genotype,
        col = c("brown", "brown", "darkcyan", "darkcyan"),
        main = "Interaction of genotypes at different timepoints")
dev.off()

# Save as SVG
svg("leaf_thick_initial_and_wilting.svg", width = 8, height = 6)
boxplot(LTtpt$LT ~ LTtpt$timepoint * LTtpt$genotype,
        col = c("brown", "brown", "darkcyan", "darkcyan"),
        main = "Interaction of genotypes at different timepoints")
dev.off()




# Save as PDF with increased text size
pdf("leaf_thick_initial_and_wilting.pdf", width = 14, height = 9)
par(cex.lab = 2, cex.axis = 2)
boxplot(LTtpt$LT ~ LTtpt$timepoint * LTtpt$genotype,
        col = c("brown", "brown", "darkcyan", "darkcyan"),
        xlab = "Timepoints: Genotype", ylab = "Leaf Thickness")
dev.off()

# Save as PNG with increased text size
png("leaf_thick_initial_and_wilting.png", width = 1400, height = 900, res = 150)
par(cex.lab = 2, cex.axis = 2)
boxplot(LTtpt$LT ~ LTtpt$timepoint * LTtpt$genotype,
        col = c("brown", "brown", "darkcyan", "darkcyan"),
        xlab = "Timepoints: Genotype", ylab = "Leaf Thickness")
dev.off()

# Save as SVG with increased text size
svg("leaf_thick_initial_and_wilting.svg", width = 12, height = 8)
par(cex.lab = 2, cex.axis = 2)
boxplot(LTtpt$LT ~ LTtpt$timepoint * LTtpt$genotype,
        col = c("brown", "brown", "darkcyan", "darkcyan"),
        xlab = "Timepoints: Genotype", ylab = "Leaf Thickness")
dev.off()
################################
###### LT two timepoints #######
################################
LTtpt <- read_csv("LTtimepoint.csv")
#plant<-aov(LT~genotype*timepoint, data = LTtpt)
plant=glm(LT~genotype*LTatwilt + tray, data=data, family="quasipoisson")
summary(plant)
##get Fx, DF, and P-value.
anova(plant, test = "F")
exp(0.1051)



plant_nb <- MASS::glm.nb(LT ~ genotype*timepoint + tray, data=LTtpt)

plant_nb2 <- glm(LT~ genotype*timepoint + tray, data = LTtpt,family= negative.binomial(theta=plant_nb[["theta"]]))

summary(plant_nb2)

anova(plant_nb2, test = "F")


##########################################
###### interacting two models ############
##########################################

mod16=glm(log(LT)~timepoint*genotype + tray, data=LTtpt)
summary(mod16)

##get Fx, DF, and P-value.
anova(mod16, test = "F")
exp(0.102)

mod17=glm(log(LT)~timepoint+genotype + tray, data=LTtpt)
summary(mod17)

##get Fx, DF, and P-value.
anova(mod17, test = "F")
exp(0.001)

anova(mod16, mod17, test="F")
plot(density(LTtpt$LT))

###################################
###### LT and genotype ############
###################################
#mod18=glm(log(LT)~genotype + tray, data=data, family="gaussian")
#summary(mod18)

##get Fx, DF, and P-value.
#anova(mod18, test = "F")
#exp(0.003)

#############################################
###### ggplot InitialLT_LTatwilt ############
#############################################

e <- ggplot(data=data, aes(x=LTatwilt, y=LT, col=genotype)) + geom_point(alpha = 0.5)  + 
  geom_smooth(method= "lm", alpha = 0.6) + theme(axis.text.y = element_text(size = 24), 
                                                 axis.text.x = element_text(size = 24), axis.title = element_text(size = 24),
                                                 legend.text = element_text(size = 24, face = "italic"), 
                                                 panel.background = element_rect(fill = "white"), 
                                                 panel.grid.major = element_line(color = "lightgrey"), 
                                                 panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 24, face = "bold"), 
                                                 legend.position = c(0.25, 0.85), legend.background = element_blank(), 
                                                 legend.title = element_text(size = 24)) + xlab("Leaf thickness at wilting(mm2)")+ylab("Initial Leaf thickness(mm2)") # + ylim(0.00, 0.07)


e 

ggsave("InitialLT_LTatwilt_border.png", plot = e, width = 10, height = 7, dpi = 300)
ggsave("InitialLT_LTatwilt_border.pdf", plot = e, width = 10, height = 7, dpi = 300)
ggsave("InitialLT_LTatwilt_border.svg", plot = e, width = 10, height = 7, dpi = 300)

##############################
###### LTatWilt and Smw ############
##############################
mod19=glm(LTatwilt~Smw*genotype + tray, data=data, family="quasipoisson")
summary(mod19)

##get Fx, DF, and P-value.
anova(mod19, test = "F")
exp(0.132717)


#######################################
###### ggplot smw_LTatwilt ############
#######################################

f <- ggplot(data=data, aes(x=LTatwilt, y=Smw, col=genotype)) + geom_point(alpha = 0.5)  + 
  geom_smooth(method= "lm", alpha = 0.6) + theme(axis.text.y = element_text(size = 24), 
                                    axis.text.x = element_text(size = 24), axis.title = element_text(size = 24),
                                    legend.text = element_text(size = 24, face = "italic"), 
                                    panel.background = element_rect(fill = "white"), 
                                    panel.grid.major = element_line(color = "lightgrey"), 
                                    panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 24, face = "bold"), 
                                    legend.position = c(0.85, 0.85), legend.background = element_blank(), 
                                    legend.title = element_text(size = 24)) + xlab("Leaf thickness at wilting")+ylab("Soil moisture at wilting") # + ylim(0.00, 0.07)
f
ggsave("smw_LTatwilt_border.png", plot = f, width = 7, height = 5, dpi = 300)
ggsave("smw_LTatwilt_border.svg", plot = f, width = 7, height = 5, dpi = 300)
ggsave("smw_LTatwilt_border.pdf", plot = f, width = 7, height = 5, dpi = 300)

####################################
###### degree of damage ############
####################################
mod20=glm(DoD~genotype + tray, data=data, family="quasipoisson")
summary(mod20)

##get Fx, DF, and P-value.
anova(mod20, test = "F")
exp(-0.437)

summary_mod20 <- summary(mod20)
capture.output(summary_mod20, file = "degree_of_damage_glm.txt")

##############################
###### ggplot DoD ############
##############################
ggplot(data=data, aes(x=DoD, y=Smw, col=genotype)) + geom_smooth(method= "lm") + theme(text = element_text(size = 20)) + xlab("Degree of Damage")+ylab("Soil moisture at wilting") # + ylim(0.00, 0.07)
g <- ggplot(data=data,aes(x=genotype,y=DoD, fill = genotype)) + 
  geom_violin(trim=T, alpha = 0.3)+geom_boxplot(width = 0.1, alpha = 0.5)+
  xlab("Genotype")+ylab("Degree of damage") + scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) +
  geom_jitter(width = 0.2) + theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), 
                        axis.title = element_text(size = 24), legend.text = element_text("none"), 
                        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), 
                        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"), 
                        legend.position = "none") +  ylim(0.00, 8.0) # + theme(legend.position = "top")
g

ggsave("degree_damage_jitter.png", plot = g, width = 7, height = 5, dpi = 300)
ggsave("degree_damage_jitter.svg", plot = g, width = 7, height = 5, dpi = 300)
ggsave("degree_damage_jitter.pdf", plot = g, width = 7, height = 5, dpi = 300)


#################################
###### LT at wilting ############
#################################
mod21=glm(LTatwilt~genotype + tray , data=data, family="quasipoisson")
summary(mod21)

##get Fx, DF, and P-value.
anova(mod21, test = "F")
exp(0.19)

#################################################################
###### Survival with interaction of RA times of genotypes #######
#################################################################
mod22=glm(surv~RA+genotype +tray, data=data, family = "quasipoisson")
summary(mod22)

##get Fx, DF, and P-value.
anova(mod22, test = "F")
exp(0.57)

data$surv=as.factor(data$surv)
mod23=glm(RA~genotype*surv + tray, data = data, family = "quasipoisson")
summary(mod23)

##get Fx, DF, and P-value.
anova(mod23, test = "F")
exp(0.28)

plant.aov=aov(RA~genotype*surv, data=data)
TukeyHSD(plant.aov)
mod24=glm(surv~RA*genotype, data=data, family = "binomial")
summary(mod24)

##get Fx, DF, and P-value.
anova(mod24, test = "F")
exp(0.04)

anova(mod24,mod22,mod23, test="F")
plot(mod23)

##correlation b/w Days to recovery and soil moisture at willting
#ggplot(data=data, aes(x=Smw, y=DaystoRecov, col=genotype)) + 
#  geom_point() + geom_smooth(method= "lm") + theme_bw()+xlab("Soil moisture at wilting") 
#+ylab("Days to recovery") + theme(text = element_text(size = 20)) + theme(legend.position = "top")


d <- ggplot(data=data,aes(x=Smw,y=DaystoRecov, color  = genotype))+ geom_smooth(method= "lm", alpha = 0.5) +
  geom_point(alpha = 0.8) +ylab("Days to recovery") + xlab("Soil moisture at wilting") +
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), 
        axis.title = element_text(size = 24), legend.text = element_text(size = 24),
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"),
        legend.position = c(0.85, 0.85), legend.background = element_blank(), legend.title = element_text(size = 24))

d 

ggsave("Smw_DaystoRecov.png", plot = d, width = 7, height = 5, dpi = 300)
ggsave("Smw_DaystoRecov.pdf", plot = d, width = 7, height = 5, dpi = 300)
ggsave("Smw_DaystoRecov.svg", plot = d, width = 7, height = 5, dpi = 300)


#################################################################
###### Interactin between gneotypes for stomata density ########
#################################################################
mod25=glm(SD~genotype*SS + tray, data = data, family = "quasipoisson")
summary(mod25)


##get Fx, DF, and P-value.
anova(mod25, test = "F")
exp(-2.839)

plot(mod25)
boxplot(data$SD ~ data$genotype, col = c("red", "cyan"), main ="Stomata density")
#p <- ggplot(data=data,aes(x=genotype,y=SD, fill = genotype))+geom_boxplot(width = 0.5, alpha = 1)+ylab("Number of Stotmata/mm2") + theme(text = element_text(size = 16)) + theme(legend.position = "top") + scale_x_discrete(labels = c("A. nemorensis", "A.sagittata"))
#ggsave("Stomata_density.png", plot = p, width = 7, height = 5, dpi = 300)
plant.aov=aov(SD~genotype, data=data)
TukeyHSD(plant.aov)

mod26=glm(SD~genotype*SS, data=data, family = "quasipoisson")
summary(mod26)

##get Fx, DF, and P-value.
anova(mod26, test = "F")
exp(1.0096)

mod27=glm(SD~SS+genotype, data=data, family = "quasipoisson")
summary(mod27)

##get Fx, DF, and P-value.
anova(mod27, test = "F")
exp(0.002946)

anova(mod27,mod25,mod26, test="F")
plot(mod25)

i <- ggplot(data=data,aes(x=genotype,y=SD, fill = genotype))+geom_jitter(width=0.1, height = 0)+geom_boxplot(width = 0.5, alpha = 0.5)+ ylab("Stomata density/mm2") + theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), axis.title = element_text(size = 24), legend.text = element_text(size = 14, face = "italic"), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"), legend.position = "top")+ scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) 

i <- ggplot(data=data,aes(x=genotype,y=SD, fill = genotype))+
  geom_jitter(width=0.1, height = 0, alpha = 0.5)+geom_boxplot(width = 0.3, alpha = 0.8)+ylab("Stomata density/millimeter2)") +
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), 
        axis.title = element_text(size = 24), legend.text = element_text("none"),
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"),
        legend.position = "none")+ scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) 

i
ggsave("stomata_density.png", plot = i, width = 7, height = 5, dpi = 300)
ggsave("stomata_density.pdf", plot = i, width = 7, height = 5, dpi = 300)
ggsave("stomata_density.svg", plot = i, width = 7, height = 5, dpi = 300)
#############################################################
###### Interactin between gneotypes for stomata size ########
#############################################################
mod28=glm(SS~genotype*SD + tray, data = data, family = "quasipoisson")
summary(mod28)
##get Fx, DF, and P-value.
anova(mod28, test = "F")
exp(0.118)


plot(mod28)
boxplot(data$SS ~ data$genotype)
boxplot(data$SS ~ data$genotype, col = c("red", "cyan"), main ="Stomata density")
j <- ggplot(data=data,aes(x=genotype,y=SS, fill = genotype))+
  geom_jitter(width=0.1, height = 0, alpha = 0.5)+geom_boxplot(width = 0.3, alpha = 0.8)+ylab("Stomata size(micrometer)") +
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), 
        axis.title = element_text(size = 24), legend.text = element_text("none"),
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"),
        legend.position = "none")+ scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) 

j
ggsave("Stomata_size.png", plot = j, width = 7, height = 5, dpi = 300)
ggsave("Stomata_size.pdf", plot = j, width = 7, height = 5, dpi = 300)
ggsave("Stomata_size.svg", plot = j, width = 7, height = 5, dpi = 300)


############################################
###### Interactin between LL and SS ########
############################################
mod29=glm(log(LL)~SS*genotype + tray, data = data, family = "quasipoisson")
summary(mod29)

##get Fx, DF, and P-value.
anova(mod29, test = "F")
exp(-0.330)

k <- ggplot(data=data,aes(x=genotype,y=SS, fill = genotype))+geom_jitter(width=0.1, height = 0)+
  geom_boxplot(width = 0.5, alpha = 0.5)+ylab("Leaf length") + 
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), 
        axis.title = element_text(size = 24), legend.text = element_text(size = 14, face = "italic"), 
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"), 
        legend.position = "top")+ scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) 
k
#ggsave("leaf_length.png", plot = k, width = 7, height = 5, dpi = 300)

############################################################
###### Interactin between genotypes for LdW and Lfw ########
############################################################
##Interaction b/w genotypes for LdW and Lfw
mod30=glm(LdW~LfW*genotype + tray, data = data, family = "quasipoisson")
summary(mod30)

##get Fx, DF, and P-value.
anova(mod30, test = "F")
exp(-2.02)

#ggplot(data=data,aes(x=genotype,y=LfW, fill = genotype))+geom_jitter(width=0.1, height = 0)+geom_boxplot(width = 0.5, alpha = 0.5)+ylab("Leaf length") + theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), axis.title = element_text(size = 24), legend.text = element_text(size = 14, face = "italic"), panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"), legend.position = "top")+ scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) 

#########################################
###### Soil Moisture at Wilting #########
#########################################
#mod31=glm(Smw~ genotype + tray, data=data, family="quasipoisson")
#summary(mod31)

mod31_nb <- MASS::glm.nb(Smw~ genotype + tray, data=data)

mod31_nb2 <- glm(Smw~ genotype + tray, data=data,family= negative.binomial(theta=mod31_nb[["theta"]]))

summary(mod31_nb2)

anova(mod31_nb2, test = "F")


##get Fx, DF, and P-value.
anova(mod31, test = "F")
exp(0.104)

boxplot(data$Smw ~ data$genotype)
h <- ggplot(data=data,aes(x=genotype,y=Smw, fill = genotype))+ 
  geom_violin(fill = "white") + geom_jitter(width=0.1, height = 0.01)+
  geom_boxplot(width = 0.1, alpha = 0.2)+ ylim(0, 0.12) +ylab("Soil moisture at wilting") + 
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), 
        axis.title = element_text(size = 24), legend.text = element_text(size = 14, face = "italic"), 
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), 
        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"), 
        legend.position = "top")+ scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) + 
  scale_fill_discrete(name = "genotype")
h
#ggsave("smw_atwilting.png", plot = h, width = 7, height = 5, dpi = 300)

#############################################################################################################
###### Interaction of wounding plants and species for the days to recovery and the degree of damage #########
#############################################################################################################
wo <- read_csv("wound.csv", na = "NA")
#mod32=glm(DaystoRecov~DoD*genotype + tray, data=data, family="quasipoisson")
#summary(mod32)
##get Fx, DF, and P-value.
#anova(mod32, test = "F")
#exp(-0.34572)

mod32_nb <- MASS::glm.nb(DaystoRecov~DoD*genotype + tray, data=data)

mod32_nb2 <- glm(DaystoRecov~DoD*genotype + tray, data=data,family= negative.binomial(theta=mod32_nb[["theta"]]))

summary(mod32_nb2)

anova(mod32_nb2, test = "F")





v <- ggplot(data=wo, aes(x=DoD, y=DaystoRecov, col=genotype)) + geom_point(size = 5, shape = 20, alpha = 0.4) + 
  geom_smooth(method= "lm", size = 3, alpha = 0.5) + xlab("Degree of damage") +ylab("Days to recovery") + 
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"),
        axis.title = element_text(size = 24), legend.text = element_text(size = 24, face = "italic"), 
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"), 
        legend.position = c(0.25, 0.85), legend.background = element_blank(), legend.title = element_text(size = 24))
v

ggsave("daystorecov_dod_wounds.png", plot = v, width = 7, height = 5, dpi = 300)
ggsave("daystorecov_dod_wounds.pdf", plot = v, width = 7, height = 5, dpi = 300)
ggsave("daystorecov_dod_wounds.svg", plot = v, width = 7, height = 5, dpi = 300)
######################################
###### LTat wilting and genotype #####
######################################
mod33=glm(log(LTatwilt)~genotype + tray, data=data, family="gaussian")
summary(mod33)

##get Fx, DF, and P-value.
anova(mod33, test = "F")
exp(0.191)
#########################################
###### Soil Moisture at Wilting vs Smw #########
#########################################
mod34=glm(Smw ~ DaystoWilt + tray, data=data, family="quasipoisson")
summary(mod34)
##get Fx, DF, and P-value.
anova(mod34, test = "F")
exp(-0.183764)
boxplot(data$Smw ~ data$DaystoWilt * data$genotype)

d <- ggplot(data=data,aes(y=Smw,x=DaystoWilt, color  = genotype))+ geom_smooth(method= "lm", alpha = 0.5) +
  geom_point(alpha = 0.8) +xlab("Days to wilting") + ylab("Soil moisture at wilting") +
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), 
        axis.title = element_text(size = 24), legend.text = element_text(size = 24),
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"),
        legend.position = c(0.85, 0.85), legend.background = element_blank(), legend.title = element_text(size = 24))

d 


ggsave("smw_atwilting_DtoW.png", plot = d, width = 7, height = 5, dpi = 300)
ggsave("smw_atwilting_DtoW.pdf", plot = d, width = 7, height = 5, dpi = 300)
ggsave("smw_atwilting_DtoW.svg", plot = d, width = 7, height = 5, dpi = 300)






h <- ggplot(data=data,aes(x=genotype,y=Smw, fill = genotype))+ geom_violin(alpha = 0.1) + 
  geom_jitter(width=0.1, height = 0.01, alpha = 0.5)+geom_boxplot(width = 0.1, alpha = 0.8)+ 
  ylim(0, 0.12) +ylab("Soil moisture at wilting") +
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"), 
        axis.title = element_text(size = 24), legend.text = element_text("none"), 
        panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "lightgrey"), 
        panel.grid.minor = element_line(color = "lightgrey"), plot.title = element_text(size = 20, face = "bold"),
        legend.position = "none")+ scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) + 
  scale_fill_discrete(name = "none")

h
ggsave("smw_atwilting.png", plot = h, width = 7, height = 5, dpi = 300)
ggsave("smw_atwilting.pdf", plot = h, width = 7, height = 5, dpi = 300)
ggsave("smw_atwilting.svg", plot = h, width = 7, height = 5, dpi = 300)
#########################################
###### wound and non wound together #########
#########################################
setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/")
wo_non <- read_csv("wound_nonwound.csv", na = "NA")
mod35=glm(surv ~ un_touched * genotype + tray, data=wo_non, family="quasipoisson")
summary(mod35)
##get Fx, DF, and P-value.
anova(mod35, test = "F")
exp(- 0.39264)
boxplot(wo_non$surv ~ wo_non$un_touched * data$genotype)
x <- ggplot(data=wo_non,aes(x=genotype,y=surv, color = genotype)) + geom_violin(alpha = 0.3) + 
  geom_jitter(width = 0.02, shape = 21)+ ylim(0,1.05) +ylab("Recovery") + 
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24,
        face = "italic"), axis.title = element_text(size = 24), 
        legend.text = element_text(size = 14, face = "italic"), 
        panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_line(color = "lightgrey"), 
        panel.grid.minor = element_line(color = "lightgrey"), 
        plot.title = element_text(size = 20, face = "bold"), 
        legend.position = "none") + 
  scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) + 
  scale_fill_discrete(name = "genotype") + border(color = "darkgrey")
x
ggsave("wound_nonwound_together.png", plot = x, width = 7, height = 5, dpi = 300)
ggsave("wound_nonwound_together.pdf", plot = x, width = 7, height = 5, dpi = 300)
ggsave("wound_nonwound_together.svg", plot = x, width = 7, height = 5, dpi = 300)
#########################################
###### only non-wounding #########
#########################################
non_wo <- read_csv("nonwound_only.csv")
mod36=glm(surv ~ genotype + tray, data=non_wo, family="quasipoisson")
summary(mod36)

##get Fx, DF, and P-value.
anova(mod36, test = "F")
exp(0.67387)

boxplot(data$Smw ~ data$DaystoWilt * data$genotype)
h <- ggplot(data=non_wo,aes(x=genotype,y=surv, color = genotype)) + geom_violin(alpha = 0.3) + 
  geom_jitter(width = 0.02, shape = 21)+ ylim(0,1.05) +ylab("Recovery") + 
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24,
                                                                          face = "italic"), axis.title = element_text(size = 24), 
        legend.text = element_text(size = 14, face = "italic"), 
        panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_line(color = "lightgrey"), 
        panel.grid.minor = element_line(color = "lightgrey"), 
        plot.title = element_text(size = 20, face = "bold"), 
        legend.position = "none") + 
  scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) + 
  scale_fill_discrete(name = "genotype") + border(color = "darkgrey")
  
h
ggsave("nonwound_only.png", plot = h, width = 7, height = 5, dpi = 300)
ggsave("nonwound_only.pdf", plot = h, width = 7, height = 5, dpi = 300)
ggsave("nonwound_only.svg", plot = h, width = 7, height = 5, dpi = 300)
#########################################
###### only RA #########
#########################################

###Rosette area 
mod37=glm(RA~genotype + tray, data=data, family="quasipoisson")
summary(mod37)

##get Fx, DF, and P-value.
anova(mod37, test = "F")
exp(-0.240269)


e <- ggplot(data=data,aes(x=genotype,y=RA, fill = genotype)) +
  geom_jitter(width=0.1, height = 0.02, alpha = 0.6)+geom_boxplot(width = 0.4, alpha = 0.5) +
  ylab("rosette area (mm2)") + 
  theme(axis.text.y = element_text(size = 24), axis.text.x = element_text(size = 24, face = "italic"),
        axis.title = element_text(size = 24),legend.title = element_text(size = 24), 
        legend.text = element_text(size = 24, face = "italic"), panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "lightgrey"), panel.grid.minor = element_line(color = "lightgrey"), 
        plot.title = element_text(size = 24, face = "bold"), legend.position = "none") +
  scale_x_discrete(labels = c("A. nemorensis", "A.sagittata")) + 
  scale_fill_discrete(name = "genotype")
e

ggsave("rosette_area_mm2.png", plot = e, width = 7, height = 5, dpi = 300)
ggsave("rosette_area_mm2.pdf", plot = e, width = 7, height = 5, dpi = 300)
ggsave("rosette_area_mm2.svg", plot = e, width = 7, height = 5, dpi = 300)
#################################################################
###### required packages for phenotype correlation network ######
#################################################################
library(DiagrammeR)
#install.packages("ggdag")
library(ggdag)
#dagify(RA~genotype, data= data) %>% ggdag()
library(psych)
library(qgraph)
library(vegan)
library(heatmaply)
library(plotly)
library(ggcorrplot)
library(circlize)
library(psych)
############################################
###### phenotype correlation networks ######
############################################
network_sag <- read_csv("phenetwork_sag.csv", na = "NA")
network_sag$survival=as.numeric(network_sag$survival)
cormat1 <- cor(network_sag[6:13], method = "pearson", use = "na.or.complete")
qgraph(cormat1,graph = "cor", layout = "circle", sampleSize = nrow(network_sag), alpha = 1.2, cut = 0.1, bonf = T, title = "Phenotype Network of A. sagittata", title.cex = 2, vsize = 9, usePCH = T, details = T,threshold="bonferroni", height = 30, width = 40, label.scale = T, theme ="TeamFortress") ## Add "filetype='png'" to save the plot in Png

network_nem <- read_csv("phenetwork_nem.csv", na = "NA")
#network_nem$genotype=as.numeric(network_nem$genotype)
cormat2 <- cor(network_nem[6:13], method = "pearson", use = "na.or.complete")
qgraph(cormat2,graph = "cor", layout = "circle", sampleSize = nrow(network_nem), alpha = 1.2, cut = 0.1, bonf = T, title = "Phenotype Network of A. nemorensis", title.cex = 2, vsize = 7, usePCH = T, details = T, threshold="bonferroni", height = 30, width = 40, label.scale = T, theme ='TeamFortress') ## Add "filetype='png'" to save the plot in Png

# Assuming `cor_table1` and `cor_table2` are your two correlation tables
dist_mat1 <- 1 - cormat1
dist_mat2 <- 1 - cormat2
## Note: We subtract the correlation values from 1 to convert them into distances.
mantel_result1 <- mantel(cormat1, cormat2, method = "pearson")
print(mantel_result1)
## or
mantel_result2 <- mantel(dist_mat1, dist_mat2, method = "pearson")
print(mantel_result2)
###############################################################################
###### phenotype correlation heatmap for A. sagittata and A. nemorensis #######
###############################################################################
hmp <- read_csv("heatmap_sag.csv", na = "NA")
cormat <- cor(hmp[3:15], method = "pearson", use = "na.or.complete")
m <- ggcorrplot(cormat,
                colors = c("#0571B0", "white", "#CA0020"),
                outline.col = "white",
                hc.order = TRUE,
                lab = T,
                lab_col = "#CCCCCC",
                lab_size = 3,
                sig.level = 0.05,
                tl.cex = 9,
                tl.col = "#636363",
                tl.srt = 45)

m
#ggsave("pearson_corr_sag.png", plot = m, width = 7, height = 5, dpi = 300)

hmp1 <- read_csv("heatmap_nem.csv", na = "NA")
cormat <- cor(hmp1[2:13], method = "pearson", use = "na.or.complete")
n <- ggcorrplot(cormat,
                colors = c("#0571B0", "white", "#CA0020"),
                outline.col = "white",
                hc.order = TRUE,
                lab = T,
                lab_col = "#CCCCCC",
                lab_size = 3,
                sig.level = 0.05,
                tl.cex = 9,
                tl.col = "#636363",
                tl.srt = 45)

n
#ggsave("pearson_corr_nem.png", plot = n, width = 7, height = 5, dpi = 300)


