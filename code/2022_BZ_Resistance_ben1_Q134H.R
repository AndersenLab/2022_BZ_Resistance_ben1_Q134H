library(devtools)
devtools::install_github("AndersenLab/easysorter")
devtools::install_github("kassambara/ggpubr")
library(easysorter)
library(dplyr)
library(broom)
library(ggplot2)
library(lme4)
library(multcompView)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(glue)
#Set variable for data directories
dirs <- c("d:/20220509_ben1a_assay","d:/20220509_ben1b_assay")
#EasySorter pipeline processing
raw <- easysorter::read_data(dirs)
raw_nocontam<- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, directories = TRUE, quantiles = TRUE)
summedplate <- subset(summedraw,select=-c(norm.n))
biopruned <- bioprune(summedplate)
out_pruned <- prune_outliers(biopruned)
Q134Hdata<- regress(out_pruned)

#Filter for meadian.EXT
Q134HmedianEXT <- filter(Q134Hdata,trait=="median.EXT")

#Filter for paper figure
Q134HFigure <- filter(Q134HmedianEXT,strain != "ECA3358")

#Figure

#Sets colors for figure
colsnewben1 <- c("N2" = "orange","ECA882"="grey","ECA3353"="blue")

statsfbzd <- Q134HFigure %>%
  dplyr::filter(condition == "Albendazole")%>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group2 == 'N2')

#Creates plot of results
summary_plot <- Q134HFigure %>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::mutate(strain=factor(strain, levels = c("N2", "ECA882","ECA3353")))%>%
  ggplot+
  aes(x=strain, y=phenotype)+
  geom_jitter(width=0.1, size = 0.3)+
  geom_boxplot(aes(fill = strain),alpha = 0.4,outlier.shape = NA)+
  xlab(" ")+
  ylab("Regressed Animal Length")+
  scale_x_discrete(labels=c("N2" = "WT", "ECA882" = "Deletion","ECA3353" = "Q134H"))+
  stat_pvalue_manual(statsfbzd, label = "p.adj.signif",y.position = c(220),xmax="group1", remove.bracket = TRUE)+
  scale_fill_manual(name= " ", labels = c("N2" = "WT", "ECA882" = "Deletion","ECA3353" = "Q134H"),values = colsnewben1)+
  cowplot::theme_cowplot(12)+
  theme(plot.background = element_rect(fill="white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())
ggsave(filename="20220513Q134Hfigure.png",plot=last_plot(),scale=1, height=2.5, width=2.5, units="in",dpi=300)



#Supplemental Figure

#Sets colors for figure
colsnewben2 <- c("N2" = "orange","ECA882"="grey","ECA3353"="blue","ECA3358"="red")

#Tukey's HSD
statsfbzd2 <- Q134HmedianEXT %>%
  dplyr::filter(condition == "Albendazole")%>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group2 == 'N2')

#Creates plot of results
summary_plot2 <- Q134HmedianEXT %>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::mutate(strain=factor(strain, levels = c("N2", "ECA882","ECA3353","ECA3358")))%>%
  ggplot+
  aes(x=strain, y=phenotype)+
  geom_jitter(width=0.1, size = .3)+
  geom_boxplot(aes(fill = strain),alpha = 0.4,outlier.shape = NA)+
  xlab(" ")+
  ylab("Regressed Animal Length")+
  scale_x_discrete(labels=c("N2" = "WT", "ECA882" = "Deletion","ECA3353" = substitute(paste(italic("ean243"))),"ECA3358"= substitute(paste(italic("ean244")))))+
  stat_pvalue_manual(statsfbzd2, label = "p.adj.signif",y.position = c(220),xmax="group1", remove.bracket = TRUE)+
  scale_fill_manual(name = " ", labels = c("N2" = "WT", "ECA882" = "Deletion","ECA3353" = substitute(paste(italic("ean243"))),"ECA3358"= substitute(paste(italic("ean244")))),values = colsnewben2)+
  cowplot::theme_cowplot(10)+
  theme(plot.background = element_rect(fill="white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank()) 
ggsave(filename="20220513Q134HSupplementalfigure.png",plot=last_plot(),scale=1, height=2.5, width=2.5, units="in",dpi=300)



