## Statistical analyses and plots for the paper:
# Performance of DNA metabarcoding, barcoding and morphology approaches in the identification of insect biodiversity
# by Romana Salis and Nikolaj Gubonin

# CONTENTS:
### Literature search:
## Fig 4C
## Order summary
## LM number of species / specimens vs congruence between morphology and barcoding

### Barcoding:
## chi square test of the matching success
## Fig 2

### Barcoding Vietnamese Butterfly subset:
## Similarity thresholds: Chi square test
## Fig 5A
## Vietnam Butterfly families: Chi square test
## Fig 5B
## Fig 5 combined
## Intraspecific distance vs km distance LM
## Fig S6
## Intraspecific distance vs sampling area LM and plot

### BINs vs morphospecies:
## Metabarcoding LM
## Vietnamese butterfly families LM
## Barcoding datasets LM
## Fig 3

### Metabarcoding:
## Fig S8 venn diagram
## Replicates vs individuals LM
## Reads vs individuals LM
## Fig_S7
## Species accumulation curve for the number of morphospecies
## Species accumulation curve for the number of BINs


library(ggplot2)
library(scales)
library(dplyr)
library(plyr)

setwd("~/Desktop/barcoding_Nikolaj/Insect_barcoding_analyses")

### literature search
df <- read.csv("literature_search_210624.csv")
head(df)
df_success_db <- subset(df, success_db != "NA")
counts = df %>% group_by(order) %>% tally
counts$n <- paste("n =", counts$n, sep=" ")
counts_db = df_success_db %>% group_by(order) %>% tally
counts_db$n <- paste("n =", counts_db$n, sep=" ")
our_data <- subset(df[103:107,])

## Figure 4C
ggplot(df, aes(x=order, y=success, fill=order)) + 
  geom_boxplot() + scale_y_continuous(labels = scales::percent)+ 
  geom_jitter(data=our_data, aes(x=order, y=success), position=position_jitter(0), color="red", size=2.5, pch=20)+
  labs(x = "", y = "Morphospecies identification success via tree- and distance-based methods")+
  coord_flip()+
  geom_text(data = counts, aes(y = 0, label = n), nudge_y=1.1, size=3)+
  theme_classic()+
  theme(legend.position="none")
ggsave("Fig_4C.pdf", width = 15, height = 9, units = "cm")
ggplot(df_success_db, aes(x=order, y=success, fill=order)) + 
  geom_boxplot() + scale_y_continuous(labels = scales::percent)+ 
  geom_jitter(data=our_data, aes(x=order, y=success), position=position_jitter(0), color="red", size=2.5, pch=20)+
  labs(x = "", y = "Morphospecies identification success via barcode database comparison")+
  coord_flip()+
  geom_text(data = counts_db, aes(y = 0, label = n), nudge_y=1.1, size=3)+
  theme_classic()+
  theme(legend.position="none")

## Order summary
sumdata <- ddply(df, "order", summarise,
                 N    = length(success),
                 mean = mean(success),
                 sd   = sd(success),
                 se   = sd / sqrt(N)
)
sumdata

## LM number of species / specimens vs congruence between morphology and barcoding
df$log_species <- log(df$species)
df$log_specimen <- log(df$specimen)
plot(species ~ success, data=df)
plot(specimen ~ success, data=df)
plot(log_species ~ success, data=df)
plot(log_specimen ~ success, data=df)
species.lm <- lm(success ~ log_species, data = df)
summary(species.lm)
anova(species.lm)
specimen.lm <- lm(success ~ log_specimen, data = df)
summary(specimen.lm)
anova(specimen.lm)


### Barcoding
## chi square test of the matching success
library(tidyverse)
data_barcoding <- read.csv("barcoding_datasets.csv")
contingency_table <- data_barcoding %>%
  select(Dataset, Category, count) %>%
  tidyr::spread(key = Category, value = count)
rownames(contingency_table) <- contingency_table$Dataset
contingency_table <- contingency_table[, -1]
chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)
data_barcoding_expected <- as.data.frame(chi_square_result$expected)
data_barcoding_Exp <- data_barcoding_expected%>% rownames_to_column() %>% pivot_longer(cols=-rowname) %>% 
  mutate(facet_cat = "Expected") 
colnames(data_barcoding_Exp) <- c("Dataset", "Category","count","facet_cat")
data_barcoding_observed <- data_barcoding[,-4]
data_barcoding_Obs <- data_barcoding_observed %>% mutate(facet_cat = "Actual")
ggplot(mapping=aes(fill=Category,y=count, x=Dataset,group=Category))+
  theme_light() +
  geom_col(data = data_barcoding_Obs,position="dodge") +
  geom_col(data = data_barcoding_Exp,position="dodge",fill="white",color="black",size=1,alpha=.3)+
  theme(legend.title=element_blank())+
  scale_fill_brewer(palette="Set2") +
  ylab("Counts (specimens)")+
  xlab("")+
  #labs(title=paste("Datasets, 97% threshold, match against BOLD database"))+
  theme(plot.title = element_text(color = "black", size = 12,
                                  face = "bold", hjust = 0.5))

## Fig 2
devtools::install_github("coolbutuseless/ggpattern")
library(ggpattern)
ggplot(data_barcoding, aes(fill=Category, y=percent, x=Dataset, pattern = Category)) + 
  geom_bar_pattern(position="stack", stat="identity", 
                   color = "black", # Adding black borders
                   pattern = ifelse(data_barcoding$Category == "No match", "crosshatch", "none"),
                   pattern_density = 0.05, 
                   pattern_angle = 45,
                   pattern_fill = "black",
                   pattern_spacing = 0.025,
                   width = 0.5) +
  scale_fill_manual(values = c("Correct" = "black", 
                               "Incorrect" = "white", 
                               "Multiple" = "darkgrey", 
                               "No match" = "white")) +
  scale_pattern_manual(values = c("Correct" = "none", 
                                  "Incorrect" = "none", 
                                  "Multiple" = "none", 
                                  "No match" = "crosshatch")) +
  theme_light() +
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(color = "black", size = 12, 
                                  face = "bold", hjust = 0.5)) +
  theme(legend.key.size = unit(1, 'cm'))+
  xlab("") +
  ylab("% specimens")
ggsave("Fig_2.pdf", width = 15, height = 10, units = "cm")


### Barcoding Vietnamese Butterfly subset 
## Similarity thresholds: Chi square test
data_barcoding_VBthreshold <- read.csv("barcoding_VB_thresholds.csv")
data_barcoding_VBthreshold$Dataset <- factor(data_barcoding_VBthreshold$Dataset, 
                                                 levels = c("97%", "98%", "99%", "100%"))
contingency_table_VBt <- data_barcoding_VBthreshold %>%
  select(Dataset, Category, count) %>%
  tidyr::spread(key = Category, value = count)
rownames(contingency_table_VBt) <- contingency_table_VBt$Dataset
contingency_table_VBt <- contingency_table_VBt[, -1]
chi_square_result_VBt <- chisq.test(contingency_table_VBt)
print(chi_square_result_VBt)
data_barcoding_VBthreshold_expected <- as.data.frame(chi_square_result_VBt$expected)
data_barcoding_VBthreshold_Exp <- data_barcoding_VBthreshold_expected%>% rownames_to_column() %>% pivot_longer(cols=-rowname) %>% 
  mutate(facet_cat = "Expected") 
colnames(data_barcoding_VBthreshold_Exp) <- c("Dataset", "Category","count","facet_cat")
data_barcoding_VBthreshold_observed <- data_barcoding_VBthreshold[,-4]
data_barcoding_VBthreshold_Obs <- data_barcoding_VBthreshold_observed %>% mutate(facet_cat = "Actual")
data_barcoding_VBthreshold_Obs$Dataset <- factor(data_barcoding_VBthreshold_Obs$Dataset, 
                                                 levels = c("97%", "98%", "99%", "100%"))
data_barcoding_VBthreshold_Exp$Dataset <- factor(data_barcoding_VBthreshold_Exp$Dataset, 
                                                 levels = c("97%", "98%", "99%", "100%"))
ggplot(mapping=aes(fill=Category,y=count, x=Dataset,group=Category))+
  theme_light() +
  geom_col(data = data_barcoding_VBthreshold_Obs,position="dodge") +
  geom_col(data = data_barcoding_VBthreshold_Exp,position="dodge",fill="white",color="black",size=1,alpha=.3)+
  theme(legend.title=element_blank())+
  scale_fill_brewer(palette="Set2") +
  ylab("Counts (specimens)")+
  xlab("Threshold")+
  theme(plot.title = element_text(color = "black", size = 12,
                                  face = "bold", hjust = 0.5))

## Fig S5A
P_tp <- ggplot(data_barcoding_VBthreshold, aes(fill=Category, y=percent, x=Dataset, pattern = Category)) + 
  geom_bar_pattern(position="stack", stat="identity", 
                   color = "black", # Adding black borders
                   pattern = ifelse(data_barcoding_VBthreshold$Category == "No match", "crosshatch", "none"),
                   pattern_density = 0.05, 
                   pattern_angle = 45,
                   pattern_fill = "black",
                   pattern_spacing = 0.025,
                   width = 0.5) +
  scale_fill_manual(values = c("Correct" = "black", 
                               "Incorrect" = "white", 
                               "Multiple" = "darkgrey", 
                               "No match" = "white")) +
  scale_pattern_manual(values = c("Correct" = "none", 
                                  "Incorrect" = "none", 
                                  "Multiple" = "none", 
                                  "No match" = "crosshatch")) +
  # scale_pattern_angle_manual(values = c(0, 45, 0)) +
  theme_light() +
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(color = "black", size = 12, 
                                  face = "bold", hjust = 0.5)) +
  theme(legend.key.size = unit(1, 'cm'))+
  xlab("") +
  ylab("% specimens")
P_tp
ggsave("Fig_S5A_percent.pdf", width = 15, height = 10, units = "cm")
# exclude 100%
data_barcoding_VBn100 <- data_barcoding_VBthreshold %>%
  filter(Dataset != "100%")
contingency_table_VBn100 <- data_barcoding_VBn100 %>%
  select(Dataset, Category, count) %>%
  tidyr::spread(key = Category, value = count)
rownames(contingency_table_VBn100) <- contingency_table_VBn100$Dataset
contingency_table_VBn100 <- contingency_table_VBn100[, -1]
chi_square_result_VBn100 <- chisq.test(contingency_table_VBn100)
print(chi_square_result_VBn100)

## Vietnam Butterfly families: Chi square test
data_barcoding_VBfamilies <- read.csv("barcoding_VB_families.csv")
data_barcoding_VBfamilies$Dataset <- factor(data_barcoding_VBfamilies$Dataset, 
                                                levels = c("Nymphalidae", "Papilionidae","Pieridae","Lycaenidae","Hesperiidae","Riodinidae"))
contingency_table_VBf <- data_barcoding_VBfamilies %>%
  select(Dataset, Category, count) %>%
  tidyr::spread(key = Category, value = count)
rownames(contingency_table_VBf) <- contingency_table_VBf$Dataset
contingency_table_VBf <- contingency_table_VBf[, -1]
chi_square_result_VBf <- chisq.test(contingency_table_VBf)
print(chi_square_result_VBf)
data_barcoding_VBfamilies_expected <- as.data.frame(chi_square_result_VBf$expected)
data_barcoding_VBfamilies_Exp <- data_barcoding_VBfamilies_expected%>% rownames_to_column() %>% pivot_longer(cols=-rowname) %>% 
  mutate(facet_cat = "Expected") 
colnames(data_barcoding_VBfamilies_Exp) <- c("Dataset", "Category","count","facet_cat")
data_barcoding_VBfamilies_observed <- data_barcoding_VBfamilies[,-4]
data_barcoding_VBfamilies_Obs <- data_barcoding_VBfamilies_observed %>% mutate(facet_cat = "Actual")
data_barcoding_VBfamilies_Obs$Dataset <- factor(data_barcoding_VBfamilies_Obs$Dataset, 
                                                 levels = c("Nymphalidae", "Papilionidae","Pieridae","Lycaenidae","Hesperiidae","Riodinidae"))
data_barcoding_VBfamilies_Exp$Dataset <- factor(data_barcoding_VBfamilies_Exp$Dataset, 
                                                 levels = c("Nymphalidae", "Papilionidae","Pieridae","Lycaenidae","Hesperiidae","Riodinidae"))
ggplot(mapping=aes(fill=Category,y=count, x=Dataset,group=Category))+
  theme_light() +
  geom_col(data = data_barcoding_VBfamilies_Obs,position="dodge") +
  geom_col(data = data_barcoding_VBfamilies_Exp,position="dodge",fill="white",color="black",size=1,alpha=.3)+
  theme(legend.title=element_blank())+
  scale_fill_brewer(palette="Set2") +
  ylab("Counts (specimens)")+
  xlab("Family")+
  theme(plot.title = element_text(color = "black", size = 12,
                                  face = "bold", hjust = 0.5))

## Fig S5B
P_fp <- ggplot(data_barcoding_VBfamilies, aes(fill=Category, y=percent, x=Dataset, pattern = Category)) + 
  geom_bar_pattern(position="stack", stat="identity", 
                   color = "black", # Adding black borders
                   pattern = ifelse(data_barcoding_VBfamilies$Category == "No match", "crosshatch", "none"),
                   pattern_density = 0.05, 
                   pattern_angle = 45,
                   pattern_fill = "black",
                   pattern_spacing = 0.025,
                   width = 0.5) +
  scale_fill_manual(values = c("Correct" = "black", 
                               "Incorrect" = "white", 
                               "Multiple" = "darkgrey", 
                               "No match" = "white")) +
  scale_pattern_manual(values = c("Correct" = "none", 
                                  "Incorrect" = "none", 
                                  "Multiple" = "none", 
                                  "No match" = "crosshatch")) +
  # scale_pattern_angle_manual(values = c(0, 45, 0)) +
  theme_light() +
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(color = "black", size = 12, 
                                  face = "bold", hjust = 0.5)) +
  theme(legend.key.size = unit(1, 'cm'))+
  xlab("") +
  ylab("% specimens")
P_fp
ggsave("Fig_S5B_percent.pdf", width = 15, height = 10, units = "cm")

## Fig 5 combined
library(patchwork)
# Combine the plots with a shared legend at the bottom
Fig_S5_combined_plot_count <- P_tc + P_fc + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
Fig_S5_combined_plot_count

Fig_S5_combined_plot_percent <- P_tp + P_fp + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
Fig_S5_combined_plot_percent
ggsave("Fig_S5_combined_plot_percent.pdf", width = 30, height = 15, units = "cm")

## Intraspecific distance vs km distance LM
library(lme4)
library(lmerTest)
library(MuMIn)
packageVersion("lme4")
packageVersion("lmerTest")
packageVersion("MuMIn")

vietnam_distance <- read.csv("vietkordfinal_x.csv")
head(vietnam_distance)

distpx <- lmer(Matchp~dist_km + (1|Currentid), data = vietnam_distance)
r.squaredGLMM(distpx)
summary(distpx)
anova(distpx)

## Fig S6
ggplot(vietnam_distance,aes(dist_km, Matchp)) +
  geom_point(size=0.5) +
  geom_smooth(method='lm', color="black",size=1) +
  ylab("Intraspecific COI similarity (%)")+
  xlab("Geographic distance (km)")+
  theme(plot.title = element_text(color = "black", size = 12,
                                  face = "bold", hjust = 0.5), legend.position = "none")+
  theme_classic()
ggsave("Fig_S6.pdf", width = 15, height = 10, units = "cm")
ggsave("Fig_S6.jpeg", width = 15, height = 10, units = "cm")

## Intraspecific distance vs sampling area LM and plot
df_area <- read.csv("vietnam_area.csv")
head(df_area)
lm_area <- lm(meanintra  ~  AreaGIS_km2, data=df_area)
summary(lm_area)
anova(lm_area)
ggplot(df_area,aes(AreaGIS_km2, meanintra)) +
  geom_point() +
  geom_smooth(method='lm', color="black") +
  theme_classic()


### BINs vs morphospecies
BINvSp <- read.csv("BINsvSp.csv")
## Metabarcoding LM
BINvSp_mb <- subset(BINvSp, method == "metabarcoding")
mb <- lm(BINs  ~  species, data=BINvSp_mb)
summary(mb)
anova(mb)

## Vietnamese butterfly families LM
BINvSp_vb <- subset(BINvSp, method == "barcoding_VB")
vb<- lm(BINs  ~  species, data=BINvSp_vb)
summary(vb)
anova(vb)

## Barcoding datasets LM
BINvSp_b <- subset(BINvSp, method == "barcoding")
b<- lm(BINs  ~  species, data=BINvSp_b)
summary(b)
anova(b)

## Fig 3 
pmb <- ggplot(BINvSp_mb,aes(species, BINs)) +
  geom_point() +
  geom_smooth(method='lm', color="black") +
  xlab("Number of morphologically identified species")+
  ylab("Number of BINs")+
  theme_classic()
pvb <- ggplot(BINvSp_vb,aes(species, BINs)) +
  geom_point() +
  geom_smooth(method='lm', color="black") +
  xlab("Number of morphologically identified species")+
  ylab("Number of BINs")+
  theme_classic()
pb <- ggplot(BINvSp_b,aes(species, BINs)) +
  geom_point() +
  geom_smooth(method='lm', color="black") +
  xlab("Number of morphologically identified species")+
  ylab("Number of BINs")+
  theme_classic()
Fig_3 <- pb + pvb + pmb +
  plot_annotation(tag_levels = 'A', tag_suffix = ')')+
  plot_layout(axis_titles = "collect") & 
  theme(legend.position = "bottom")+
  plot_annotation(tag_levels = 'A')
Fig_3
ggsave("Fig_3_updated.pdf", width = 30, height = 10, units = "cm")


### Metabarcoding
## Fig S8 venn diagram
install.packages("VennDiagram")               
library("VennDiagram")    
grid.newpage()                                        
draw.pairwise.venn(area1 = 97,                        
                   area2 = 97,
                   cross.area = 88,
                   fill = c("#0073C2FF", "#EFC000FF"),
                   lty = "blank",
                   category = c("Morphospecies", "BINs"))

## Replicates vs individuals LM
MetaB <- read.csv("Metabarcoding_Results_CommonSpecies.csv") # use dataset with common species only
head(MetaB)
MetaB$log_ind <- log(MetaB$Individuals)
rep_v_ind <- lm( log_ind ~   Replicates7, data=MetaB)
summary(rep_v_ind)
anova(rep_v_ind)

## Reads vs individuals LM
MetaB$log_reads <- log(MetaB$Reads7)
read_v_ind <- lm(log_ind   ~ log_reads , data=MetaB)
summary(read_v_ind)
anova(read_v_ind)

## Fig_S7
pr <- ggplot(rep_v_ind,aes(Replicates7,log_ind)) +
  geom_point() +
  geom_smooth(method='lm', color="black") +
  ylab("Log number of individuals per species")+
  xlab("Detection (number of replicates)")+
  #ylim(0,7)+
  theme_classic()
pi <- ggplot(read_v_ind,aes(log_reads, log_ind)) +
  geom_point() +
  geom_smooth(method='lm', color="black") +
  ylab("Log number of individuals per species")+
  xlab("Log number of reads")+
  theme_classic()
Fig_S7 <- pr + pi +
  plot_annotation(tag_levels = 'A', tag_suffix = ')')+
  plot_annotation(tag_levels = 'A')
Fig_S7
ggsave("Fig_S9.pdf", width = 20, height = 10, units = "cm")

## Species accumulation curve for the number of morphospecies
library(vegan)
MetaB <- read.csv("Metabarcoding_Results_CommonSpecies.csv")#metabarcoding BIN dataset with just species common with morphology
head(MetaB)
replicate_columns <- MetaB[, grep("Rep", names(MetaB))]
replicate_columns <- replicate_columns[, -10]# Extract relevant columns for replicates
replicate_columns <- replicate_columns[, -9]
replicate_columns <- replicate_columns[, -8]
replicate_columns
replicate_columns_t <- t(replicate_columns)# Transpose the data for vegan (samples in rows, species in columns)
spec_accum_result <- specaccum(replicate_columns_t)# Calculate the species accumulation curve
spec_accum_result
plot(spec_accum_result, xlab="Number of replicates", ylab="Cumulative number of species",
     main="", ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col="lightblue")# Plot the species accumulation curve
lines(spec_accum_result$sites, spec_accum_result$richness, type="b", col="black", lwd=2)# Add error bars for standard deviation

MetaB_BINSp <- read.csv("Metabarcoding_Results_BINSpecies.csv")#full metabarcoding BIN dataset (include unique BINs)
head(MetaB_BINSp)
replicate_columns <- MetaB_BINSp[, grep("Rep", names(MetaB_BINSp))]
replicate_columns <- replicate_columns[, -10]# Extract relevant columns for replicates
replicate_columns <- replicate_columns[, -9]
replicate_columns <- replicate_columns[, -8]
replicate_columns <- replicate_columns[-98,]
replicate_columns
replicate_columns_t <- t(replicate_columns)
spec_accum_result <- specaccum(replicate_columns_t)
spec_accum_result
plot(spec_accum_result, xlab="Number of replicates", ylab="Cumulative number of species",
     main="", ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col="lightblue")
lines(spec_accum_result$sites, spec_accum_result$richness, type="b", col="black", lwd=2)

## Species accumulation curve for number of BINs
MetaB_BINs <- read.csv("Metabarcoding_Results_BINs.csv")
head(MetaB_BINs)
replicate_columns <- MetaB_BINs[, grep("Rep", names(MetaB_BINs))]
replicate_columns <- replicate_columns[, -8]
replicate_columns
replicate_columns_t <- t(replicate_columns)
spec_accum_result <- specaccum(replicate_columns_t)
spec_accum_result
plot(spec_accum_result, xlab="Number of replicates", ylab="Cumulative number of BINs",
     main="", ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col="lightblue")
lines(spec_accum_result$sites, spec_accum_result$richness, type="b", col="black", lwd=2)
