setwd("C:/Users/jipko/OneDrive - UvA/RP Bioclock/Data Analysis")
# Libraries Needed
library(readxl) 
library(vegan)
library(ggplot2)


library(lme4)
library(MASS)
library(brglm2)

FAM <- read_excel("Sample Sheet Week 19-20 May FAM.xlsx")
View(FAM)
str(FAM)
#### Sub-setting the data + adding Diversity indexes ####

FAM <- FAM[,-(35:37)]
FAMshan <- diversity(FAM[, c(13:31)]) # Shannon index is based on randomness present at a site and considers both species richness and equitability in distribution in a sample
FAMsimp <- diversity(FAM[, c(13:31)], index = "simpson") # Simpson index is considered more as a dominance index as it accounts proportion of species in a sample.
FAM <- cbind(FAM, FAMshan, FAMsimp) # binds indexes to the dataset
FAM$`Trap Type` <- as.factor(FAM$`Trap Type`) # force variable to be read as factor
FAM$Area <- as.factor(FAM$Area)
FAM$`Trap Location` <- as.factor(FAM$`Trap Location`)
FAM$Date <- as.factor(FAM$Date)
FAM$`Trap Type` <- relevel(FAM$`Trap Type`, ref = "Light Trap") # re-levels light trap as reference level for visualizations 
str(FAM)
View(FAM)
#### Visualizing the data ####

# Boxplots Total Caught & Shannon index
data_ggp <- data.frame(x = FAM$`Trap Type`, 
                       y = c(FAM$`Macro Total`, FAM$`Micro Total`, FAM$`Total Indiv.`),
                       group = c(rep("Macro", nrow(FAM)),
                                 rep("Micro", nrow(FAM)),
                                 rep("Total", nrow(FAM)))) 
data_ggp$x <- relevel(data_ggp$x, ref = "Light Trap")
# Reshape data frame

BoxplotCounts <- ggplot(FAM, aes(x= FAM$`Trap Type`, y = FAM$`Total Indiv.`)) +
  geom_boxplot() + 
  xlab("Trap Type") + 
  ylab("Number of individuals") + 
  scale_colour_discrete("Moth Category") + 
  theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))
BoxplotCounts

Boxplotshann <- ggplot(FAM, aes(y = FAM$FAMshan, x = FAM$`Trap Type`, fill = FAM$Area)) + 
  geom_boxplot() + 
  xlab("Trap Type") + 
  ylab("Shannon-Weiner-index") + 
  theme(legend.background = element_rect(fill="gray90", size=.5 , linetype="dotted")) + theme(text = element_text(size = 18)) +
  scale_fill_manual(values = c("#01579b", "#ffea00")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank())
Boxplotshann
summary(Boxplotshann)
text = element_text(size = 20)
# Look for Correlations
cor(FAM[, c(8,9,11,12)]) 
# 1) Temperature and Wind speed are negatively correlated
# 2) Wind speed and Cloud cover are negatively correlated

#### Generalized Linear Modeling ####

hist(FAM$`Total Indiv.`) # data is very much zero-inflated
# Zero-inflated poisson regression is used to model count data that has an excess of zero counts, like I do. We can do this with the Generalized Linear model
# Model selection Count Data #
# The function stepAIC will give us the model with the best fit in terms of lowest information loss.
# First make the complex model with our variables of interest; Area, Traptype, Traplocation and Date. Due to co-linearity/ perfectly correlating variables Temperature, Wind & Cloud cover could not be included in the first model since they are co-linear with trapping Date. Trap Location can be included but results in a lot of NA's. NA in this case means that the coefficient is not estimable. This can happen due to exact collinearity, as I've mentioned. But, it can also happen due to not having enough observations to estimate the relevant parameters (e.g. if p>n). If you predictors are categorical and you're adding interaction terms, an NA can also mean that there are no observations with that combination of levels of the factors
# Since the data is very much zero inflated it consists of very low means, which will be problematic in understanding our model estimates and standard errors. This could eventually lead to complete seperation problems and large standard errors and estimates. Therefor I used a bias reduction method from the brglmfit
# Websites helping to understand
# https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faqwhat-is-complete-or-quasi-complete-separation-in-logisticprobit-regression-and-how-do-we-deal-with-them/
# https://stats.stackexchange.com/questions/399403/trouble-modeling-zero-inflated-data-estimates-and-standard-errors-are-off-with 

#### Count Model ####

InteractGLM <- glm(FAM$`Total Indiv.` ~ FAM$`Trap Type` * FAM$DateN * FAM$`Light Polution` * FAM$`Trap Location`, family = "poisson", method="brglmFit")
summary(InteractGLM)
# Let stepAIC do its work
stepAIC(InteractGLM, direction = "both")
# Error due to NA's caused by Trap Location variable
InteractGLM2 <- glm(FAM$`Total Indiv.` ~ FAM$`Trap Type` * FAM$DateN * FAM$`Light Polution`, family = "poisson", method="brglmFit")
summary(InteractGLM2)
stepAIC(InteractGLM2, direction = "both")
# We take the called model as our best fit.
FinalGLM <- glm(formula = FAM$`Total Indiv.` ~ FAM$`Trap Type` * FAM$DateN, family = "poisson", method="brglmFit")
summary(FinalGLM)
# 1) Light traps catch significantly more individuals compared to other trap types
# 2) Trap Area does not seem to make a difference in trap catches.
# 3) Significantly more trap catches over the days? -> Why?

#### PCA on Abiotic Variables ####
POFAM <- read_excel("C:/Users/jipko/OneDrive - UvA/RP Bioclock/Data Analysis/Data for PCA Env_Spec.xlsx")
POFAM$Species <- as.factor(POFAM$Species)
str(POFAM)
View(POFAM)
pcaPO<-rda(POFAM[,c(8,9,11,12)], scale = TRUE)
library(ggfortify)
library(cluster)
pca_res <- prcomp(POFAM[,c(8,9,10,12,13,14)], scale = TRUE)
autoplot(pca_res, data = POFAM, loadings = TRUE, loadings.colour = 'blue', label = TRUE,loadings.label = TRUE, frame = FALSE, frame.colour = "Area", loadings.label.size = 3 + geom_jitter())
autoplot(pca_res, data = POFAM, colour = 'Date2', loadings = FALSE, frame = TRUE, frame.colour = "Date2",loadings.colour = 'blue', label = TRUE, loadings.label = TRUE, loadings.label.size = 3 + geom_jitter( width = 100000000, height = 100000000 ))
summary(pca_res)
pca_res$rotation

POFAM2 <- read_excel("Data for PCA Env_Spec.xlsx")
str(POFAM2)
row.names(POFAM2) <- POFAM2$ID
View(POFAM2)
pca_res2 <- prcomp(POFAM2[,c(8,9,10,12,13)], scale = TRUE)
autoplot(pca_res2, data = POFAM2, loadings = TRUE, loadings.colour = 'blue', label = FALSE,loadings.label = TRUE, frame = FALSE, loadings.label.size = 3)
autoplot(pca_res2, data = POFAM2, colour = 'Date2', loadings = TRUE, frame = TRUE, frame.colour = "Date2",loadings.colour = 'blue', loadings.label = TRUE, ggrepel.max.overlaps = Inf, loadings.label.size = 3)
summary(pca_res)
pca_res$rotation
ggplot(pca_res2)
res.ind <- get_pca_ind(pca_res)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 


#### Significance of Variables / ANOVA & Old School PCA ####
POFAM <- cbind(POFAM, scores(pcaPO)$sites) 
View(POFAM)
POFAMtest <- aov(POFAM$PC1 ~ POFAM$Date * POFAM$Area * POFAM$`Trap Type`)
summary(POFAMtest)

FinalGLM2 <- glm(formula = FAM2$`Total Indiv.` ~ FAM2$`Light Polution` * FAM2$`Wind m/s`, family = "poisson", method="brglmFit")
stepAIC(FinalGLM2, direction = "both")
y <- glm(formula = FAM2$`Total Indiv.` ~ FAM2$`Wind m/s` + FAM2$`Light Polution`, 
         family = "poisson", method = "brglmFit")
summary
x <- glm(formula = FAM2$`Total Indiv.` ~ FAM2$`Wind m/s` * FAM2$`Light Polution`, 
         family = "poisson", method = "brglmFit")
summary(x)


#### Diversity Model ####

# Since our data has excess zeros in terms of counts. It also has excess zero's in terms of diversity indexes; except for one sticky trap only light traps caught more then 1 species. For this reason our diversity analysis will only be done on the light traps
FAM2 <- as.data.frame(FAM[c(1:4, 35:38, 69:72,103:106),])
View(FAM2)
ShannGLM <- glm(FAM2$FAMshan ~ FAM2$`Light Polution` * FAM2$DateN)
summary(ShannGLM)
stepAIC(ShannGLM, direction = "both")
ShannGLM2 <- glm(FAM2$FAMshan ~ FAM2$`Light Polution` + FAM2$DateN)
summary(ShannGLM2)
# 1) There is no significant difference in moth diversity between lightened and darkened areas for Light traps.
# 2) There is no significant difference in moth diversity between days for Light traps.
# 3) Since the other trap types have hardly caught a diverse range of species, light traps are in that sense superior.

#### PERMANOVA ####
# Permutational multivariate analysis of variance (PERMANOVA) is a non-parametric multivariate statistical test. It is used to compare groups of objects and test the null hypothesis that the centroids and dispersion of the groups as defined by measure space are equivalent for all groups.
# PERMANOVA assumes no distribution, allows for differences in between-group variation, is insensitive to multicollinearity, allows for multiple variables and is insensitive to many zeros.
# https://archetypalecology.wordpress.com/2018/02/21/permutational-multivariate-analysis-of-variance-permanova-in-r-preliminary/
FAMdist <- vegdist(FAM2[,c(13:31)], method = "bray")
FAMdist
FAMdiv <- adonis2(FAM2[,c(13:31)] ~ FAM2$DateN * FAM2$Area, permutations = 999, method = "bray")
FAMdiv
# The different areas and trapping dates (Light vs Dark, and day 1:4) indeed have a significant effect on species composition (p < 0.05).

#### PCA ####
View(FAM2)
pcaFAM2 <- rda(FAM2[,c(8:31)], scale = TRUE)
plot(pcaFAM2)
biplot(pcaFAM2)
scores(pcaFAM2)$species
summary(pcaFAM2)
# Better visualizations
library(factoextra)
library(FactoMineR)
FAM3 <- FAM2[ , which(apply(FAM2, 2, var) != 0)]
View(FAM3)
pcaFAM3 <- prcomp(FAM3[,c(6:23)], scale = TRUE) 
summary(pcaFAM3)
fviz_eig(pcaFAM3)
fviz_pca_var(pcaFAM3,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
# Visualize contributions to dimension 1
fviz_contrib(pcaFAM3, choice = "var", axes = 1, top = 20)
# Visualize contributions to dimension 2
fviz_contrib(pcaFAM3, choice = "var", axes = 2, top = 20)
# Table with contributions
res.var <- get_pca_var(pcaFAM3)
contrib <- data.frame(res.var$contrib)# Contributions to the PCs
View(contrib)
# Grouped PCA Visualizations
fviz_pca_ind(pcaFAM3, # Grouped by Area
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = FAM2$Area, # color by groups
             palette = c("#01579b", "#ffea00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")
# The first three principal components (Dim1, Dim2, Dim3) explained 18,5%, 14,1 % and 13,2 % respectively, captured 45,8 % of variance in the moth community data. Dim1 was dominated by Geomitridae (Geometer Moths), which loaded in the negative end. The positive axis was loaded with Sphingidae (sphinx moths) (FIG. ^C). Dim2 was dominated negatively by Pyralidae and positively by Notodontidae. Dim3 was dominated negatively by Gracillariidae and positively by Tortricidae and Erebidae.

#### Significance of Variables / ANOVA & PCA ####

library(dplyr)
get_pca_ind(pcaFAM3)$coord
score <- as_tibble(get_pca_ind(pcaFAM3)$coord) #extract individual scores
View(score)
FAM2 <- cbind(FAM2, score[1:4]) 
View(FAM2)
PCAGLM <- aov(FAM2$Dim.1 ~ FAM2$Area * FAM2$Date) # Simple linear model
stepAIC(PCAGLM, direction = "both")
PCAGLM <- aov(FAM2$Dim.1 ~ FAM2$Area + FAM2$Date)
summary(PCAGLM)                                    # D1 Date Significantly Different
PCAGLM2 <- aov(FAM2$Dim.2 ~ FAM2$Area * FAM2$Date) 
stepAIC(PCAGLM2, direction = "both")
PCAGLM2 <- aov(FAM2$Dim.2 ~ FAM2$Area + FAM2$Date)
summary(PCAGLM2)                                   # D2 Area Significantly Different

PCAGLM3 <- aov(FAM2$Dim.3 ~ FAM2$Area * FAM2$Date) 
stepAIC(PCAGLM, direction = "both")
PCAGLM3 <- aov(FAM2$Dim.3 ~ FAM2$Area + FAM2$Date)
summary(PCAGLM3)                                   # D3 no effect
PCAGLM4 <- aov(FAM2$Dim.4 ~ FAM2$Area * FAM2$Date) 
stepAIC(PCAGLM, direction = "both")                
PCAGLM4 <- aov(FAM2$Dim.4 ~ FAM2$Area + FAM2$Date)
summary(PCAGLM4)                                   # D4 no effect