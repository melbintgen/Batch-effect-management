# Title: Batch effect management practice

#----------------------- Packages installation and loading --------------------#

# CRAN
cran.pkgs <- c('pheatmap', 'vegan', 'ruv', 'ggplot2', 
               'performance', 'gridExtra')

# install.packages(cran.pkgs)

# Bioconductor
bioc.pkgs <- c('mixOmics', 'sva', 'limma', 'Biobase', 'metagenomeSeq', 
               'PLSDAbatch', 'TreeSummarizedExperiment')

# if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(bioc.pkgs)  

# load packages 
suppressMessages(suppressWarnings(sapply(c(cran.pkgs, bioc.pkgs), require, 
                                         character.only = TRUE)))

# print package versions
sapply(c(cran.pkgs, bioc.pkgs), package.version)


#----------------------------- Data pre-processing ----------------------------#

## Pre-filtering
# AD data
data('AD_data') 
ad.count <- assays(AD_data$FullData)$Count
dim(ad.count)

ad.filter.res <- PreFL(data = ad.count)
ad.filter <- ad.filter.res$data.filter
dim(ad.filter)

# zero proportion before filtering
ad.filter.res$zero.prob
# zero proportion after filtering
sum(ad.filter == 0)/(nrow(ad.filter) * ncol(ad.filter))

# extract the metadata
ad.metadata <- rowData(AD_data$FullData)

# extract the batch info
ad.batch = factor(ad.metadata$sequencing_run_date, 
                  levels = unique(ad.metadata$sequencing_run_date))

# extract the treatment info
ad.trt = as.factor(ad.metadata$initial_phenol_concentration.regroup)

# add names on each 
names(ad.batch) <- names(ad.trt) <- rownames(ad.metadata)

#### Exercise 1: How many samples are there within each treatment and batch group?

length(ad.batch)
summary(ad.batch)

length(ad.trt)
summary(ad.trt)

table(ad.batch, ad.trt)

#######################################

## Transformation
ad.clr <- logratio.transfo(X = ad.filter, logratio = 'CLR', offset = 1) 
class(ad.clr) = 'matrix'

#--------------------------- Batch effect detection ---------------------------#

## Principal component analysis (PCA)
ad.pca.before <- pca(ad.clr, ncomp = 3, scale = TRUE)

Scatter_Density(object = ad.pca.before, batch = ad.batch, trt = ad.trt, 
                title = 'AD data', trt.legend.title = 'Phenol conc.')

#### Exercise 2: Interpret the PCA plot created above

## Boxplots 
ad.OTU.name <- selectVar(ad.pca.before, comp = 1)$name[1]

ad.OTU_batch <- data.frame(value = ad.clr[,ad.OTU.name], batch = ad.batch)
head(ad.OTU_batch)

box_plot(df = ad.OTU_batch, title = paste(ad.OTU.name, '(AD data)'), 
         x.angle = 30)

## Density plots
density_plot(df = ad.OTU_batch, title = paste(ad.OTU.name, '(AD data)'))

## statistical tests
# reference level: 14/04/2016
ad.batch <- relevel(x = ad.batch, ref = '14/04/2016')

ad.OTU.lm <- linear_regres(data = ad.clr[,ad.OTU.name], 
                           trt = ad.trt, 
                           batch.fix = ad.batch, 
                           type = 'linear model')
summary(ad.OTU.lm$model$data)

# reference level: 21/09/2017
ad.batch <- relevel(x = ad.batch, ref = '21/09/2017')

ad.OTU.lm <- linear_regres(data = ad.clr[,ad.OTU.name], 
                           trt = ad.trt, 
                           batch.fix = ad.batch, 
                           type = 'linear model')
summary(ad.OTU.lm$model$data)

#### Exercise 3: Interpret the result when "21/09/2017" was set as the reference level

## Heatmap
# scale the clr data on both OTUs and samples
ad.clr.s <- scale(ad.clr, center = TRUE, scale = TRUE)
ad.clr.ss <- scale(t(ad.clr.s), center = TRUE, scale = TRUE)

ad.anno_col <- data.frame(Batch = ad.batch, Treatment = ad.trt)
ad.anno_colors <- list(Batch = color.mixo(seq_len(5)), 
                       Treatment = pb_color(seq_len(2)))
names(ad.anno_colors$Batch) = levels(ad.batch)
names(ad.anno_colors$Treatment) = levels(ad.trt)

pheatmap(ad.clr.ss, 
         cluster_rows = FALSE, 
         fontsize_row = 4, 
         fontsize_col = 6,
         fontsize = 8,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = ad.anno_col,
         annotation_colors = ad.anno_colors,
         border_color = 'NA',
         main = 'AD data - Scaled')

##########################################

## Partial redundancy analysis (pRDA)
ad.factors.df <- data.frame(trt = ad.trt, batch = ad.batch)
head(ad.factors.df)

ad.rda.before <- varpart(ad.clr, ~ trt, ~ batch, 
                         data = ad.factors.df, scale = TRUE)
ad.rda.before$part$indfract

#---------------------------- Managing batch effects --------------------------#

## Accounting for batch effects
### Linear regression: linear model (LM) and linear mixed model (LMM)
ad.lm <- linear_regres(data = ad.clr, 
                       trt = ad.trt, 
                       batch.fix = ad.batch, 
                       type = 'linear model',
                       p.adjust.method = 'fdr')

# p values adjusted for batch effects
ad.p.adj <- ad.lm$adj.p 

check_model(ad.lm$model$OTU12)

head(ad.lm$adj.R2)

head(ad.lm$AIC)

### Remove unwanted variation in 4 steps (RUV4) 
# empirical negative controls
ad.empir.p <- c()
for(e in seq_len(ncol(ad.clr))){
  ad.empir.lm <- lm(ad.clr[,e] ~ ad.trt)
  ad.empir.p[e] <- summary(ad.empir.lm)$coefficients[2,4]
}
ad.empir.p.adj <- p.adjust(p = ad.empir.p, method = 'fdr')
ad.nc <- ad.empir.p.adj > 0.05

# estimate k
ad.k.res <- getK(Y = ad.clr, X = ad.trt, ctl = ad.nc)
ad.k <- ad.k.res$k

# RUV4
ad.ruv4 <- RUV4(Y = ad.clr, X = ad.trt, ctl = ad.nc, k = ad.k) 
ad.ruv4.p <- ad.ruv4$p
ad.ruv4.p.adj <- p.adjust(ad.ruv4.p, method = "fdr")

#######################################

## Correcting for batch effects
### ComBat
# the treatment design matrix
ad.mod <- model.matrix( ~ ad.trt)
ad.ComBat <- t(ComBat(t(ad.clr), batch = ad.batch, 
                      mod = ad.mod, par.prior = FALSE))

### PLSDA-batch
# estimate the number of treatment components
ad.trt.tune <- plsda(X = ad.clr, Y = ad.trt, ncomp = 5)
ad.trt.tune$prop_expl_var #1

# estimate the number of batch components
ad.batch.tune <- PLSDA_batch(X = ad.clr, 
                             Y.trt = ad.trt, Y.bat = ad.batch,
                             ncomp.trt = 1, ncomp.bat = 10)
ad.batch.tune$explained_variance.bat 
sum(ad.batch.tune$explained_variance.bat$Y[seq_len(4)]) #4

ad.PLSDA_batch.res <- PLSDA_batch(X = ad.clr, 
                                  Y.trt = ad.trt, Y.bat = ad.batch,
                                  ncomp.trt = 1, ncomp.bat = 4)
ad.PLSDA_batch <- ad.PLSDA_batch.res$X.nobatch

### Remove Unwanted Variation-III (RUVIII)
ad.replicates <- ad.metadata$sample_name.data.extraction
head(table(ad.replicates, ad.batch))

ad.replicates.matrix <- replicate.matrix(ad.replicates)

ad.RUVIII <- RUVIII(Y = ad.clr, M = ad.replicates.matrix, 
                    ctl = ad.nc, k = ad.k)
rownames(ad.RUVIII) <- rownames(ad.clr)

#---------------------- Assessing batch effect correction ---------------------#

## Methods that detect batch effects
### PCA
ad.pca.before <- pca(ad.clr, ncomp = 3, scale = TRUE)
ad.pca.ComBat <- pca(ad.ComBat, ncomp = 3, scale = TRUE)
ad.pca.PLSDA_batch <- pca(ad.PLSDA_batch, ncomp = 3, scale = TRUE)
ad.pca.RUVIII <- pca(ad.RUVIII, ncomp = 3, scale = TRUE)

ad.batch = factor(ad.metadata$sequencing_run_date, 
                  levels = unique(ad.metadata$sequencing_run_date))

ad.pca.before.plot <- Scatter_Density(object = ad.pca.before, 
                                      batch = ad.batch, 
                                      trt = ad.trt, 
                                      title = 'Before correction')
ad.pca.ComBat.plot <- Scatter_Density(object = ad.pca.ComBat, 
                                      batch = ad.batch, 
                                      trt = ad.trt, 
                                      title = 'ComBat')
ad.pca.PLSDA_batch.plot <- Scatter_Density(object = ad.pca.PLSDA_batch, 
                                           batch = ad.batch, 
                                           trt = ad.trt, 
                                           title = 'PLSDA-batch')
ad.pca.RUVIII.plot <- Scatter_Density(object = ad.pca.RUVIII, 
                                      batch = ad.batch, 
                                      trt = ad.trt, 
                                      title = 'RUVIII')

grid.arrange(ad.pca.before.plot, ad.pca.ComBat.plot, 
             ad.pca.PLSDA_batch.plot, 
             ad.pca.RUVIII.plot, ncol = 2)

#### Exercise 4: Interpret the PCA plots created above

### pRDA
# arrange data before and after batch effect correction into a list
ad.corrected.list <- list(`Before correction` = ad.clr, 
                          ComBat = ad.ComBat, 
                          `PLSDA-batch` = ad.PLSDA_batch, 
                          RUVIII = ad.RUVIII)

ad.prop.df <- data.frame(Treatment = NA, Batch = NA, 
                         Intersection = NA, 
                         Residuals = NA) 

# run rda in a loop
for(i in seq_len(length(ad.corrected.list))){
  rda.res = varpart(ad.corrected.list[[i]], ~ trt, ~ batch,
                    data = ad.factors.df, scale = TRUE)
  ad.prop.df[i, ] <- rda.res$part$indfract$Adj.R.squared}

rownames(ad.prop.df) = names(ad.corrected.list)

# change the order of explained variance
ad.prop.df <- ad.prop.df[, c(1,3,2,4)]

# remove values less than zero, and recalculate the proportions
ad.prop.df[ad.prop.df < 0] = 0
ad.prop.df <- as.data.frame(t(apply(ad.prop.df, 1, 
                                    function(x){x/sum(x)})))

partVar_plot(prop.df = ad.prop.df)

########################################

## Other methods
### R2

# scale
ad.corr_scale.list <- lapply(ad.corrected.list, 
                             function(x){apply(x, 2, scale)})

# perform one-way ANOVA on each variable within each dataset
ad.r_values.list <- list()
for(i in seq_len(length(ad.corr_scale.list))){
  # for each dataset
  ad.r_values <- data.frame(trt = NA, batch = NA)
  for(c in seq_len(ncol(ad.corr_scale.list[[i]]))){
    # for each variable
    ad.fit.res.trt <- lm(ad.corr_scale.list[[i]][,c] ~ ad.trt)
    ad.r_values[c,1] <- summary(ad.fit.res.trt)$r.squared
    ad.fit.res.batch <- lm(ad.corr_scale.list[[i]][,c] ~ ad.batch)
    ad.r_values[c,2] <- summary(ad.fit.res.batch)$r.squared
  }
  ad.r_values.list[[i]] <- ad.r_values
}
names(ad.r_values.list) <- names(ad.corr_scale.list)
head(ad.r_values.list$`Before correction`)

# generate boxplots for each dataset
ad.boxp.list <- list()
for(i in seq_len(length(ad.r_values.list))){
  ad.boxp.list[[i]] <- 
    data.frame(r2 = c(ad.r_values.list[[i]][ ,'trt'],
                      ad.r_values.list[[i]][ ,'batch']), 
               Effects = as.factor(rep(c('Treatment','Batch'), 
                                       each = 231)))
}
names(ad.boxp.list) <- names(ad.r_values.list)
head(ad.boxp.list$`Before correction`)

ad.r2.boxp <- rbind(ad.boxp.list$`Before correction`,
                    ad.boxp.list$ComBat,
                    ad.boxp.list$`PLSDA-batch`,
                    ad.boxp.list$RUVIII)

ad.r2.boxp$methods <- rep(c('Before correction', 
                            'ComBat','PLSDA-batch',
                            'RUVIII'), each = 462)

ad.r2.boxp$methods <- factor(ad.r2.boxp$methods, 
                             levels = unique(ad.r2.boxp$methods))
head(ad.r2.boxp)

ggplot(ad.r2.boxp, aes(x = Effects, y = r2, fill = Effects)) +
  geom_boxplot(alpha = 0.80) +
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "right") + facet_grid( ~ methods) + 
  scale_fill_manual(values=pb_color(c(12,14))) 

### Alignment scores

# calculate the alignment scores
ad.scores <- c()
names(ad.batch) <- rownames(ad.clr)
for(i in seq_len(length(ad.corrected.list))){
  res <- alignment_score(data = ad.corrected.list[[i]], 
                         batch = ad.batch, 
                         var = 0.95, 
                         k = 8, 
                         ncomp = 50)
  ad.scores <- c(ad.scores, res)
}
head(ad.scores)

# rearrange the data for ggplot
ad.scores.df <- data.frame(scores = ad.scores, 
                           methods = names(ad.corrected.list))

ad.scores.df$methods <- factor(ad.scores.df$methods, 
                               levels = rev(names(ad.corrected.list)))

head(ad.scores.df)

ggplot() + geom_col(aes(x = ad.scores.df$methods, 
                        y = ad.scores.df$scores)) + 
  geom_text(aes(x = ad.scores.df$methods, 
                y = ad.scores.df$scores/2, 
                label = round(ad.scores.df$scores, 3)), 
            size = 3, col = 'white') + 
  coord_flip() + theme_bw() + ylab('Alignment Scores') + 
  xlab('') + ylim(0,0.85)

#------------------------------- To go further --------------------------------#
# see file "Batch_effect_management.Rmd"


