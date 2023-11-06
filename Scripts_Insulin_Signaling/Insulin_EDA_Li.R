library(edgeR)
library(here)
library(reshape)
library(RColorBrewer)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(wesanderson)
library(dplyr)

# This code was used to evaluate the activation of insulin and Igf signaling 
# in different populations at different time points in the BM
# data from Single-Cell Analysis of Neonatal HSC Ontogeny Reveals Gradual and 
# Uncoordinated Transcriptional Reprogramming that Begins before Birth were used

#### load the Li data and subset HSC ####
data_li <- read.delim(here("PublicData/Public_Data_Li/Li_GSE128759_gene_level_read_counts.txt"), 
                   row.names = "ensembl_gene_id")

dim(data_li)
colnames(data_li)

# split the data into HSC and HPC


# subset HSC E14, 16, and 18
data_hpc <- data_li[, c(1:21)]
colnames(data_hpc)

data_hsc <- data_li[, c(22:47)]
colnames(data_hsc)

# adapt column names to get sample information
# for HPC
col_names_hpc <- c()

col_tp_hpc <- c(rep("Adult", 4), rep("E16", 4), rep("P7", 4), rep("P0", 4),
                rep("P14", 5))
col_tp_hpc

col_rep_hpc <- c(1:4, 1:4, 1:4, 1:4, 1:5)
col_rep_hpc

for (i in 1:length(col_tp_hpc)) {
    col_names_hpc[i] <- paste0(col_tp_hpc[i], "_HPC_", col_rep_hpc[i], sep="")
}

col_names_hpc
colnames(data_hpc) <- col_names_hpc
head(data_hpc)

# for HSC
col_names_hsc <- c()

col_tp_hsc <- c(rep("E18", 3), rep("E16", 4), rep("E14", 3),
                rep("P7", 4), rep("P14", 4), rep("P0", 4), rep("Adult", 4))
col_tp_hsc

col_rep_hsc <- c(1:3, 1:4, 1:3, 1:4, 1:4, 1:4, 1:4)
col_rep_hsc

for (i in 1:length(col_tp_hsc)) {
  col_names_hsc[i] <- paste0(col_tp_hsc[i], "_HSC_", col_rep_hsc[i], sep="")
}

col_names_hsc
colnames(data_hsc) <- col_names_hsc
head(data_hsc)


#### preprocess data ####
# remove NAs
data_hsc_filt <- na.omit(data_hsc)
dim(data_hsc)
dim(data_hsc_filt)

data_hpc_filt <- na.omit(data_hpc)
dim(data_hpc)
dim(data_hpc_filt)


# check for duplicate rows
table(duplicated(rownames(data_hpc_filt)))
table(duplicated(rownames(data_hsc_filt)))


# remove outlier
# FL_E16_HSC_1
data_hsc_filt <- data_hsc_filt[,!colnames(data_hsc_filt) == "E16_HSC_1"]

# define experimental groups 
my_fun <- function(i) {
  remove_last_n <- 2
  substring(i, first = 1, nchar(i) - remove_last_n)
}

# for HSC
groups_hsc <- sapply(colnames(data_hsc_filt), my_fun)
groups_hsc <- unname(groups_hsc)
groups_hsc <- factor(groups_hsc)
groups_hsc

# for HPC
groups_hpc <- sapply(colnames(data_hpc_filt), my_fun)
groups_hpc <- unname(groups_hpc)
groups_hpc <- factor(groups_hpc)
groups_hpc


# build DGEL object and run EDA of only Li data#
dgel_hsc <- DGEList(counts = data_hsc_filt, group = groups_hsc)
dgel_hpc <- DGEList(counts = data_hpc_filt, group = groups_hpc)


#### filter counts ####
keep_hsc <- filterByExpr(dgel_hsc)
keep_hpc <- filterByExpr(dgel_hpc)


dgel_hsc <- dgel_hsc[keep_hsc,,keep.lib.sizes=FALSE]
dim(dgel_hsc$counts)

dgel_hpc <- dgel_hpc[keep_hpc,,keep.lib.sizes=FALSE]
dim(dgel_hpc$counts)

#### convert ENSEMBL to gene SYMBOL ####
# HSC
symbols_hsc <- select(org.Mm.eg.db, keys=rownames(dgel_hsc$counts), 
                  columns="SYMBOL", keytype="ENSEMBL")

# remove NAs and duplicates
symbols_hsc <- symbols_hsc[!is.na(symbols_hsc$SYMBOL),]

table(duplicated(symbols_hsc$ENSEMBL))
symbols_hsc <- symbols_hsc[!duplicated(symbols_hsc$ENSEMBL),]

# intersect with rownames in dgel and subset
overlap_hsc <- intersect(rownames(dgel_hsc$counts), symbols_hsc$ENSEMBL)
length(overlap_hsc)

dgel_hsc$counts <- dgel_hsc$counts[overlap_hsc,]

# rename rows in dgel
rownames(dgel_hsc$counts) <- symbols_hsc$SYMBOL
head(dgel_hsc$counts)

# HPC
symbols_hpc <- select(org.Mm.eg.db, keys=rownames(dgel_hpc$counts), 
                      columns="SYMBOL", keytype="ENSEMBL")

# remove NAs and duplicates
symbols_hpc <- symbols_hpc[!is.na(symbols_hpc$SYMBOL),]

table(duplicated(symbols_hpc$ENSEMBL))
symbols_hpc <- symbols_hpc[!duplicated(symbols_hpc$ENSEMBL),]

# intersect with rownames in dgel and subset
overlap_hpc <- intersect(rownames(dgel_hpc$counts), symbols_hpc$ENSEMBL)
length(overlap_hpc)

dgel_hpc$counts <- dgel_hpc$counts[overlap_hpc,]

# rename rows in dgel
rownames(dgel_hpc$counts) <- symbols_hpc$SYMBOL
head(dgel_hpc$counts)


#### normalise counts
dgel_hsc <- normLibSizes(dgel_hsc)
dgel_hpc <- normLibSizes(dgel_hpc)


design_hsc <- model.matrix(~groups_hsc)
dgel_hsc <- estimateDisp(dgel_hsc,design_hsc)

design_hpc <- model.matrix(~groups_hpc)
dgel_hpc <- estimateDisp(dgel_hpc,design_hpc)


#### run EDA
# HSC
pseudoCounts_hsc <- log2(dgel_hsc$counts + 1)
boxplot(pseudoCounts_hsc, col = "gray", las = 3, cex.names = 1)

colConditions_hsc <- wes_palette("BottleRocket1")
colConditions_hsc <- colConditions_hsc[match(groups_hsc,levels(groups_hsc))]
plotMDS(pseudoCounts_hsc, gene.selection="common", col = colConditions_hsc)

# HPC
pseudoCounts_hpc <- log2(dgel_hpc$counts + 1)
boxplot(pseudoCounts_hpc, col = "gray", las = 3, cex.names = 1)

colConditions_hpc <- wes_palette("Zissou1")
colConditions_hpc <- colConditions_hpc[match(groups_hpc,levels(groups_hpc))]
plotMDS(pseudoCounts_hpc, gene.selection="common", col = colConditions_hpc)

#### plot expression of insulin signaling genes
ins_genes <- c("Igf1", "Igf2", "Igf1r", "Igf2r", "Insr")

# HSC
# calculate CPM and plot expression
cpm_log_hsc <- cpm(dgel_hsc,log=TRUE)

group.dif_hsc<-cpm_log_hsc[rownames(cpm_log_hsc) %in% ins_genes,]
# edit group names
colnames(group.dif_hsc)<-dgel_hsc$samples$group
group.dif_hsc<-melt(group.dif_hsc)
names(group.dif_hsc)<-c("ID","group","cpm")
# reorder groups
groups_ord_hsc <- c("E14_HSC", "E16_HSC", "E18_HSC", "P0_HSC",
                    "P7_HSC", "P14_HSC", "Adult_HSC")
group.dif_hsc <-  group.dif_hsc %>%
  mutate(group =  factor(group, levels = groups_ord_hsc)) %>%
  arrange(group) 

# HSC
par(mfrow=c(3,2))
for(i in 1:length(ins_genes)){
  cpm_range <- group.dif_hsc[group.dif_hsc$ID==ins_genes[i],]$cpm
  boxplot(cpm~group,group.dif_hsc[group.dif_hsc$ID==ins_genes[i],],
          main=ins_genes[i],ylab="count (cpm log)",frame.plot=F,col=wes_palette("BottleRocket1"),
          las=2, xlab="", ylim=c(min(cpm_range),max(cpm_range)))
}

# HPC
# calculate CPM and plot expression
cpm_log_hpc <- cpm(dgel_hpc,log=TRUE)

group.dif_hpc<-cpm_log_hpc[rownames(cpm_log_hpc) %in% ins_genes,]
# edit group names
colnames(group.dif_hpc)<-dgel_hpc$samples$group
group.dif_hpc<-melt(group.dif_hpc)
names(group.dif_hpc)<-c("ID","group","cpm")
# reorder groups
groups_ord_hpc <- c("E16_HPC", "P0_HPC", "P7_HPC", "P14_HPC", "Adult_HPC")
group.dif_hpc <-  group.dif_hpc %>%
  mutate(group =  factor(group, levels = groups_ord_hpc)) %>%
  arrange(group) 

# HPC 
par(mfrow=c(3,2))
for(i in 1:length(ins_genes)){
  cpm_range <- group.dif_hpc[group.dif_hpc$ID==ins_genes[i],]$cpm
  boxplot(cpm~group,group.dif_hpc[group.dif_hpc$ID==ins_genes[i],],
          main=ins_genes[i],ylab="count (cpm log)",frame.plot=F,col=wes_palette("Zissou1"),
          las=2, xlab="", ylim=c(min(cpm_range),max(cpm_range)))
}
