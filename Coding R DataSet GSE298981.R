#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE298981", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13667", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "01010101"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Kontrol","Perlakuan"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

# calculate precision weights and show plot of mean-variance trend
v <- vooma(gset, design, plot=T)
# OR weights by group
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # attach gene annotations

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.Title","Gene.Symbol","GB_LIST","Entrez.Gene","Gene.Ontology.Biological.Process"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=2)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE298981", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE298981", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)


# Load library yang dibutuhkan
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")

library(ggplot2)
library(pheatmap)

# 1. Membuat data untuk Bar Chart (Analisis Jalur Fungsional)
pathway_data <- data.frame(
  Pathway = c("Epithelial Barrier & Keratin", "Cytoskeleton Remodeling", 
              "Immune Regulation (TL1A pathway)", "Acute Inflammation", 
              "Cellular Metabolism", "Protein Synthesis"),
  Count = c(8, 7, 6, 9, 5, 5), # Jumlah gen perkiraan
  Status = c("Down-regulated", "Up-regulated", "Up-regulated", 
             "Down-regulated", "Down-regulated", "Up-regulated")
)

# Membuat Bar Chart
ggplot(pathway_data, aes(x = reorder(Pathway, Count), y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Memutar grafik agar lebih mudah dibaca
  scale_fill_manual(values = c("Down-regulated" = "#2c7bb6", "Up-regulated" = "#d7191c")) +
  labs(title = "Functional Pathway Analysis (Top Affected)",
       subtitle = "Based on DEG Analysis & Gudino et al. (2025)",
       x = "Biological Processes",
       y = "Number of Genes") +
  theme_minimal()

# 1. Install dan Load Library
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# 2. Menyiapkan Data (Berdasarkan 40 gen teratas dari file Anda)
# memasukkan nilai Log2 Fold Change sebagai representasi ekspresi
genes <- c(
  # 20 Up-regulated (Merah)
  "ARHGEF37", "MLEC", "EIF5", "UNC119B", "SH3BP2", "KLF13", "EEF2K", "PDXK", 
  "TFCP2L1", "FOXK1", "PLCD3", "CBX6", "C1orf63", "RAB6B", "KHNYN", "ARHGAP1", 
  "CNIH4", "ZNF622", "VPS28", "TMEM115",
  # 20 Down-regulated (Biru)
  "KRT6A", "NAMPT", "HSPA8", "DUSP1", "SERPINA3", "SOCS3", "ALB", "SNORD14D", 
  "SC4MOL", "IER5L", "KCNS3", "PLA2G2A", "TNC", "CD24", "HBEGF", "C15orf48", 
  "CCL20", "LIF", "CXCL1", "IL1B"
)

# Nilai Log2 Fold Change (Data disederhanakan untuk visualisasi)
log2FC <- c(
  3.03, 2.61, 2.88, 3.18, 2.40, 2.12, 2.48, 2.49, 2.34, 2.37, 
  2.00, 2.62, 2.95, 2.33, 2.29, 2.04, 2.11, 2.15, 2.18, 2.20, # Up
  -5.26, -5.09, -4.28, -3.99, -3.71, -3.68, -3.67, -3.53, -3.53, -3.28, 
  -3.22, -3.08, -2.95, -2.88, -2.75, -2.66, -2.55, -2.48, -2.40, -2.35 # Down
)

# Membuat Matrix
heatmap_matrix <- matrix(log2FC, ncol = 1)
rownames(heatmap_matrix) <- genes
colnames(heatmap_matrix) <- "Log2 Fold Change"

# 3. Membuat Heatmap
pheatmap(heatmap_matrix,
         color = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(100),
         cluster_rows = TRUE,    # Mengelompokkan gen berdasarkan kemiripan
         cluster_cols = FALSE,   # Tidak perlu mengelompokkan kolom karena hanya ada satu
         main = "Heatmap Top 40 DEGs (PFD vs Control)",
         display_numbers = TRUE, # Menampilkan angka Log2FC di dalam kotak
         fontsize_row = 8,
         angle_col = 0)
