"R version-4.1.3"
# load packages and functions used in current analysis
suppressWarnings(suppressPackageStartupMessages({
    library(gwasvcf)
    library(TwoSampleMR)
    library(gassocplot)
    library(gwasglue)
}))

#
setwd('MM_tutorial')
source('Functions.R', echo = F, encoding = 'UTF-8')

# load previously save R data
load("MMMR.RData")

# Figure 1 ----
"subset data of genes LYZ, C1QB, BTN3A2"
subExpo <- c("eqtl-a-ENSG00000090382", "eqtl-a-ENSG00000173369", "eqtl-a-ENSG00000186470")
res_sub <- MR_res[MR_res$id.exposure %in% subExpo, ]
dat_sub <- dat2[dat2$id.exposure %in% subExpo, ]

# Figure 1A
pdf("MR_forestplot.pdf", width = 12, height = 15)
forest_mr(res_sub, wrap_str = 30, exponentiate = F, combine = T, rmID = T)
dev.off()

# Figure 1B
pdf("MR_Scatter.pdf", width = 15, height = 15)
MR_plot(dat_sub, res_sub, ncol = 3)
dev.off()

# Figure S3
pdf("MR_leaveoneout.pdf", width = 15, height = 15)
MR_plot(dat_sub, res_sub, ncol = 3, type = "leaveoneout_plot")
dev.off()


# Figure 2 ----

# genes LYZ, C1QB, BTN3A2
Genes <- c("LYZ", "C1QB", "BTN3A2")

pdf("Surv_OS_PFS.pdf", width = 9, height = 3)
panSurv("SKCM", Genes, ncol = 3, Surv = "OS", time = 365 * 50)

panSurv("SKCM", Genes, ncol = 3, Surv = "PFS", time = 365 * 50)
dev.off()

# Figure S1, S2
Genes <- c(
    "LYZ", "PKIA", "FPR2", "C1QB", "RUFY1", "BTN3A2",
    "DPM2", "KRT73", "PLVAP", "SERPINB8", "SLC38A11", "RTKN2"
)

pdf("GepiaSurv_OS_match.pdf", width = 9, height = 12)
panSurv("SKCM", Genes, ncol = 3, Surv = "OS", time = 365 * 50)
dev.off()


pdf("GepiaSurv_PFS.pdf", width = 9, height = 12)
panSurv("SKCM", Genes, ncol = 3, Surv = "PFS", time = 365 * 50)
dev.off()


# Figure 3 ----
LYZ <- read.csv("./CTDdb/LYZ_Chem.csv", header = T, row.names = 1)
C1QB <- read.csv("./CTDdb/C1QB_Chem.csv", header = T, row.names = 1)
BTN3A2 <- read.csv("./CTDdb/BTN3A2_Chem.csv", header = T, row.names = 1)

allGenes <- rbind(LYZ, C1QB, BTN3A2)
xlsx::write.xlsx(allGenes, "all_Chem.xlsx")

# Figure 3B
allGenes = xlsx::read.xlsx("all_Chem.xlsx")
allGenes$Effect <- unlist(lapply(allGenes$InteractionActions, function(x) strsplit(x, "\\|")[[1]][1]))
allGenes$Effect <- gsub("\\^", "\n", allGenes$Effect)

allHomo <- allGenes[grep("^Homo", allGenes$Organism), ]
dim(allHomo)
table(allHomo$Effect)

# format alluvial structure
sankeyDat <- allHomo[, c("GeneSymbol", "GeneForms", "Effect", "ChemicalID")]
sankeyDat <- Addper("GeneSymbol", sankeyDat, top = 3)
sankeyDat <- Addper("GeneForms", sankeyDat, top = 3)
sankeyDat <- Addper("ChemicalID", sankeyDat, top = 25, sepL = "")
sankeyDat <- Addper("Effect", sankeyDat, top = 7)

sankeyDat$Group <- sankeyDat$GeneSymbol
dat_lodes <- to_lodes_form(as.data.frame(sankeyDat), axes = 1:4)


pdf("allHomo_Alluvium.pdf", width = 15, height = 12)
ggplot(
    data = dat_lodes,
    aes(
        x = x, stratum = stratum, alluvium = alluvium,
        label = stratum
    )
) +
    geom_alluvium(aes(fill = Group),
        decreasing = NA, curve_type = "arctangent",
        width = 1 / 3, show.legend = T
    ) +
    geom_stratum() +
    geom_text(stat = "stratum", size = 3) +
    theme_minimal(base_size = 12) +
    labs(
        title = "Chemicals interaction-actions on gene/protein of LYZ, C1QB & BTN3A2",
        subtitle = "(Alluvial diagrams of chemicals effect on gene/protein)",
        caption = "Note: dataset acquired from CTDbase."
    ) +
    xlab("")
dev.off()

