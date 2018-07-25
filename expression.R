## ----message=FALSE-------------------------------------------------------
# source("requirements.R")
source("functions.R")
library(magrittr)
library(dplyr)
library(survival)
library(GGally)

## ------------------------------------------------------------------------
es <- read.gct("data/lgg.gct")

## ------------------------------------------------------------------------
es

## ------------------------------------------------------------------------
es <- es[, es$sample_type == "Primary solid Tumor"]

## ------------------------------------------------------------------------
exprs(es)[is.na(exprs(es)) | exprs(es) < 0] <- 0

## ------------------------------------------------------------------------
fData(es)$mean_expression <- apply(exprs(es), 1, mean)
es <- es[head(order(fData(es)$mean_expression, decreasing=TRUE), 10000)]


## ----ggplot2-1, fig.height=3.5, fig.width=7,dev='svg', message=FALSE-----
library(ggplot2)
ggplot(data=mtcars, aes(x=mpg, y=hp, color = as.factor(gear))) + 
    geom_point()

## ----ggplot2-2, fig.height=3.5, fig.width=7,dev='svg'--------------------
library(ggplot2)
ggplot(data=mtcars, aes(x=mpg, y=hp, color = as.factor(gear))) + 
    geom_point() +
    theme_bw()

## ----ggplot2-3, fig.height=3.5, fig.width=7,dev='svg'--------------------
library(ggplot2)
ggplot(data=mtcars, aes(x=mpg, y=hp, color = as.factor(gear))) + 
    geom_point() +
    theme_bw() +
    ggtitle("mtcars: mpg vs hp")

## ----pca, cache=TRUE, cache.path="./expression-cache/"-------------------
pca <- prcomp(t(exprs(es)))
str(pca, list.len=5)

## ------------------------------------------------------------------------
str(pca$sdev)
pca_explained <- (pca$sdev)^2/sum(pca$sdev^2)*100
head(pca_explained)

## ----fig.height=3, fig.width=7,dev='svg'---------------------------------
pca_table <- cbind(as.data.frame(pca$x), 
                   pData(es))
pp <- ggplot(data = pca_table) + 
    geom_point(mapping=aes(x=PC1, y=PC2, color=histological_type)) +
    xlab(sprintf("PC1 (%.1f%%)", pca_explained[1])) +
    ylab(sprintf("PC2 (%.1f%%)", pca_explained[2]))
pp

## ----fig.height=3, fig.width=7,dev='svg'---------------------------------
ggsave(pp, file="./pca12.png", width=6, height=4)

## ------------------------------------------------------------------------
es_sub <- es[, es$histological_type %in% c("astrocytoma", "oligodendroglioma")]

## ---- message=FALSE------------------------------------------------------
library(limma)

## ------------------------------------------------------------------------
es.design <- model.matrix(~0+histological_type, data=pData(es_sub))

head(es.design)

## ------------------------------------------------------------------------
fit <- lmFit(es_sub, es.design)

fit2 <- contrasts.fit(fit, makeContrasts(histological_typeoligodendroglioma-
                                             histological_typeastrocytoma, 
                                         levels=es.design))
fit2 <- eBayes(fit2)

## ------------------------------------------------------------------------
de <- topTable(fit2, adjust.method="BH", number=Inf)
head(de)

## ----message=FALSE-------------------------------------------------------
library(AnnotationDbi)
library(org.Hs.eg.db)

de$entrez <- mapIds(org.Hs.eg.db, keys = de$id,
                    keytype = "SYMBOL", column = "ENTREZID")
head(de)

## ------------------------------------------------------------------------
de_sub <- de %>%
    arrange(desc(AveExpr)) %>%
    filter(!duplicated(entrez) & !is.na(entrez))

stats <- setNames(de_sub$t, de_sub$entrez)

str(stats)

## ---- message=FALSE------------------------------------------------------
library(fgsea)
library(data.table)

pathways <- reactomePathways(names(stats))
str(pathways, list.len=5)

## ----fgsea,cache=TRUE, cache.path="expression-cache/"--------------------
fr <- fgsea(pathways, stats = stats, nperm=10000, minSize=15, maxSize=500)
head(fr)

## ------------------------------------------------------------------------
frSig <- fr[padj < 0.01][order(NES)]
head(frSig)

## ------------------------------------------------------------------------
tail(frSig)

## ----gseaPlot, fig.height=3.5, fig.width=7,dev='svg'---------------------
plotEnrichment(pathways$`Cholesterol biosynthesis`, stats) +
    ggtitle("Cholesterol biosynthesis")

## ----pheatmap, fig.height=3.5, fig.width=7,dev='svg'---------------------
library(pheatmap)
m <- matrix(rnorm(200), 20, 10)
pheatmap(m)

## ------------------------------------------------------------------------
pathwayEntrez <- pathways$`Cholesterol biosynthesis`
pathwaySymbol <- mapIds(org.Hs.eg.db, pathwayEntrez,
                        column="SYMBOL", keytype="ENTREZID")


es_cholesterol <- es_sub[fData(es_sub)$id %in% pathwaySymbol, ]

heatmap_table <- exprs(es_cholesterol)
heatmap_table <- t(apply(heatmap_table, 1, scales::rescale))

## ----cholesterol-1, fig.height=3.5, fig.width=7,dev='svg'----------------

pheatmap(heatmap_table, 
         show_colnames = FALSE)

## ------------------------------------------------------------------------
es_cholesterol <- es_cholesterol[, order(apply(exprs(es_cholesterol), 2, mean))]

heatmap_table <- exprs(es_cholesterol)
heatmap_table <- t(apply(heatmap_table, 1, scales::rescale))
rownames(heatmap_table) <- fData(es_cholesterol)$id

annotation_col <- data.frame(row.names=colnames(es_cholesterol),
                             histological_type=es_cholesterol$histological_type,
                             stringsAsFactors = FALSE)



## ----cholesterol-2, fig.height=3.5, fig.width=10,dev='png'---------------
pheatmap(heatmap_table, 
         show_colnames = FALSE, cluster_cols = FALSE,
         annotation_col = annotation_col)

## ------------------------------------------------------------------------
es$pathway_score <- apply(exprs(es)[fData(es)$id %in% pathwaySymbol, ],
                          2, mean)
load("clinical.rda")
clinical$pathway_score <- 
    es$pathway_score[match(clinical$bcr_patient_barcode,
                           es$bcr_patient_barcode)]


## ----cholesterol-surv-1, fig.height=3.5, fig.width=7,dev='png'-----------
median(clinical$pathway_score)
clinical %>% 
    survfit(Surv(years_to_event, event_status) ~ pathway_score > 10.08,
            data=.) %>%
    ggsurv

## ----cholesterol-surv-2, fig.height=3.5, fig.width=7,dev='png'-----------
clinical %>% 
    filter(histological_type == "astrocytoma") %>%
    survfit(Surv(years_to_event, event_status) ~ pathway_score > 10.08, 
            data=.) %>%
    ggsurv

## ----cholesterol-surv-3, fig.height=3.5, fig.width=7,dev='png'-----------
clinical %>% 
    filter(histological_type == "oligodendroglioma") %>%
    survfit(Surv(years_to_event, event_status) ~ pathway_score > 10.08,
            data=.) %>%
    ggsurv

## ------------------------------------------------------------------------
clinical %>% 
    filter(histological_type == "oligodendroglioma") %>%
    survdiff(Surv(years_to_event, event_status) ~ pathway_score > 10.08, data=.)

## ------------------------------------------------------------------------
clinical %>% 
    filter(histological_type == "oligodendroglioma") %>%
    coxph(Surv(years_to_event, event_status) ~ pathway_score, data=.)


