#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(tidyverse)
})

# ---- args ------------------------------------------------------------------
setwd("X:/BMG/watsonlab/OBX_capture/nf_flairr_demux")
input_csv  <- "demux_out/barcode_sharing_summary_05152026.csv"
output_pdf <- "demux_out/umi_sharing_heatmap_05152026.pdf"
mask_diag  <- T   # set TRUE to drop diagonal so cross-talk pops more

# ---- load ------------------------------------------------------------------
df <- fread("demux_out/barcode_sharing_summary_05152026.csv")
meta <- fread("../metadata/Kinnex_FLAIRR-seq_biosample_pool_meta_05_13_2026.csv")

# ---- filter to shared UMIs & binarize --------------------------------------
count_cols <- grep("^count_", names(df), value = TRUE)
df_cross   <- df[barcode_count > 1]

binary_matrix <- as.matrix(df_cross[, ..count_cols])
mode(binary_matrix) <- "integer"
binary_matrix <- (binary_matrix > 0) + 0L

# ---- co-occurrence ---------------------------------------------------------
co <- crossprod(binary_matrix)        # t(M) %*% M
diag(co) <- 0L

# drop samples with no cross-sample sharing
keep <- rowSums(co) > 0
co   <- co[keep, keep]

# strip count_ prefix
rownames(co) <- sub("^count_", "", rownames(co))
colnames(co) <- sub("^count_", "", colnames(co))

# ---- long format for ggplot ------------------------------------------------
long <- as.data.table(as.table(co))
setnames(long, c("sample_a", "sample_b", "shared"))

# preserve matrix order on axes
long[, sample_a := factor(sample_a, levels = rownames(co))]
long[, sample_b := factor(sample_b, levels = rev(rownames(co)))]  

long1 <- long %>% 
  left_join(meta %>% 
              select(sample_id, "pool_sampleA"="pool"), 
            by=c("sample_a"="sample_id")) %>% 
  left_join(meta %>% 
              select(sample_id, "pool_sampleB"="pool"), 
            by=c("sample_b"="sample_id")) %>% 
  mutate(same_pool=case_when(pool_sampleA == pool_sampleB ~ "T",
                             TRUE ~ "F"))

# ---- plot ------------------------------------------------------------------
p <- ggplot(long, 
            aes(sample_a, 
                sample_b, 
                fill = shared)) +
  geom_raster() +
  #facet_wrap(~same_pool) +
  #geom_text(aes(label = shared), size = 2) +
  scale_fill_gradient(low = "white", 
                      high = "#b2182b",
                      name = "Number of\nShared UMI Combos",
                      labels = comma) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid  = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Count of Shared UMI Combos Between Samples")

ggsave("pairwise_sharing_heatmap_annotated.pdf",
       p, width = 22, height = 18, dpi = 150)
