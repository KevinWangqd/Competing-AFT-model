library(GEOquery)

# 2. Read the series matrix file
gse <- getGEO(filename = "GSE72094_series_matrix.txt.gz")

# 3. Extract expression data (genes × samples)
expr <- exprs(gse)
expr
# 4. Extract sample metadata (clinical data)
pheno <- pData(gse)

# 5. Inspect
dim(expr)       # number of genes × samples
head(expr[,1:5])

# Extract relevant columns
surv_df <- pheno[, c(
  "survival_time_in_days:ch1", 
  "vital_status:ch1", 
  "age_at_diagnosis:ch1", 
  "gender:ch1", 
  "Stage:ch1", 
  "smoking_status:ch1", 
  "kras_status:ch1", 
  "stk11_status:ch1", 
  "tp53_status:ch1"
)]

# Rename columns for convenience
colnames(surv_df) <- c(
  "time_days", "status", "age", "gender", "stage",
  "smoking_status", "kras", "stk11", "tp53"
)

surv_df <- surv_df[which(surv_df$time_days!="NA"),]
surv_df$time_day <- as.numeric(surv_df$time_day)

# Recode vital_status: dead = 1 (event), alive = 0 (censored)
surv_df$status <- ifelse(tolower(surv_df$status) == "dead", 1, 0)
surv_df$gender

surv_df$stage_main <- substr(surv_df$stage, 1, 1)  # take first character
surv_df$stage_main <- factor(surv_df$stage_main, levels = c("1","2","3","4"))
write.csv(surv_df, "surv.csv")
# Inspect the first few rows
head(surv_df)





surv_data <- read.csv("final_data.csv") 
summary(surv_data)

library(GEOquery)



gpl <- getGEO("GPL15048", destdir = ".")
feature <- Table(gpl)
colnames(feature)

genes <- c(
  "CCNA2", "AURKA", "AURKB", "FEN1",
  "CD44", "CCND3", "NCALD", "MACF1",
  "MIR296", "RAMP2-AS1", "LRC4", "NLRC4",
  "PLEKHN1", "RASIP1", "SPP1",
  "GPT2", "FAM220A", "SGPL1", "PCOLCE2",
  "GABPB1-IT1"
)

# Extract probe IDs for each gene
probe_ids_list <- lapply(genes, function(gene) {
  feature$ID[grepl(gene, feature$GeneSymbol)]
})

names(probe_ids_list) <- genes
probe_ids_list

all_probes <- unlist(probe_ids_list)
all_probes <- c('merck-CR747584_s_at'	,
                'merck-DB326922_at'	,
                'merck-NM_004512_at',	
                'merck-NM_006277_s_at'	,
                'merck-NM_014065_s_at',	
                'merck-NM_019601_at'	,
                'merck-NM_024552_x_at',
                'merck-NM_198517_at',
                'merck2-NM_015947_at')
# Keep only probes that are present in the expression matrix
all_probes_available <- all_probes[all_probes %in% rownames(expr)]

# Subset the expression matrix
expr_subset <- expr[all_probes_available, , drop = FALSE]

# Check the result
dim(expr_subset)
head(expr_subset)

surv_data$ID

expr_t <- t(expr_subset)  # now rows = samples, cols = probes

# 2. Convert to data frame and add sample IDs as a column
expr_df <- as.data.frame(expr_t)
expr_df$ID <- rownames(expr_df)

# 3. Merge with survival data by sample ID
merged_data <- merge(surv_data, expr_df, by = "ID")

# 4. Check
dim(merged_data)
head(merged_data)
write.csv(merged_data,"final_data_new.csv")
