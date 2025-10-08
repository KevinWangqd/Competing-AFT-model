library(GEOquery)

# 2. Read the series matrix file
gse <- getGEO(filename = "GSE76427_series_matrix.txt.gz")

# 3. Extract expression data (genes × samples)
expr <- exprs(gse)
expr
# 4. Extract sample metadata (clinical data)
pheno <- pData(gse)

# 5. Inspect
dim(expr)       # number of genes × samples
head(expr[,1:5])
colnames(pheno)



# Extract relevant columns
pheno_useful <- pheno[, c(
  "geo_accession",
  "patient id:ch1",
  "age (years):ch1",
  "gender (1=m, 2=f):ch1",
  "bclc_staging:ch1",
  "diagnosis:ch1",
  "tissue:ch1",
  "tnm_staging_clinical:ch1",
  "duryears_os:ch1",
  "event_os:ch1",
  "duryears_rfs:ch1",
  "event_rfs:ch1"
)]
head(pheno_useful)


colnames(pheno_useful) <- c(
  "GEO_ID", "Patient_ID", "Age", "Gender", "BCLC_Stage", "Diagnosis",
  "Tissue", "TNM_Stage", "OS_time", "OS_event", "RFS_time", "RFS_event"
)

# Keep only samples with OS_time not NA
pheno_filtered <- pheno_useful[!is.na(pheno_useful$OS_time), ]

head(pheno_filtered)



#CCDC107, CXCL12, GIGYF1, GMNN, and IFFO1)



library(GEOquery)

gpl <- getGEO("GPL10558", destdir = ".", AnnotGPL = TRUE)

# The object is an ExpressionSet-like list; extract the table
gpl_table <- Table(gpl)

# Check columns
head(gpl_table)
colnames(gpl_table)
gpl_table$`Gene symbol`

genes <- c("CCDC107", "CXCL12", "GIGYF1", "GMNN", "IFFO1")
probe_ids_list <- lapply(genes, function(g){
  gpl_table$ID[gpl_table$`Gene symbol` == g]
})
names(probe_ids_list) <- genes
probe_ids_list
probe_ids <- unlist(probe_ids_list)

# 5. Subset expression data (expr) for these probes
expr_subset <- expr[probe_ids, , drop = FALSE]

# 6. Optionally, rename rows to gene symbols if multiple probes exist
rownames(expr_subset) <- sapply(rownames(expr_subset), function(id) {
  g <- names(probe_ids_list)[sapply(probe_ids_list, function(x) id %in% x)]
  if(length(g) > 0) g else id
})

# 7. Merge expression data with phenotype
# Ensure that sample IDs in expression (colnames) match pheno$ID
common_samples <- intersect(colnames(expr_subset), pheno_filtered$GEO_ID)

expr_subset <- expr_subset[, common_samples, drop = FALSE]
expr_subset <- t(expr_subset)

pheno_matched <- pheno_filtered[match(rownames(expr_subset), pheno_filtered$GEO_ID), ]

# Combine into one data frame
combined_data <- cbind(pheno_matched, expr_subset)

# Check
head(combined_data)

write.csv(combined_data, "surv_GSE76427.csv")
