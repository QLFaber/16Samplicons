#####This is for quelccaya 16S V4 region

setwd("~/Ohio_work/quelccaya/qiime2_march26")

##feature-table contains counts
##taxonomy.tsv contains taxonomy
##if you want fasta file- you can upload the rep-seqs.tsv to view.qiime2.org and download fasta

##setup metadata table
metadata<-read.csv("quelccaya_metadata.csv",header=TRUE,row.names=1)

#count table
q2counts <- read.delim(
  "feature-table.tsv",
  header = TRUE,
  row.names = 1,
  skip = 1,
  check.names = FALSE,
  comment.char = ""
)

######prepare to rename everything
# Compute total abundance per ASV (row sums)
asv_totals <- rowSums(q2counts)

# Order ASVs from most to least abundant
q2counts <- q2counts[order(asv_totals, decreasing = TRUE), ]

# Number of ASVs (rows in q2counts)
num_asvs <- nrow(q2counts)


###import and set up taxonomy table
q2taxonomy<-read.delim("taxonomy.tsv",header=TRUE,row.names=1)

# Split taxonomy string into ranks
tax_split <- strsplit(as.character(q2taxonomy$Taxon), ";")

# Determine maximum number of ranks in your dataset
max_ranks <- max(sapply(tax_split, length))  # usually 7

# Pad shorter vectors with NA so all rows have same length
tax_pad <- lapply(tax_split, function(x) {
  length(x) <- max_ranks  # automatically pads with NA
  return(x)
})

# Combine into a matrix
tax_mat <- do.call(rbind, tax_pad)

# Remove rank prefixes (d__, p__, etc.)
tax_mat <- gsub("^[a-z]__", "", tax_mat)

# Replace empty strings with NA
tax_mat[tax_mat == ""] <- NA

# Assign DADA2-style column names
colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:max_ranks]

# Set rownames to ASV IDs
rownames(tax_mat) <- rownames(q2taxonomy)

# Final DADA2-style taxonomy table
taxa.q2 <- tax_mat

# Check
head(taxa.q2)
taxa.q2 <- taxa.q2[rownames(q2counts), ]


####skip if you dont need fasta file
# Step 1: Read FASTA file as a table
fasta_file <- "sequences.fasta"
lines <- readLines(fasta_file)

# Separate headers and sequences
headers <- lines[grep("^>", lines)]
seqs <- lines[!grepl("^>", lines)]

# Create a data frame
fasta_df <- data.frame(header = headers, sequence = seqs, stringsAsFactors = FALSE)

# Step 2: Make sure headers match q2counts rownames (remove > if necessary)
# e.g., if headers are like >ASV123
fasta_df$header <- sub("^>", "", fasta_df$header)

# Step 3: Keep only sequences present in q2counts
common_asvs <- intersect(fasta_df$header, rownames(q2counts))
fasta_df <- fasta_df[match(rownames(q2counts), fasta_df$header), , drop = FALSE]

# Step 4: Rename sequences to ASV_1, ASV_2, ...
fasta_df$header <- paste0("ASV_", 1:num_asvs)

# Step 5: Write renamed FASTA back
out_file <- "new_sequences_renamed.fasta"
writeLines(
  paste0(">", fasta_df$header, "\n", fasta_df$sequence),
  con = out_file
)
################



##########Actually rename 
rownames(q2counts) <- paste0("ASV_", 1:num_asvs)
write.csv(q2counts,"q2counts.csv")
rownames(taxa.q2) <- paste0("ASV_", 1:num_asvs)
write.csv(taxa.q2,"q2counts.csv")




###############remove contaminants-skip if you dont have a blank###########

library(decontam)
#Telling R where blanks are
#to list your samples use this
colnames(q2counts) 
# find your negative controls/blank columns. My blank is the 68th and 69th sample in my list
vector_for_decontam <- rep(FALSE, 135)
vector_for_decontam[68:69] <- TRUE
contam_df <- isContaminant(t(q2counts), neg=vector_for_decontam)
table(contam_df$contaminant) # identified 26 contaminants for me = "TRUE" column
contam_df
# getting IDs of those identified as likely contaminants
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
# taxonomy of your contaminants
taxa.q2[row.names(q2counts) %in% contam_asvs, ]


###keeping fasta file and taxonomy the same with contaminants,,but removing from counts file!
# making new count table
q2counts_nocontam <- q2counts[!row.names(q2counts) %in% contam_asvs, ]
write.csv(q2counts_nocontam, "ASVs_counts-no-contam.csv")



######################################################
