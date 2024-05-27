library("WGCNA")
library(data.table)
library(dplyr)
#for tutorila  https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
data <- read.csv("/Volumes/USB/gene_network_feb_2024/cpm_collection/merged_final_pre_gene_network.csv",row.names=1)


# Check for missing values and remove them if necessary
datExpr <- na.omit(data)

# Choose a power for the soft-thresholding
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1, col="red")

# Choosing the power with the highest scale independence
softPower = sft$powerEstimate

# Creating the adjacency matrix
adjacency = adjacency(datExpr, power = softPower)

# Transform the adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Identifying modules via hierarchical clustering
geneTree = hclust(as.dist(dissTOM), method = "average")
# Dynamic Tree Cut to detect module
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE)

# Relating modules to external traits (if you have them)
# datTraits <- read.csv("traits.csv", row.names = 1)
# moduleTraitCor = cor(dynamicMods, datTraits, use = "p")

# Output the results
# write.csv(moduleTraitCor, file = "ModuleTraitCorrelations.csv")

## from anotehr source 

input_mat = t(datExpr)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 16
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
write.csv(module_df,
          file = "/Volumes/USB/gene_network_feb_2024/cpm_collection/gene_modules.txt",
)

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
# Get the names of the columns to pivot (excluding 'treatment')
#me_columns <- setdiff(names(MEs0), "treatment")
me_columns <- names(MEs0)



# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)
#  MEs0 is a data frame where columns are modules and rows are treatments
# Convert MEs0 from wide to long format
mME_base <- reshape2::melt(MEs0, id.vars = "treatment", variable.name = "module", value.name = "eigengene")

mME_base$module <- gsub("ME", "", mME_base$module)

module_order <- c("turquoise", "blue", "brown", "green", "yellow", "red", "grey", "black")

# Reassign the factor levels based on the preferred order
# Convert the 'module' column to a factor, setting the levels in the desired order
  # Make sure this is defined according to your module names
mME_base$module <- factor(mME_base$module, levels = module_order)



library(ggplot2)

# Replace 'module_column' with the actual name of the column that contains the module information
ggplot(mME_base, aes(x = treatment, y =module, fill = eigengene)) +
  geom_tile() +  # Use geom_tile to create the heatmap tiles
  theme_bw() +  # Use a black and white theme
  scale_fill_gradient2(  # Define a color gradient for the fill
    low = "blue",        # Color for low values
    high = "red",        # Color for high values
    mid = "white",       # Color for midpoint values
    midpoint = 0,        # The midpoint of the scale
    limits = c(-1, 1)    # The limits for the color scale
  ) +
  theme(axis.text.x = element_text(angle = 90)) +  # Rotate the x-axis text
  labs(
    title = "Module-trait Relationships",  # Add a title to the plot
    y = "Modules",                         # Label for the y-axis
    fill = "Correlation"                   # Label for the color legend
  )


### in the future analysis as normalizing the matrix using edgeR and then do the wgcna if have time 

# pick out a few modules of interest here
modules_of_interest = c("black", "turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]

subexpr = datExpr[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

## make gene_id columns for submod_df
subexpr$gene_id <- row.names(subexpr)
# Convert the subexpr data frame to long format
submod_df <- melt(data.frame(subexpr), id.vars = "gene_id")

# Rename the columns appropriately
colnames(submod_df) <- c("gene_id", "sample", "expression")

# Add module information based on gene_id
# Here we ensure that gene_id matches between submod_df and module_df
submod_df <- merge(submod_df, module_df, by.x = "gene_id", by.y = "row.names", all.x = TRUE)


submod_df %>% ggplot(., aes(x=sample, y=expression, group=gene_id)) +
  geom_line(aes(color = colors),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(colors)) +
  labs(x = "treatment",
       y = "normalized expression")
library(edgeR)

### using edge R normalization for above
dge <- DGEList(counts=datExpr, group=metadata_1$genotype, samples = metadata_1)
## filterByExpr function in edgeR is designed to filter out genes that are not expressed at a significant level across the samples in a dataset.
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
# Normalizing the Data
#Normalization is crucial for removing differences in library sizes and other technical variations.
dge <- calcNormFactors(dge)