#Requiered libraries
library(FSA)
library(dplyr)
library(car)
library(ggplot2)
library(devtools)
library(pairwiseAdonis)
library(vegan)
library(openxlsx)
library(BiocManager)
library(Biostrings)
library(tidyr)
library(stringr)
library(patchwork)
library(VennDiagram)
library(corrplot)
library(ggsignif)


#################################################################################
####16S rRNA-Diversity study#####################################################
#################################################################################


#RAREFACTION CURVE###########################################################
# Import ASV-table from DADA2 pipeline. Dataset structure should be: Rows = samples, Columns = ASVs and remove Singletons, doubletons, tripletons, and quadrupletons 
raw_ASV_table <- read.csv2("data/ASV_Table_DADA2.csv")
raw_ASV_table_noSampleNames<-raw_ASV_table[2:3196] #required to make code commands working
ASV_table_noSampleNames<-raw_ASV_table_noSampleNames[, colSums(raw_ASV_table_noSampleNames, na.rm = TRUE) > 4, drop = FALSE]


#ASV_table_noSampleNames<-ASV_table[2:3196] #required to make code commands working
SampleNames<-c("A1_AS","A2_AS","A3_AS","A4_AS","A5_AS",
               "B1_AS","B2_AS","B3_AS","B4_AS","B5_AS",
               "C1_AS","C2_AS","C3_AS","C4_AS","C5_AS",
               "D1_AS","D2_AS","D3_AS","D4_AS","D5_AS",
               "E1_AS","E2_AS","E3_AS","E4_AS","E5_AS",
               "A1_TWW","A2_TWW","A3_TWW","A4_TWW","A5_TWW",
               "B1_TWW","B2_TWW","B3_TWW","B4_TWW","B5_TWW",
               "C1_TWW","C2_TWW","C3_TWW","C4_TWW","C5_TWW",
               "D1_TWW","D2_TWW","D3_TWW","D4_TWW","D5_TWW",
               "E1_TWW","E2_TWW","E3_TWW","E4_TWW","E5_TWW")

ASV_table<-ASV_table_noSampleNames
ASV_table$Sample<-SampleNames
write.xlsx(ASV_table,"ASV_table.xlsx")

#Rarefaction curves

# Create the rarefaction curve
tiff("Rarefaction_curve.tiff", units="cm", width=30, height=25, res=300)
curve<-rarecurve(ASV_table_noSampleNames,
                 ylab = "ASV",
                 step = 50,         # step size for sample size increments
                 col = rainbow(nrow(ASV_table_noSampleNames)),  # color per sample
                 cex = 0.8,        # point size
                 label = TRUE,   # add sample labels
)

# Add legend (all samples)
legend("topright", 
       legend = SampleNames,   # sample names
       col = rainbow(nrow(ASV_table_noSampleNames)), 
       lty = 1, 
       cex = 0.6, 
       ncol = 2)   # adjust number of columns to fit better

abline(v = 2588, col = "black", lty = 2, lwd = 2)

dev.off()

#Create and export rarefied ASV-table

set.seed(123)  # for reproducibility
rarefied_2588 <- rrarefy(ASV_table_noSampleNames, sample = 2588)
rarefied_2588 <- as.data.frame(rarefied_2588)
rarefied_2588$Sample<-SampleNames

write.xlsx(rarefied_2588,"Rarefied_ASV_table_2588.xlsx")

#SELECT LEGIONELLA ASVs####################################################

#Import custom Legionella database
Leg_DB <- readDNAStringSet("./Data/legionella_database.fasta")

#Relate BLASTn-tool results with ASVs from DADA2
#Create a table with three columns: Sequencing ID, Genus, Species
headers <- names(Leg_DB)
parsed_data <- data.frame(
  Sequence_ID = sapply(strsplit(headers, " "), `[`, 1),
  Genus=sapply(strsplit(headers," "),`[`,2),
  Species = sapply(strsplit(headers, " "), `[`, 3),
  stringsAsFactors = FALSE
)

#import BLASTn-tool output
blast_out <- read.csv2("data/BLASTn_output.csv")

#merge parsed_data and blast_out by sequence_ID
mergeddata1 <- left_join(blast_out, parsed_data, by = "Sequence_ID")

#import file that related ASV ID with DADA2 ID
dada2_asvid<-read.csv2("data/dada2id_asvid.csv")

#merge mergedata1 with dada2_asvid
mergeddata2 <- left_join(dada2_asvid, mergeddata1, by = "DADA2ID")

#Select in mergeddata2 those rows with highest percentage identity for each ASV
Highest_percentage_ID <- mergeddata2 %>%
  group_by(ASV_ID) %>%
  filter(pident_corr == max(pident_corr)) %>%
  ungroup()

#Select in Highest_percentage_ID those rows percentage identity of min 93%
Highest_percentage_ID_above_93 <- Highest_percentage_ID %>%
  filter(pident_corr >93)

#Export table containing ASVs with highest percentage identity of min. 93%
write.xlsx(Highest_percentage_ID_above_93,"Highest_percentage_ID_abover_93.xlsx")

#Create a rarefied matrix containing only ASVs corresponding to Legionella
legionella_ASVs<-as.character(Highest_percentage_ID_above_93$ASV_ID)

legionella_matrix<-rarefied_2588[,colnames(rarefied_2588) %in% legionella_ASVs ]

#expoert Legionella-ASV matrix
write.xlsx(legionella_matrix,"Legionella_matrix.xlsx")


#ALPHA-DIVERSITY###########################################################



#######################################################################

# Calculate alpha diversity
alpha_div <- data.frame(
  Sample   = SampleNames,
  Observed = specnumber(legionella_matrix),
  Shannon  = diversity(legionella_matrix, index = "shannon")
)

alpha_div <- alpha_div %>%
  mutate(
    Origin = str_extract(Sample, "^[A-Za-z]+"),
    Type   = str_extract(Sample, "(?<=_)\\w+$"),
    Group  = paste(Origin, Type, sep = "_")
  )

# Convert to long format
alpha_long <- alpha_div %>%
  pivot_longer(
    cols = c(Observed, Shannon),
    names_to = "Index",
    values_to = "Value"
  )

# Kruskal-Wallis tests: AS vs TWW within each WWTP and diversity index
kruskal_AS_TWW <- alpha_long %>%
  group_by(Origin, Index) %>%
  summarise(
    p_value = kruskal.test(Value ~ Type)$p.value,
    .groups = "drop"
  )

print(kruskal_AS_TWW)

# Kruskal-Wallis tests: WWTP comparison separately for AS and TWW
kruskal_WWTP <- alpha_long %>%
  group_by(Type, Index) %>%
  summarise(
    p_value = kruskal.test(Value ~ Origin)$p.value,
    .groups = "drop"
  )

print(kruskal_WWTP)

# Dunn test function for pairwise WWTP comparisons within AS or TWW
get_sig_comparisons <- function(data, sample_type){
  
  data_sub <- data %>% filter(Type == sample_type)
  
  dunn <- dunnTest(
    Value ~ Origin,
    data = data_sub,
    method = "none"
  )$res
  
  dunn_sig <- dunn %>%
    filter(P.unadj < 0.05)
  
  if(nrow(dunn_sig) == 0){
    return(NULL)
  }
  
  comps <- strsplit(as.character(dunn_sig$Comparison), " - ")
  ymax <- max(data_sub$Value, na.rm = TRUE)
  
  data.frame(
    xmin = paste0(sapply(comps, '[', 1), "_", sample_type),
    xmax = paste0(sapply(comps, '[', 2), "_", sample_type),
    y_position = seq(
      ymax * 1.10,
      ymax * (1.10 + 0.12 * (nrow(dunn_sig) - 1)),
      length.out = nrow(dunn_sig)
    ),
    annotation = ifelse(
      dunn_sig$P.unadj < 0.001, "***",
      ifelse(dunn_sig$P.unadj < 0.01, "**", "*")
    )
  )
}

# Plot alpha diversity
plot_list <- alpha_long %>%
  group_split(Index) %>%
  lapply(function(df){
    
    current_index <- unique(df$Index)
    
    ylab <- ifelse(
      current_index == "Observed",
      "Observed richness",
      "Shannon index"
    )
    
    sig_AS  <- get_sig_comparisons(df, "AS")
    sig_TWW <- get_sig_comparisons(df, "TWW")
    
    # Manual significance bracket for D_AS vs D_TWW
    # Change annotation if needed, e.g. "**" or "***"
    sig_D_AS_TWW <- if (unique(df$Index) == "Observed") {
      data.frame(
        xmin = "D_AS",
        xmax = "D_TWW",
        y_position = max(df$Value, na.rm = TRUE) * 1.75,
        annotation = "*"
      )
    } else {
      NULL
    }
    
    p <- ggplot(df, aes(x = Group, y = Value, fill = Type)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 1) +
      labs(
        title = NULL,
        x = "WWTP",
        y = ylab
      ) +
      theme_minimal(base_size = 18) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
      ) +
      annotate(
        "text",
        x = c(1.5, 3.5, 5.5, 7.5, 9.5),
        y = -Inf,
        label = c("A", "B", "C", "D", "E"),
        vjust = 2,
        size = 5
      )
    
    if(!is.null(sig_AS)){
      p <- p +
        geom_signif(
          data = sig_AS,
          aes(
            xmin = xmin,
            xmax = xmax,
            annotations = annotation,
            y_position = y_position
          ),
          manual = TRUE,
          inherit.aes = FALSE,
          textsize = 5,
          vjust = 0.75
        )
    }
    
    if(!is.null(sig_TWW)){
      sig_TWW$y_position <- sig_TWW$y_position + max(df$Value, na.rm = TRUE) * 0.25
      
      p <- p +
        geom_signif(
          data = sig_TWW,
          aes(
            xmin = xmin,
            xmax = xmax,
            annotations = annotation,
            y_position = y_position
          ),
          manual = TRUE,
          inherit.aes = FALSE,
          textsize = 5,
          vjust = 0.75
        )
    }
    
    if(!is.null(sig_D_AS_TWW)){
      p <- p +
        geom_signif(
          data = sig_D_AS_TWW,
          aes(
            xmin = xmin,
            xmax = xmax,
            annotations = annotation,
            y_position = y_position
          ),
          manual = TRUE,
          inherit.aes = FALSE,
          textsize = 5,
          vjust = 0.75
        )
    }
    
    return(p)
  })

# Save plot
tiff("Alpha_Legionella.tiff",
     units = "cm", width = 30, height = 25, res = 300)

wrap_plots(plot_list, ncol = 2) %>% print()

dev.off()

#Correlation of ALPHA-DIVERSITY####


mean <- alpha_long %>%
  group_by(Origin, Type, Index) %>%
  summarize(MedianValue = mean(Value, na.rm = TRUE))

Val_corr<-read.csv2("data/For_alpha_corr_avg.csv")
num_data <- Val_corr[, -1]

# Compute the Spearman correlation matrix
cor_matrix <- cor(num_data, method = "spearman", use = "pairwise.complete.obs")

# Custom function for significance testing with Spearman correlations
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(mat[, i], mat[, j], method = "spearman", ...)
      p.mat[i, j] <- p.mat[j, i] <- test$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(list(p = p.mat))
}

# Perform a correlation test to get p-values
cor_test <- cor.mtest(num_data)

# Create a logical matrix for significant correlations (e.g., p < 0.05)
significant <- cor_test$p < 0.05

# Define the new variable names
new_var_names <- c("Observed AS", "Observed TWW", "Shannon AS", "Shannon TWW","qPCR InAST","qPCR AS","qPCR TWW", "#Aerated processes", "Temperature",
                   "Infuent COD", "Influent TN","Influent TP", "Influent ammonium", "Influent nitrate","Influent nitrite","Influent Ortho-phosphate","Influent protein",
                   "Influent iron(III)", "Influent iron(II)", "Influent iron(II)/(III)")

# Assign the new names to the columns and rows of the correlation matrix

colnames(cor_matrix) <- new_var_names
rownames(cor_matrix) <- new_var_names
tiff("Correlation.tiff", units="cm", width=30, height=25, res=300)
corrplot(cor_matrix, method = "circle", type = "upper", 
         order = "original", p.mat = cor_test$p, sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 1.2, cl.cex=1.2,mar=c(0,0,0,0)
)
mtext(expression("Spearman's rank correlation coefficient"),
      side = 4, line = -1, las = 3, cex=1.5, at=7.5)
dev.off()
#BETA-DIVERSITY###########################################################

#Prepare data, all zero-rows must be removed
legionella_matrix_no0 <- legionella_matrix[rowSums(legionella_matrix) != 0, ]

# Calculate Bray-Curtis distance matrix
bray_curtis_dist <- vegdist(legionella_matrix_no0, method = "bray")

#Dendrogram
# Perform hierarchical clustering
hc <- hclust(bray_curtis_dist, method = "average")

# Plot the dendrogram
legionella_matrix_withNames<-legionella_matrix
legionella_matrix_withNames$Sample<-SampleNames
legionella_matrix_withNames_no0 <- legionella_matrix_withNames[rowSums(legionella_matrix) != 0, ]


labels_vec <- legionella_matrix_withNames_no0$Sample

# Sanity check: number of labels must equal number of rows in asv_table_no0
if(length(labels_vec) != nrow(legionella_matrix_no0)){
  stop("Length mismatch: length(labels_vec) = ", length(labels_vec),
       " but nrow(legionella_matrix_no0) = ", nrow(legionella_matrix_no0),
       ". Make sure labels_vec is in the same order and corresponds to the rows used to compute bray_curtis_dist.")
}

# Ensure the rownames of the data used for distance match your labels (this keeps everything consistent)
rownames(legionella_matrix_no0) <- labels_vec

# Clustering based on bray-curtis distanece
hc <- hclust(bray_curtis_dist, method = "average")

# Force hclust to carry the proper labels (important if anything else changed them)
hc$labels <- labels_vec

# Build color mapping by first letter
first_letter <- substr(labels_vec, 1, 1)
unique_letters <- sort(unique(first_letter))   # sort for reproducibility
# Use a palette; change to your own if you prefer
palette_colors <- palette_colors <- c(
  A = "orange",
  B = "lightblue",
  C = "lightgreen",
  D = "gold",
  E = "darkblue"
)
label_colors <- palette_colors[first_letter]

# Function to color leaves (edge and label)
color_branches <- function(x) {
  if (is.leaf(x)) {
    lab <- attr(x, "label")
    this_letter <- substr(lab, 1, 1)
    col <- palette_colors[this_letter]
    # color the edge leading to this leaf
    attr(x, "edgePar") <- c(attr(x, "edgePar"), list(col = col, lwd = 2))
    # color the label
    attr(x, "nodePar") <- c(attr(x, "nodePar"), list(lab.col = col, pch = NA))
  }
  x
}

# Convert to dendrogram and apply coloring
dend <- as.dendrogram(hc)
dend_colored <- dendrapply(dend, color_branches)

# Plot to file
tiff("Beta_dendrogram_Legionella.tiff", units="cm", width=30, height=25, res=300)
plot(dend_colored, main = "", xlab = "Sample", ylab = "Bray–Curtis distance",
     cex = 1.2, cex.axis = 1.2, cex.lab = 1.5)

dev.off()

#Significance test, PERMANOVA
#comparing sample groups between each other
group_factor1 <- c("A_AS", "A_AS", "A_AS","A_AS","A_AS",
                   "B_AS","B_AS","B_AS","B_AS","B_AS",
                   "C_AS","C_AS","C_AS","C_AS","C_AS",
                   "D_AS","D_AS","D_AS","D_AS","D_AS",
                   "E_AS","E_AS","E_AS","E_AS","E_AS",
                   "A_TWW","A_TWW","A_TWW","A_TWW","A_TWW",
                   "B_TWW","B_TWW","B_TWW","B_TWW","B_TWW",
                   "C_TWW","C_TWW","C_TWW","C_TWW","C_TWW",
                   "D_TWW","D_TWW","D_TWW","D_TWW","D_TWW",
                   "E_TWW","E_TWW","E_TWW","E_TWW","E_TWW")

permanova_result <- adonis2(bray_curtis_dist ~ group_factor1)
permanova_pairwise_results <- pairwise.adonis(bray_curtis_dist,factors=group_factor1,p.adjust.m = "fdr")
permanova_pairwise_results_nonsig1 <-permanova_pairwise_results%>%filter(p.adjusted>0.05) 

#comparing WWTPSs between eachother
group_factor2 <- c("A", "A", "A","A","A",
                   "B","B","B","B","B",
                   "C","C","C","C","C",
                   "D","D","D","D","D",
                   "E","E","E","E","E",
                   "A","A","A","A","A",
                   "B","B","B","B","B",
                   "C","C","C","C","C",
                   "D","D","D","D","D",
                   "E","E","E","E","E")

permanova_result_wwtp <- adonis2(bray_curtis_dist ~ group_factor2)
permanova_pairwise_results_wwtp <- pairwise.adonis(bray_curtis_dist,factors=group_factor2,p.adjust.m = "fdr")
permanova_pairwise_results_nonsig_wwtp <-permanova_pairwise_results_wwtp%>%filter(p.adjusted>0.05) 

#ZETA-DIVERSITY's Venn-diagrams ###########################################################

#Transform legionella_matrix into binary table (1 = present, 0= not present)
legionella_matrix_bi<-as.data.frame(lapply(legionella_matrix, function(x) as.integer(x > 0)))

legionella_matrix_bi_withNames<-legionella_matrix_bi
legionella_matrix_bi_withNames$Sample<-SampleNames
write.xlsx(legionella_matrix_bi_withNames,"legionella_matrix_bi.xlsx")

#Create binary matrix combining AS and TWW (1 =present in AS and/or TWW, 0 = not present in AS nor TWW)
legionella_matrix_bi_combined <- matrix(0, nrow = 25, ncol = ncol(legionella_matrix))
colnames(legionella_matrix_bi_combined) <- colnames(legionella_matrix)
for (i in 1:25) {
  # Check if there is any value > 0 in either row i or row i+24 in the same column
  legionella_matrix_bi_combined[i, ] <- as.numeric(
    apply(legionella_matrix[c(i, i + 25), ], 2, function(x) any(x > 0))
  )
}

SampleNames_combined<-c("A1", "A2", "A3","A4","A5",
                        "B1","B2","B3","B4","B5",
                        "C1","C2","C3","C4","C5",
                        "D1","D2","D3","D4","D5",
                        "E1","E2","E3","E4","E5")


legionella_matrix_bi_combined_withNames <- cbind(legionella_matrix_bi_combined, Sample = SampleNames_combined)
write.xlsx(legionella_matrix_bi_combined_withNames,"legionella_matrix_bi_combined.xlsx")

#WWTP A####
#AS+TWW
tiff("Venn_A_AS_TWW.tiff", units="cm", width=30, height=25, res=300)
select_rows<-c(1:5)
df<-legionella_matrix_bi_combined[select_rows,]

venn_list <- apply(df, 1, function(x) {
  names(x)[x == 1]
},simplify=FALSE)

venn_list <- as.list(venn_list)

# Name each list element by sample name (rownames of df)
names(venn_list) <- SampleNames_combined[select_rows]

#Venn-plot
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  alpha = 0.5,
  margin=0.1,
  cex = 1.5,
  cat.cex = 1.2,
  cat.dist =0.3,
)
grid::grid.draw(venn.plot)

dev.off()

#WWTP B####
#AS+TWW
tiff("Venn_B_AS_TWW.tiff", units="cm", width=30, height=25, res=300)
select_rows<-c(6:10)
df<-legionella_matrix_bi_combined[select_rows,]

venn_list <- apply(df, 1, function(x) {
  names(x)[x == 1]
},simplify=FALSE)

venn_list <- as.list(venn_list)

# Name each list element by sample name (rownames of df)
names(venn_list) <- SampleNames_combined[select_rows]

#Venn-plot
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  alpha = 0.5,
  margin=0.1,
  cex = 1.5,
  cat.cex = 1.2,
  cat.dist =0.3,
)
grid::grid.draw(venn.plot)

dev.off()

#WWTP C####

#AS+TWW
tiff("Venn_C_AS_TWW.tiff", units="cm", width=30, height=25, res=300)
select_rows<-c(11:15)
df<-legionella_matrix_bi_combined[select_rows,]

venn_list <- apply(df, 1, function(x) {
  names(x)[x == 1]
},simplify=FALSE)

# Apply() returned a matrix-like object, convert to list
venn_list <- as.list(venn_list)

# Name each list element by sample name (rownames of df)
names(venn_list) <- SampleNames_combined[select_rows]

#Venn-plot
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  alpha = 0.5,
  margin=0.1,
  cex = 1.5,
  cat.cex = 1.2,
  cat.dist =0.3,
)
grid::grid.draw(venn.plot)

dev.off()

#WWTP D####

#AS+TWW
tiff("Venn_D_AS_TWW.tiff", units="cm", width=30, height=25, res=300)
select_rows<-c(16:20)
df<-legionella_matrix_bi_combined[select_rows,]

venn_list <- apply(df, 1, function(x) {
  names(x)[x == 1]
},simplify=FALSE)

venn_list <- as.list(venn_list)

# Name each list element by sample name (rownames of df)
names(venn_list) <- SampleNames_combined[select_rows]

#Venn-plot
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  alpha = 0.5,
  margin=0.1,
  cex = 1.5,
  cat.cex = 1.2,
  cat.dist =0.3,
)
grid::grid.draw(venn.plot)

dev.off()

#WWTP E####

#AS+TWW
tiff("Venn_E_AS_TWW.tiff", units="cm", width=30, height=25, res=300)
select_rows<-c(21:25)
df<-legionella_matrix_bi_combined[select_rows,]

venn_list <- apply(df, 1, function(x) {
  names(x)[x == 1]
},simplify=FALSE)

# Apply() returned a matrix-like object, convert to list
venn_list <- as.list(venn_list)

# Name each list element by sample name (rownames of df)
names(venn_list) <- SampleNames_combined[select_rows]

#Venn-plot
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  alpha = 0.5,
  margin=0.1,
  cex = 1.5,
  cat.cex = 1.2,
  cat.dist =0.3,
)
grid::grid.draw(venn.plot)

dev.off()

#No Winter###
#WWTPA no winter####

#AS+TWW
tiff("Venn_A_AS_TWW_noWinter.tiff", units="cm", width=30, height=25, res=300)
select_rows<-c(1:4)
df<-legionella_matrix_bi_combined[select_rows,]

venn_list <- apply(df, 1, function(x) {
  names(x)[x == 1]
},simplify=FALSE)

venn_list <- as.list(venn_list)

# Name each list element by sample name (rownames of df)
names(venn_list) <- SampleNames_combined[select_rows]

#Venn-plot
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  alpha = 0.5,
  margin=0.1,
  cex = 1.5,
  cat.cex = 1.2,
  cat.dist =0.3,
)
grid::grid.draw(venn.plot)

dev.off()
#WWTP D no winter####
#AS+TWW
tiff("Venn_D_AS_TWW_noWinter.tiff", units="cm", width=30, height=25, res=300)
select_rows<-c(16:19)
df<-legionella_matrix_bi_combined[select_rows,]

venn_list <- apply(df, 1, function(x) {
  names(x)[x == 1]
},simplify=FALSE)

venn_list <- as.list(venn_list)

# Name each list element by sample name (rownames of df)
names(venn_list) <- SampleNames_combined[select_rows]

#Venn-plot
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  alpha = 0.5,
  margin=0.1,
  cex = 1.5,
  cat.cex = 1.2,
  cat.dist =0.3,
)
grid::grid.draw(venn.plot)

dev.off()
#WWTP E only Samples from 2023####

#AS+TWW
tiff("Venn_E_AS_TWW_2023.tiff", units="cm", width=30, height=25, res=300)
select_rows<-c(21:24)
df<-legionella_matrix_bi_combined[select_rows,]

venn_list <- apply(df, 1, function(x) {
  names(x)[x == 1]
},simplify=FALSE)

# Apply() returned a matrix-like object, convert to list
venn_list <- as.list(venn_list)

# Name each list element by sample name (rownames of df)
names(venn_list) <- SampleNames_combined[select_rows]

#Venn-plot
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  alpha = 0.5,
  margin=0.1,
  cex = 1.5,
  cat.cex = 1.2,
  cat.dist =0.3,
)
grid::grid.draw(venn.plot)

dev.off()


#Create an Excel containg ASV_ID, species with highest pident and pident####
ASV_ID<-unique(legionella_ASVs)
Leg_ASVs<-as.data.frame(ASV_ID)
Excel_ASV_species <-Highest_percentage_ID_above_93 %>%
  group_by(ASV_ID) %>%
  filter(pident_corr == max(pident_corr)) %>%
  ungroup() %>%
  group_by(ASV_ID) %>%
  summarise(Species = paste(unique(Species),collapse = ", "),pident=dplyr::first(pident_corr)) %>%
  right_join(Leg_ASVs, by = "ASV_ID")

write.xlsx(Excel_ASV_species,"Leg_ASVs_species.xlsx")

#Create Table with Species, ASVs and frequency per WWTP####

Species_assigned <-Highest_percentage_ID_above_93 %>%
  group_by(ASV_ID) %>%
  filter(pident_corr == max(pident_corr, na.rm = TRUE)) %>%
  summarise(
    max_pident = max(pident_corr, na.rm = TRUE),
    Species = ifelse(
      max_pident < 97,
      "Legionella spp.",
      paste(unique(Species), collapse = ", ")
    ),
    .groups = "drop"
  ) %>%
  select(ASV_ID, Species)

legionella_matrix_all<-ASV_table_noSampleNames[,colnames(ASV_table_noSampleNames) %in% legionella_ASVs ]
legionella_matrix_bi_all<-as.data.frame(lapply(legionella_matrix_all, function(x) as.integer(x > 0)))

legionella_matrix_bi_all_combined <- matrix(0, nrow = 25, ncol = ncol(legionella_matrix))
colnames(legionella_matrix_bi_all_combined) <- colnames(legionella_matrix_all)
for (i in 1:25) {
  # Check if there is any value > 0 in either row i or row i+24 in the same column
  legionella_matrix_bi_all_combined[i, ] <- as.numeric(
    apply(legionella_matrix_all[c(i, i + 25), ], 2, function(x) any(x > 0))
  )
}

mat_long <- legionella_matrix_bi_all_combined %>%
  as.data.frame() %>%
  mutate(Sample = row_number()) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "ASV_ID",
    values_to = "Present"
  )

mat_long <- mat_long %>%
  mutate(Group = case_when(
    Sample %in% 1:5 ~ "A",
    Sample %in% 6:10 ~ "B",
    Sample %in% 11:15 ~ "C",
    Sample %in% 16:20 ~ "D",
    Sample %in% 21:25 ~ "E"
  ))


counts_group <- mat_long %>%
  group_by(Group, ASV_ID) %>%
  summarise(Count = sum(Present), .groups = "drop")


counts_species <- counts_group %>%
  left_join(Species_assigned, by = "ASV_ID")


counts_species <- counts_species %>%
  mutate(Label = ifelse(Count > 0,
                        paste0(ASV_ID, " (", Count, ")"),
                        NA))

frequency_table_wwtp <- counts_species %>%
  select(Species, Group, Label) %>%
  pivot_wider(
    names_from = Group,
    values_from = Label,
    values_fn = ~ paste(na.omit(.x), collapse = ", "),
    values_fill = ""
  ) %>%
  arrange(Species)


write.xlsx(frequency_table_wwtp,"Table_frequency_WWTP.xlsx")

#Create a table with species frequency in A, B, C, D, E and over all
mat_species <- mat_long %>%
  left_join(Species_assigned, by = "ASV_ID")

species_freq <- mat_species %>%
  group_by(Species, Group, Sample) %>% 
  summarise(Present = max(Present), .groups = "drop") %>%  # species present in sample if any ASV present
  group_by(Species, Group) %>%
  summarise(Frequency = sum(Present), .groups = "drop") %>%
  pivot_wider(
    names_from = Group,
    values_from = Frequency,
    values_fill = 0
  )

# Add total (overall) frequency
species_freq <- species_freq %>%
  mutate(overall = A + B + C + D + E)
