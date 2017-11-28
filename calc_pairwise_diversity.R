library(Biostrings)
library(ape)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsignif)
source("two_sample_stat_jp.R")

patient <- '601'
rebound_timepoints <- list("601"="W29", "610"="W31", "613"="W37")
timepoints <- c("D14","W23", rebound_timepoints[[patient]])

fasta_file <- paste("~/data_alpha/home/jpai/julio/pairwise_distance/",patient,"_Align_092017.fasta",sep="")
sequences <- readDNAStringSet(fasta_file)
seq_length <- unique(width(sequences))
stopifnot(length(seq_length) == 1)

names(sequences) <- gsub("-","_",names(sequences))
seq_ids <- names(sequences)

all_dists <- data.frame(timepoint=character(0), genetic_distance=numeric(0))
all_dist_dfs <- list()

for (timepoint in timepoints) {
  print(timepoint)
  grouped_seq_ids <- seq_ids[grepl(paste("_",timepoint,"_",sep = ""), seq_ids)]
  grouped_sequences <- as.DNAbin(sequences[grouped_seq_ids])

  d <- dist.dna(grouped_sequences, model="JC69", as.matrix = TRUE)
  stopifnot(rownames(d) == colnames(d))
  rownames(d) <- c(1:nrow(d))
  colnames(d) <- c(1:ncol(d))
  df <- melt(d, varnames=c("row","col"))[melt(upper.tri(d))$value,]
  df$length <- seq_length
  
  all_dist_dfs[[timepoint]] <- df
  all_dists <- rbind(all_dists, data.frame(timepoint=timepoint, genetic_distance=df$value))
}

# compute p-values using two sample-statistics Z test
pvalues <- list()
for (comb in combn(timepoints, 2, simplify = FALSE)) {
  key <- paste(comb[1], comb[2], sep="_")
  pvalues[[key]] <- two_sample_test(all_dist_dfs[comb])
}

my_comparisons <- combn(timepoints, 2, simplify = FALSE)
pvalue_annot <- sapply(pvalues, function(x) { formatC(x$p.ztest, digits=3) } )
names(pvalue_annot) <- NULL

ggboxplot(all_dists, x = "timepoint", y = "genetic_distance", color = "timepoint", palette = c("darkgreen","navy","red2"), 
          add = "point", width = 0.1, add.params = list(alpha=0.3)) + 
  coord_cartesian(ylim = c(0, 0.1), expand = FALSE) + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  ggtitle(paste("pairwise genetic distance (patient ", patient,")",sep = "")) + ylab("pairwise diversity (nt)") +
  geom_signif(comparisons = my_comparisons, textsize = 4, annotations = pvalue_annot, y_position = c(0.083, 0.087, 0.095))