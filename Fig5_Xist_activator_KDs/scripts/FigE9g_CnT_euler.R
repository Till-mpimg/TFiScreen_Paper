library(tidyverse)
library(GenomicRanges)
library(ChIPpeakAnno)
library(eulerr)


theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6), plot.title = element_text(hjust = 0.5),
                  axis.title = element_blank(), panel.border = element_rect(size = 0.5, color = "black", fill = NA),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  legend.text = element_text(size = 6), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  strip.background = element_blank(), strip.text = element_text(size = 6),
                  legend.key.size = unit(0.5,"line")))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

marks <- c("H3K4me3", "H3K27ac", "H3K4me1")
targets <- c("Zic3", "Nfrkb", "Otx2", "Foxd3")
chr_names <- paste0("chr", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, "X"))


lapply(marks, function (x) {
  
  files <- paste0("./output_files/", targets, "_", x, "_diffbind.bed") 
  mark <- x
  
  diffbind <- lapply(files, readPeakFile)
  diffbind_fil <- lapply(diffbind, function(x) {
    return(x[seqnames(x) %in% chr_names])
  })
  
  diffbind_gained <- lapply(diffbind_fil, function(x) {
    return(x[x$V9 > 0])
  }) 
  names(diffbind_gained) <- paste0(targets, "_gained")
  
  
  diffbind_lost <- lapply(diffbind_fil, function(x) {
    y = x[x$V9 < 0]
    y$V12 = abs(y$V9)
    return(y)
  })
  names(diffbind_lost) <- paste0(targets, "_lost") 
  
  #Calculate and plot overlaps for loss of peaks
  overlaps_lost <- findOverlapsOfPeaks(diffbind_lost)
  venn_cnt_lost <- overlaps_lost$venn_cnt[,1:5]
  sets_lost <- venn_cnt_lost[, -ncol(venn_cnt_lost)]
  counts_lost <- venn_cnt_lost[,5]
  
  combi_names_lost <- apply(sets_lost, 1, function(x) paste(names(x) [x == 1], collapse = "&"))
  
  euler_input_lost <- tapply(counts_lost, combi_names_lost, sum)
  euler_input_lost <- euler_input_lost[names(euler_input_lost) != ""]
  
  
  p1 <- plot(euler(euler_input_lost), quantities = TRUE)
  
  #Plot euler diagram and write to file
  pdf(paste0("./output_files/FigE9g_Euler_lost_", mark, ".pdf"), width = 2, height = 2,
      useDingbats=FALSE)
  print(p1)
  dev.off()
 
  
  #Calculate and plot overlaps for gain of peaks
  overlaps_gained <- findOverlapsOfPeaks(diffbind_gained)
  venn_cnt_gained <- overlaps_gained$venn_cnt[,1:5]
  sets_gained <- venn_cnt_gained[, -ncol(venn_cnt_gained)]
  counts_gained <- venn_cnt_gained[,5]
  
  combi_names_gained <- apply(sets_gained, 1, function(x) paste(names(x) [x == 1], collapse = "&"))
  
  euler_input_gained <- tapply(counts_gained, combi_names_gained, sum)
  euler_input_gained <- euler_input_gained[names(euler_input_gained) != ""]
  
  p2 <- plot(euler(euler_input_gained), quantities = TRUE)
  
  #Plot euler diagram and write to file
  pdf(paste0("./output_files/FigE9g_Euler_gained_", mark, ".pdf"), width = 2, height = 2,
      useDingbats=FALSE)
  print(p2)
  dev.off()   
  
})



