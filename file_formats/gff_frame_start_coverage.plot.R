#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)
if (length(args) != 3) {
  write("\nusage: gff_frame_start_coverage.plot.R infile out.png title\n", stderr())
  quit(status = 1)
}

df <- read.csv(args[1], sep = "\t")
windows <- df['end'] - df['start']
window <- windows[1,1] + 1

if (window == 1000) {
  ylabel = paste("Counts/", 1, "kb", sep='')
} else if (window > 1000) {
  ylabel = paste("Counts/", window/1000, "kb", sep='')
} else {
  ylabel = paste("Counts/", window, "bp", sep='')
}

df <- select(df, X.chr, strand, cnt_f0, cnt_f1, cnt_f2)

df_m <- melt(df, id.vars = c("X.chr", "strand"))

p <- ggplot(df_m, aes(variable, value, fill = strand)) +
  geom_violin(adjust=1, position = position_dodge(width = 0.75)) +
  scale_x_discrete(labels = c('0','1','2')) +
  xlab('Frame') +
  ylab(ylabel) +
  ggtitle(args[3]) + 
  facet_grid(. ~ X.chr) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key = element_blank(),
    strip.background = element_rect(
      colour = "white", fill = "white",
      size = 0.2
    )
  )

ggsave(p, file = args[2], width = 8, height = 4)
