#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)
if (length(args) != 2) {
  write("\nusage: gff_frame_start_coverage.plot.R infile out.png\n", stderr())
  quit(status=1)
}

df <- read.csv(args[1], sep = "\t")
# print(head(df))
df <- select(df, X.chr, strand, cnt_f0, cnt_f1, cnt_f2)
# print(head(df))

df_m <- melt(df, id.vars = c("X.chr", "strand"))
# print(head(df_m))

p <- ggplot(df_m) +
  geom_boxplot(aes(variable, value, fill = strand),
               position = position_dodge(width = 0.75)) +
  # geom_point(aes(variable, value, color=strand))+
  scale_x_discrete(labels = c('0','1','2')) +
  xlab('Frame') +
  ylab('Counts') +
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

