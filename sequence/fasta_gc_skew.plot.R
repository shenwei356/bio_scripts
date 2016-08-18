#!/usr/bin/env Rscript
# https://github.com/shenwei356/bio_scripts
library(methods)
library(proto)
library(argparse)
library(ggplot2)
library(reshape2)
library(scales)
library(ggthemes)

shenwei356.theme <-  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(colour = "black", size = 0.8),
    axis.line.y = element_line(colour = "black", size = 0.8),
    axis.ticks.x = element_line(size = 0.8),
    axis.ticks.y = element_line(size = 0.8),
    #     axis.text.x = element_text(
    #       angle = 30, hjust = 1, vjust = 1
    #     ),
    
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "plain"),
    legend.background = element_rect(fill = "transparent"),
    strip.background = element_rect(
      colour = "white",
      fill = "white",
      size = 0.2
    ),
    strip.text.x = element_text(size = 16),
    
    text = element_text(
      size = 14,
      family = "arial",
      face = "bold"
    ),
    plot.title = element_text(
      size = 16,
      family = "arial",
      face = "bold"
    )
  )


parser <-
  ArgumentParser(description = "Plot GC and GC Skew with the result produced by fasta_gc_skew.py",
                 formatter_class = "argparse.RawTextHelpFormatter")

parser$add_argument("infile", type = "character",
                    help = "gcskew file produced by fasta_gc_skew.py")
parser$add_argument("outfile", type = "character",
                    help = "outfile")
parser$add_argument(
  "-xi",
  "--x-interval",
  type = "integer",
  default = 1000000,
  help = "x axix interval [1,000,000]"
)
parser$add_argument("-n",
                    type = "integer",
                    default = 10,
                    help = "divide the normalized accum_gcskew by n so it looks better [10]")
parser$add_argument(
  "--width",
  metavar = "width",
  type = "integer",
  default = 20,
  help = "output image width [20]"
)
parser$add_argument(
  "--height",
  metavar = "height",
  type = "integer",
  default = 5,
  help = "output image height [5]"
)
parser$add_argument(
  "-g",
  "--gc-content",
  action = "store_true",
  dest = "gc_content",
  help = "only plot GC Content"
)
parser$add_argument("-s",
                    "--gc-skew",
                    action = "store_true",
                    dest = "gc_skew",
                    help = "only plot GC Skew")
parser$add_argument(
  "-t",
  "--title",
  metavar = "title",
  type = "character",
  default = "GC Content/GC Skew",
  help = "title"
)

args <- parser$parse_args()

if (args$title == "") {
  args$title = NULL
}

df <- read.csv(args$infile, sep = "\t")
df['accum_gcskew'] = df['accum_gcskew'] / max(df['accum_gcskew']) / args$n

if (args$gc_content && !args$gc_skew) {
  df['gcskew'] = NULL
  df['accum_gcskew'] = NULL
}
if (!args$gc_content && args$gc_skew) {
  df['gc'] = NULL
}

df_m <- melt(df, id.vars = c("chr", "loc"))

p <- ggplot(df_m) +
  geom_line(aes(loc, value, color = variable)) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_size(range = c(0.1)) +
  scale_colour_stata() +
  facet_grid(chr ~ .) +
  ylab(NULL) +
  xlab("Position(bp)") +
  scale_x_continuous(breaks = seq(0, max(df$loc), by = args$x_interval),
                     labels = comma) +
  ggtitle(args$title) +
  shenwei356.theme +
  theme(legend.position = "top")

ggsave(
  p,
  file = args$outfile,
  width = args$width,
  height = args$height
)
