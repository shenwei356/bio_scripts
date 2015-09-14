#!/usr/bin/env Rscript
# https://github.com/shenwei356/bio_scripts

library(methods)
library(proto)
library(argparse)
library(ggplot2)
library(reshape2)
library(scales)

#-----------------------------------------------------------------------------

description <- paste(
  "Plot distribution.",
  "Infile should be a tsv file of two columns (group and \"value\")", sep = ""
)

parser <-
  ArgumentParser(description = description,
                 formatter_class = "argparse.RawTextHelpFormatter")

#-----------------------------------------------------------------------------

parser$add_argument("infile", type = "character",
                    help = "infile")
parser$add_argument("outfile", type = "character",
                    help = "outfile")

parser$add_argument(
  "-bw", "--binwidth", type = "double",
  default = 0.1, help = "binwidth"
)

parser$add_argument("--xlab", type = "character", default = "Value",
                    help = "xlabel")
parser$add_argument("--ylab", type = "character", default = "Density",
                    help = "ylabel")
parser$add_argument("--width", type = "integer", default = 6,
                    help = "output image width [20]")
parser$add_argument("--height", type = "integer", default = 3,
                    help = "output image height [5]")

parser$add_argument(
  "-t", "--title", metavar = "title", type = "character",
  default = "", help = "title"
)

#-----------------------------------------------------------------------------

args <- parser$parse_args()

if (args$title == "") {
  args$title = ""
}

#-----------------------------------------------------------------------------

df <- read.csv(args$infile, sep = "\t")

p <- ggplot(df, aes(x = value, fill = group, colour = group)) +
  geom_histogram(
    aes(y = ..density..), alpha = .3, position = "identity", binwidth = args$binwidth
  ) +
  geom_density(alpha = .2) +
  ylab(args$ylab) +
  xlab(args$xlab) +
  ggtitle(args$title) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key = element_blank(),
    # legend.position = "none",
    legend.title = element_blank()
  )

ggsave(
  p, file = args$outfile, width = args$width, height = args$height
)