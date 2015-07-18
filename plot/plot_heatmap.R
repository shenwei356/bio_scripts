#!/usr/bin/env Rscript
library(proto)
library(argparse)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

description <- "Plot heatmap"
parser <-
  ArgumentParser(description = description, formatter_class = "argparse.RawTextHelpFormatter")

parser$add_argument("infile", type = "character", help = "infile (tsv, with head")
parser$add_argument("-H", "--header", action = "store_false", dest = "header",
                    help = "header")
parser$add_argument(
  "-F", "--field-seperator", metavar = "field_seperator", type = "character",
  default = "\t", help = "field seperator"
)

parser$add_argument("outname", type = "character", help = "outname")
parser$add_argument(
  "--title", metavar = "title", type = "character", default = "",
  help = "title"
)
parser$add_argument(
  "--width", metavar = "width", type = "integer", default = 5,
  help = "output image width"
)
parser$add_argument(
  "--height", metavar = "height", type = "integer", default = 5,
  help = "output image height"
)
parser$add_argument(
  "-thr", "--treeheight_row", metavar = "treeheight_row", type = "integer",
  default = 30,   help = "treeheight_row"
)
parser$add_argument(
  "-thc", "--treeheight_col", metavar = "treeheight_col", type = "integer",
  default = 30,   help = "treeheight_col"
)
parser$add_argument(
  "-fo", "--fontsize", metavar = "fontsize", type = "integer",
  default = 10,   help = "fontsize"
)
parser$add_argument(
  "-fr", "--fontsize_row", metavar = "fontsize_row", type = "integer",
  default =   10,
  help = "fontsize_row"
)
parser$add_argument(
  "-fc", "--fontsize_col", metavar = "fontsize_col", type = "integer",
  default =     10,
  help = "fontsize_col"
)

args <- parser$parse_args()

data = read.csv(args$infile, header = args$header, sep = args$field_seperator)#, row.names = 'NAME')
df <- data[,-1]
rownames(df) <- data[,1]

png(
  filename = paste(args$outname, 'png', sep = "."),
  width = args$width, height = args$height, units = "in", res = 300
)

pheatmap(
  df,
  main = args$title,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = args$treeheight_row,
  treeheight_col = args$treeheight_col,
  color =  colorRampPalette(rev(brewer.pal(10,"RdYlBu")))(256),
  border_color = "white",
  fontsize = args$fontsize,
  fontsize_row = args$fontsize_row,
  fontsize_col = args$fontsize_col
)

dev.off()