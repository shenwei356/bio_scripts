#!/usr/bin/env Rscript
library(proto)
library(argparse)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

description <- paste("Plot heatmap. ", 
                     "Infile should be a csv/tsv file with header containing",
                     " column names. Annotation for row is also supported",
                     ", please put them in the last column.", sep="")
parser <- ArgumentParser(description = description,
                         formatter_class = "argparse.RawTextHelpFormatter")

#-----------------------------------------------------------------------------

parser$add_argument("infile", type = "character",
                    help = "infile (tsv, with head")
parser$add_argument("-H", "--header", action = "store_false",
                    dest = "header", help = "header")
parser$add_argument(
  "-F", "--field-seperator", metavar = "field_seperator", type = "character",
  default = "\t", help = "field seperator"
)
parser$add_argument(
  "-a", "--with-annot", action = "store_true", dest = "with_annotation",
  help = "add annotation_row from the last column"
)
parser$add_argument(
  "-al", "--with-annot-legend", action = "store_true",
  dest = "with_annotation_legend", help = "show annotation_row_legend"
)

#-----------------------------------------------------------------------------

parser$add_argument("outname", type = "character", help = "outname")
parser$add_argument(
  "--title", metavar = "title", type = "character", default = "",
  help = "title"
)

#-----------------------------------------------------------------------------

parser$add_argument(
  "-s", "--scale", metavar = "scale", type = "character",
  default = "row", help = "scale. row | column | none [row]"
)

parser$add_argument(
  "-ncr", "--not-cluster-rows", action = "store_false", dest = "not_cluster_rows",
  help = "do not cluster_rows"
)
parser$add_argument(
  "-ncc", "--not-cluster-cols", action = "store_false", dest = "not_cluster_cols",
  help = "do not cluster_cols"
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
  default = 30, help = "treeheight_row"
)
parser$add_argument(
  "-thc", "--treeheight_col", metavar = "treeheight_col", type = "integer",
  default = 30, help = "treeheight_col"
)
parser$add_argument(
  "-fo", "--fontsize", metavar = "fontsize", type = "integer",
  default = 10, help = "fontsize"
)
parser$add_argument(
  "-fr", "--fontsize_row", metavar = "fontsize_row", type = "integer",
  default = 10, help = "fontsize_row"
)
parser$add_argument(
  "-fc", "--fontsize_col", metavar = "fontsize_col", type = "integer",
  default = 10, help = "fontsize_col"
)

#-----------------------------------------------------------------------------

args <- parser$parse_args()

if (! (args$scale == "row" || args$scale == "column" || args$scale == "none")){
  write("value of option -s/--scale should be in [row, column, none]", file=stderr())
  quit(1)
}

#-----------------------------------------------------------------------------

data = read.csv(args$infile, header = args$header, sep = args$field_seperator)
df <- data[,-1]
rownames(df) <- data[,1]

#-----------------------------------------------------------------------------

annotation_row = NA
with_annotation_legend = FALSE
if (args$with_annotation) {
  annotation_row <- df[,NCOL(df)]
  df[,NCOL(df)] <- NULL
  annotation_row = as.data.frame(annotation_row)
  rownames(annotation_row) <- rownames(df)
  colnames(annotation_row) <- " "
  with_annotation_legend = args$with_annotation_legend
}

#-----------------------------------------------------------------------------

png(
  filename = paste(args$outname, 'png', sep = "."),
  width = args$width, height = args$height, units = "in", res = 300
)

pheatmap(
  df,
  main = args$title,
  scale = args$scale,
  cluster_rows = args$not_cluster_rows,
  cluster_cols = args$not_cluster_cols,
  
  treeheight_row = args$treeheight_row,
  treeheight_col = args$treeheight_col,
  color =  colorRampPalette(rev(brewer.pal(10,"RdYlBu")))(256),
  border_color = "white",
  fontsize = args$fontsize,
  fontsize_row = args$fontsize_row,
  fontsize_col = args$fontsize_col,
  
  annotation_row = annotation_row,
  annotation_legend = with_annotation_legend
)

dev.off()