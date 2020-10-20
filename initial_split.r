#!/usr/bin/env Rscript

library("argparser",quietly=TRUE)

# Create a parser
p <- arg_parser("Split survival data into test/train sets")
# Add command line arguments
p <- add_argument(p, "file", help="number to round", default="",type="character")

argv <- parse_args(p)
file <- argv$file
if (file=="") {
  print(p)
  
}