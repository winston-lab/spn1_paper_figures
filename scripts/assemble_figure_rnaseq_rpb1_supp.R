library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(antisense_single_locus_datavis_rdata,
                splicing_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2),
                   c(1,1,1,1,1,2,2,2,2,2,2,2))

    load(antisense_single_locus_datavis_rdata)
    load(splicing_rdata)

    figure_rnaseq_rpb1_supp = arrangeGrob(rnaseq_single_locus_datavis,
                                          splicing_volcano,
                                          layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_rnaseq_rpb1_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(antisense_single_locus_datavis_rdata = snakemake@input[["antisense_single_locus_datavis"]],
     splicing_rdata = snakemake@input[["splicing"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

