library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(differential_expression_rtqpcr_rdata,
                antisense_single_locus_datavis_rdata,
                chipseq_abundance_barplots_rpb1_rdata,
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
                   c(3,3,3,3,3,2,2,2,2,2,2,2),
                   c(3,3,3,3,3,2,2,2,2,2,2,2),
                   c(3,3,3,3,3,2,2,2,2,2,2,2),
                   c(3,3,3,3,3,2,2,2,2,2,2,2),
                   c(3,3,3,3,3,2,2,2,2,2,2,2))

    load(differential_expression_rtqpcr_rdata)
    load(antisense_single_locus_datavis_rdata)
    load(chipseq_abundance_barplots_rpb1_rdata)

    figure_rnaseq_rpb1_supp = arrangeGrob(diffexp_rtqpcr_scatter,
                                          rnaseq_single_locus_datavis,
                                          chipseq_abundance_barplot,
                                          layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_rnaseq_rpb1_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(differential_expression_rtqpcr_rdata = snakemake@input[["differential_expression_rtqpcr"]],
     antisense_single_locus_datavis_rdata = snakemake@input[["antisense_single_locus_datavis"]],
     chipseq_abundance_barplots_rpb1_rdata = snakemake@input[["chipseq_abundance_barplots_rpb1"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

