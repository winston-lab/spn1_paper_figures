library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(chipseq_abundance_barplots_h3_rdata,
                h3_mods_non_h3_norm_rdata,
                h3_mods_facet_expression_rdata,
                # h3k4me3_rdata,
                # h3k36me2_rdata,
                # h3k36me3_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,NA,NA,NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,NA,NA,NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,NA,NA,NA,NA,NA,NA,NA),
                   c(2,2,2,2,2,2,2,2,2,4,4,4),
                   c(2,2,2,2,2,2,2,2,2,4,4,4),
                   c(2,2,2,2,2,2,2,2,2,4,4,4),
                   c(3,3,3,3,3,3,3,3,3,4,4,4),
                   c(3,3,3,3,3,3,3,3,3,4,4,4),
                   c(3,3,3,3,3,3,3,3,3,4,4,4),
                   c(3,3,3,3,3,3,3,3,3,4,4,4),
                   c(3,3,3,3,3,3,3,3,3,4,4,4),
                   c(3,3,3,3,3,3,3,3,3,4,4,4))

    load(chipseq_abundance_barplots_h3_rdata)
    load(h3_mods_non_h3_norm_rdata)
    load(h3_mods_facet_expression_rdata)
    # load(h3k4me3_rdata)
    # h3k4me3 = h3_mods_facet_expression_length
    # load(h3k36me2_rdata)
    # h3k36me2 = h3_mods_facet_expression_length
    # load(h3k36me3_rdata)
    # h3k36me3 = h3_mods_facet_expression_length

    figure_h3_mods = arrangeGrob(chipseq_abundance_barplot,
                                 h3_mods_non_h3_norm,
                                 h3_mods_facet_expression,
                                 nullGrob(),
                                 # h3k36me2,
                                 # h3k36me3,
                                 # h3k4me3,
                                 layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_h3_mods,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(chipseq_abundance_barplots_h3_rdata = snakemake@input[["chipseq_abundance_barplots_h3"]],
     h3_mods_non_h3_norm_rdata = snakemake@input[["h3_mods_non_h3_norm"]],
     h3_mods_facet_expression_rdata = snakemake@input[["h3_mods_facet_expression"]],
     # h3k4me3_rdata = snakemake@input[["h3k4me3"]],
     # h3k36me2_rdata = snakemake@input[["h3k36me2"]],
     # h3k36me3_rdata = snakemake@input[["h3k36me3"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

