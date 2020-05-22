library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(histone_spn1_depletion_westerns_rdata,
                chipseq_abundance_barplots_h3_rdata,
                h3_mods_non_h3_norm_rdata,
                h3_mods_facet_expression_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,2,2,2,2,2,NA,NA,NA),
                   c(1,1,1,1,2,2,2,2,2,NA,NA,NA),
                   c(1,1,1,1,2,2,2,2,2,NA,NA,NA),
                   c(1,1,1,1,3,3,3,3,3,3,3,3),
                   c(1,1,1,1,3,3,3,3,3,3,3,3),
                   c(NA,NA,NA,NA,3,3,3,3,3,3,3,3),
                   c(4,4,4,4,4,4,4,4,4,5,5,5),
                   c(4,4,4,4,4,4,4,4,4,5,5,5),
                   c(4,4,4,4,4,4,4,4,4,5,5,5),
                   c(4,4,4,4,4,4,4,4,4,5,5,5),
                   c(4,4,4,4,4,4,4,4,4,5,5,5),
                   c(4,4,4,4,4,4,4,4,4,5,5,5))

    load(histone_spn1_depletion_westerns_rdata)
    load(chipseq_abundance_barplots_h3_rdata)
    load(h3_mods_non_h3_norm_rdata)
    load(h3_mods_facet_expression_rdata)

    figure_h3_mods = arrangeGrob(histone_spn1_depletion_westerns,
                                 chipseq_abundance_barplot,
                                 h3_mods_non_h3_norm,
                                 h3_mods_facet_expression,
                                 nullGrob(),
                                 layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_h3_mods,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(histone_spn1_depletion_westerns_rdata = snakemake@input[["histone_spn1_depletion_westerns"]],
     chipseq_abundance_barplots_h3_rdata = snakemake@input[["chipseq_abundance_barplots_h3"]],
     h3_mods_non_h3_norm_rdata = snakemake@input[["h3_mods_non_h3_norm"]],
     h3_mods_facet_expression_rdata = snakemake@input[["h3_mods_facet_expression"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

