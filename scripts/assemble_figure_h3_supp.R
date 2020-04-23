library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(h3_metagene_rdata,
                h3_vs_rpb1_ma_5p500bp_rdata,
                h3_plmin_reduced_lfc_scatterplots_rdata,
                reduced_h3_matched_metagenes_rdata,
                fig_width=17.4,
                fig_height=12,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,NA,2,2,2,2,2,2),
                   c(1,1,1,1,1,NA,2,2,2,2,2,2),
                   c(1,1,1,1,1,NA,2,2,2,2,2,2),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(4,4,4,4,4,4,4,NA,NA,NA,NA,NA),
                   c(4,4,4,4,4,4,4,NA,NA,NA,NA,NA),
                   c(4,4,4,4,4,4,4,NA,NA,NA,NA,NA),
                   c(4,4,4,4,4,4,4,NA,NA,NA,NA,NA),
                   c(4,4,4,4,4,4,4,NA,NA,NA,NA,NA))

    load(h3_metagene_rdata)
    load(h3_vs_rpb1_ma_5p500bp_rdata)
    load(h3_plmin_reduced_lfc_scatterplots_rdata)
    load(reduced_h3_matched_metagenes_rdata)

    figure_h3_supp = arrangeGrob(h3_metagene,
                                 h3_vs_rpb1_ma_5p500bp,
                                 h3_plmin_reduced_lfc_scatterplots,
                                 reduced_h3_matched_metagenes,
                                 layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_h3_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(h3_metagene_rdata = snakemake@input[["h3_metagene"]],
     h3_vs_rpb1_ma_5p500bp_rdata = snakemake@input[["h3_vs_rpb1_ma_5p500bp"]],
     h3_plmin_reduced_lfc_scatterplots_rdata = snakemake@input[["h3_plmin_reduced_lfc_scatterplots"]],
     reduced_h3_matched_metagenes_rdata = snakemake@input[["reduced_h3_matched_metagenes"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

