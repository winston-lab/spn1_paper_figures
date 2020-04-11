library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(rnaseq_maplot_rdata,
                # rnaseq_single_locus_datavis_rdata,
                rnaseq_vs_rpb1_single_locus_rdata,
                rpb1_metagenes_rdata,
                rnaseq_vs_rpb1_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(3,3,3,3,3,3,NA,2,2,2,2,2),
                   c(3,3,3,3,3,3,NA,2,2,2,2,2),
                   c(3,3,3,3,3,3,NA,2,2,2,2,2),
                   c(3,3,3,3,3,3,NA,2,2,2,2,2),
                   c(4,4,4,4,4,4,4,5,5,5,5,5),
                   c(4,4,4,4,4,4,4,5,5,5,5,5),
                   c(4,4,4,4,4,4,4,5,5,5,5,5),
                   c(4,4,4,4,4,4,4,5,5,5,5,5))

    load(rnaseq_maplot_rdata)
    # load(rnaseq_single_locus_datavis_rdata)
    load(rnaseq_vs_rpb1_single_locus_rdata)
    load(rpb1_metagenes_rdata)
    load(rnaseq_vs_rpb1_rdata)

    figure_rnaseq_rpb1 = arrangeGrob(rnaseq_maplot,
                                     rpb1_metagenes,
                                     rnaseq_v_rpb1,
                                     rnaseq_vs_rpb1_single_locus,
                                     nullGrob(),
                                     # rnaseq_single_locus_datavis,
                                     layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_rnaseq_rpb1,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(rnaseq_maplot_rdata = snakemake@input[["rnaseq_maplot"]],
     # rnaseq_single_locus_datavis_rdata = snakemake@input[["rnaseq_single_locus_datavis"]],
     rnaseq_vs_rpb1_single_locus_rdata = snakemake@input[["rnaseq_vs_rpb1_single_locus"]],
     rpb1_metagenes_rdata = snakemake@input[["rpb1_metagenes"]],
     rnaseq_vs_rpb1_rdata = snakemake@input[["rnaseq_vs_rpb1"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

