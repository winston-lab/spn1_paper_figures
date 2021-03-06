library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(spn1_depletion_western_rdata,
                spn1_depletion_chipseq_barplot_rdata,
                spn1_depletion_metagene_rdata,
                spn1_depletion_scatter_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,4,4,4,4,4,4,4,4),
                   c(1,1,1,1,4,4,4,4,4,4,4,4),
                   c(1,1,1,1,4,4,4,4,4,4,4,4),
                   c(1,1,1,1,4,4,4,4,4,4,4,4),
                   c(1,1,1,1,4,4,4,4,4,4,4,4),
                   c(1,1,1,1,4,4,4,4,4,4,4,4),
                   c(2,2,2,2,3,3,3,3,3,3,3,3),
                   c(2,2,2,2,3,3,3,3,3,3,3,3),
                   c(2,2,2,2,3,3,3,3,3,3,3,3),
                   c(2,2,2,2,3,3,3,3,3,3,3,3),
                   c(2,2,2,2,3,3,3,3,3,3,3,3),
                   c(2,2,2,2,3,3,3,3,3,3,3,3))
    # layout = rbind(c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3),
    #                c(1,1,4,4,4,4,2,2,3,3,3,3))

    load(spn1_depletion_western_rdata)
    load(spn1_depletion_chipseq_barplot_rdata)
    spn1_depletion_chipseq_barplot = chipseq_abundance_barplot
    load(spn1_depletion_metagene_rdata)
    load(spn1_depletion_scatter_rdata)

    figure_spn1_depletion = arrangeGrob(spn1_depletion_western,
                                        spn1_depletion_chipseq_barplot,
                                        spn1_depletion_scatter,
                                        spn1_depletion_metagene,
                                        layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_spn1_depletion,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(spn1_depletion_western_rdata = snakemake@input[["spn1_depletion_western"]],
     spn1_depletion_chipseq_barplot_rdata = snakemake@input[["spn1_depletion_chipseq_barplot"]],
     spn1_depletion_metagene_rdata = snakemake@input[["spn1_depletion_metagene"]],
     spn1_depletion_scatter_rdata = snakemake@input[["spn1_depletion_scatter"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

