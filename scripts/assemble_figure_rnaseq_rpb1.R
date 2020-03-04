library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(rnaseq_maplot_rdata,
                rpb1_metagenes_rdata,
                rnaseq_vs_rpb1_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,3,3,3,3,3,3),
                   c(1,1,1,1,1,1,3,3,3,3,3,3),
                   c(1,1,1,1,1,1,3,3,3,3,3,3),
                   c(1,1,1,1,1,1,3,3,3,3,3,3),
                   c(1,1,1,1,1,1,3,3,3,3,3,3),
                   c(2,2,2,2,2,2,3,3,3,3,3,3),
                   c(2,2,2,2,2,2,3,3,3,3,3,3),
                   c(2,2,2,2,2,2,3,3,3,3,3,3),
                   c(2,2,2,2,2,2,4,4,4,4,4,4),
                   c(2,2,2,2,2,2,4,4,4,4,4,4),
                   c(2,2,2,2,2,2,4,4,4,4,4,4),
                   c(2,2,2,2,2,2,4,4,4,4,4,4))

    load(rnaseq_maplot_rdata)
    load(rpb1_metagenes_rdata)
    load(rnaseq_vs_rpb1_rdata)

    temp_panel = ggplot() +
        annotate(geom="text",
                 x=0,
                 y=0,
                 label="RNA-seq single loci",
                 size=2) +
        labs(tag="b") +
        theme_void() +
        theme(plot.tag=element_text(size=8,
                                    face="bold",
                                    family="FreeSans"))

    figure_rnaseq_rpb1 = arrangeGrob(rnaseq_maplot,
                                     temp_panel,
                                     rpb1_metagenes,
                                     rnaseq_v_rpb1,
                                     layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_rnaseq_rpb1,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(rnaseq_maplot_rdata = snakemake@input[["rnaseq_maplot"]],
     rpb1_metagenes_rdata = snakemake@input[["rpb1_metagenes"]],
     rnaseq_vs_rpb1_rdata = snakemake@input[["rnaseq_vs_rpb1"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

