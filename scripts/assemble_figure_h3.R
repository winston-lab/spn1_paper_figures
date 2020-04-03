library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(h3_vs_rpb1_ma_rdata,
                reduced_h3_h3_metagene_rdata,
                h3_single_locus_datavis_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3))

    load(h3_vs_rpb1_ma_rdata)
    load(reduced_h3_h3_metagene_rdata)
    load(h3_single_locus_datavis_rdata)

    # temp_panel = ggplot() +
    #     annotate(geom="text",
    #              x=0,
    #              y=0,
    #              label="single gene H3 examples",
    #              size=2) +
    #     labs(tag="c") +
    #     theme_void() +
    #     theme(plot.tag=element_text(size=9,
    #                                 face="bold",
    #                                 family="FreeSans"),
    #           plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"))

    figure_h3 = arrangeGrob(h3_vs_rpb1_ma,
                            reduced_h3_h3_metagene,
                            h3_single_locus_datavis,
                            layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_h3,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(h3_vs_rpb1_ma_rdata = snakemake@input[["h3_vs_rpb1_ma"]],
     reduced_h3_h3_metagene_rdata = snakemake@input[["reduced_h3_h3_metagene"]],
     h3_single_locus_datavis_rdata = snakemake@input[["h3_single_locus_datavis"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

