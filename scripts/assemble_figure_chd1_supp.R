library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(h3k36me3_single_locus_datavis_rdata,
                chd1_qpcr_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3),
                   c(3,3,3,3,3,3,3,3,3,3,3,3))

    load(h3k36me3_single_locus_datavis_rdata)
    load(chd1_qpcr_rdata)

    temp_panel = ggplot() +
        labs(tag="c") +
        theme_void() +
        theme(plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"),
              plot.tag=element_text(size=9,
                                    family="FreeSans",
                                    face="bold"))

    figure_chd1_supp = arrangeGrob(h3k36me3_single_locus_datavis,
                                   chd1_qpcr,
                                   temp_panel,
                                   layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_chd1_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(h3k36me3_single_locus_datavis_rdata = snakemake@input[["h3k36me3_single_locus_datavis"]],
     chd1_qpcr_rdata= snakemake@input[["chd1_qpcr"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

