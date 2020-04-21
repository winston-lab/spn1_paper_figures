library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(h3_plmin_reduced_lfc_scatterplots_rdata,
                fig_width=17.4,
                fig_height=12,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2))

    load(h3_plmin_reduced_lfc_scatterplots_rdata)

    figure_h3_supp = arrangeGrob(h3_plmin_reduced_lfc_scatterplots,
                                nullGrob(),
                                layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_h3_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(h3_plmin_reduced_lfc_scatterplots_rdata = snakemake@input[["h3_plmin_reduced_lfc_scatterplots"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

