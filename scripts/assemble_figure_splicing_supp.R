library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(rpgene_datavis_supp_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2,2))

    load(rpgene_datavis_supp_rdata)

    figure_splicing_supp = arrangeGrob(rpgene_datavis,
                                 nullGrob(),
                                 layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_splicing_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(rpgene_datavis_supp_rdata = snakemake@input[["rpgene_datavis_supp"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

