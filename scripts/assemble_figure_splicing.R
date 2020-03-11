library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(splicing_rdata,
                splicing_rtqpcr_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2))

    load(splicing_rdata)
    load(splicing_rtqpcr_rdata)

    figure_splicing = arrangeGrob(splicing_volcano,
                                        splicing_rtqpcr_barplot,
                                        layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_splicing,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(splicing_rdata = snakemake@input[["splicing"]],
     splicing_rtqpcr_rdata = snakemake@input[["splicing_rtqpcr"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

