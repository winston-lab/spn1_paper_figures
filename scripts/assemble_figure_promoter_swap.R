library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(promoter_swap_diagram_rdata,
                promoter_swap_rtqpcr_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
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

    load(promoter_swap_diagram_rdata)
    load(promoter_swap_rtqpcr_rdata)

    figure_promoter_swap = arrangeGrob(promoter_swap_diagram,
                                       promoter_swap_scatter,
                                       layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_promoter_swap,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(promoter_swap_diagram_rdata = snakemake@input[["promoter_swap_diagram"]],
     promoter_swap_rtqpcr_rdata = snakemake@input[["promoter_swap_rtqpcr"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

