library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(promoter_swap_rtqpcr_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2),
                   c(2,2,2,2,2,2,2,2,2,2,2,2))

    load(promoter_swap_rtqpcr_rdata)

    temp_panel = ggplot() +
        annotate(geom="text",
                 x=0,
                 y=0,
                 label="Schematic of promoter swap strains.",
                 size=2) +
        labs(tag="a") +
        theme_void() +
        theme(plot.tag=element_text(size=9,
                                    face="bold",
                                    family="FreeSans"),
              plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"))

    figure_promoter_swap = arrangeGrob(temp_panel,
                                       promoter_swap_scatter,
                                       layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_promoter_swap,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(promoter_swap_rtqpcr_rdata = snakemake@input[["promoter_swap_rtqpcr"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

