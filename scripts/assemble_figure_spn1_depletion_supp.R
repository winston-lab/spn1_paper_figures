library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(spn1_depletion_viability_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2),
                   c(1,1,1,1,1,1,2,2,2,2,2,2))

    load(spn1_depletion_viability_rdata)

    temp_panel = ggplot() +
        annotate(geom="text",
                 x=0,
                 y=0,
                 label="Spn1 depletion timecourse?",
                 size=2) +
        labs(tag="a") +
        theme_void() +
        theme(plot.tag=element_text(size=9,
                                    face="bold",
                                    family="FreeSans"),
              plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"))

    figure_spn1_depletion_supp = arrangeGrob(temp_panel,
                                             viability_barplot,
                                             layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_spn1_depletion_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(spn1_depletion_viability_rdata = snakemake@input[["spn1_depletion_viability"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

