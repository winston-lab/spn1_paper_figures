library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(#h3_metagene_rdata,
                h3_modification_datavis_rdata,
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
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1))

    # load(h3_metagene_rdata)
    load(h3_modification_datavis_rdata)

    # temp_panel = ggplot() +
    #     annotate(geom="text",
    #              x=0,
    #              y=0,
    #              label="H3 shift vs Spt6",
    #              size=2) +
    #     labs(tag="c") +
    #     theme_void() +
    #     theme(plot.tag=element_text(size=9,
    #                                 face="bold",
    #                                 family="FreeSans"),
    #           plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"))

    figure_h3_mods = arrangeGrob(#h3_metagene,
                                 #temp_panel,
                                 h3_modification_datavis,
                                 layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_h3_mods,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(#h3_metagene_rdata = snakemake@input[["h3_metagene"]],
     h3_modification_datavis_rdata = snakemake@input[["h3_modification_datavis"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

