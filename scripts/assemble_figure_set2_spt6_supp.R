library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(set2_spt6_v_rpb1_rdata,
                set2_spt6_rpb1norm_v_rpb1_rdata,
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

    load(set2_spt6_v_rpb1_rdata)
    load(set2_spt6_rpb1norm_v_rpb1_rdata)

    figure_set2_spt6_supp = arrangeGrob(set2_spt6_v_rpb1,
                                        set2_spt6_rpb1norm_v_rpb1,
                                        layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_set2_spt6_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(set2_spt6_v_rpb1_rdata = snakemake@input[["set2_spt6_v_rpb1"]],
     set2_spt6_rpb1norm_v_rpb1_rdata = snakemake@input[["set2_spt6_rpb1norm_v_rpb1"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

