library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(spn1_depletion_viability_rdata,
                spn1_v_rpb1_rdata,
                spn1_rpb1norm_v_rpb1_nondepleted_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,NA,NA,NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,NA,NA,NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,NA,NA,NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,NA,NA,NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,NA,NA,NA,NA,NA,NA,NA),
                   c(2,2,2,2,2,2,2,3,3,3,3,3),
                   c(2,2,2,2,2,2,2,3,3,3,3,3),
                   c(2,2,2,2,2,2,2,3,3,3,3,3),
                   c(2,2,2,2,2,2,2,3,3,3,3,3),
                   c(2,2,2,2,2,2,2,3,3,3,3,3),
                   c(2,2,2,2,2,2,2,3,3,3,3,3),
                   c(2,2,2,2,2,2,2,3,3,3,3,3))

    load(spn1_depletion_viability_rdata)
    load(spn1_v_rpb1_rdata)
    load(spn1_rpb1norm_v_rpb1_nondepleted_rdata)

    figure_spn1_depletion_supp = arrangeGrob(viability_barplot,
                                             spn1_v_rpb1,
                                             spn1_rpb1norm_v_rpb1_nondepleted,
                                             layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_spn1_depletion_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(spn1_depletion_viability_rdata = snakemake@input[["spn1_depletion_viability"]],
     spn1_v_rpb1_rdata = snakemake@input[["spn1_v_rpb1"]],
     spn1_rpb1norm_v_rpb1_nondepleted_rdata = snakemake@input[["spn1_rpb1norm_v_rpb1_nondepleted"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

