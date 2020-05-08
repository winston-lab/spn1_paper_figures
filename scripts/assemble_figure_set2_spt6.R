library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(coip_western_rdata,
                set2_metagene_rdata,
                spt6_metagene_rdata,
                set2_maplot_rdata,
                spt6_maplot_rdata,
                set2_abundance_chipseq_barplot_rdata,
                spt6_abundance_chipseq_barplot_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(2,2,3,3,3,3,3,4,4,4,4,4),
                   c(2,2,3,3,3,3,3,4,4,4,4,4),
                   c(2,2,3,3,3,3,3,4,4,4,4,4),
                   c(2,2,3,3,3,3,3,4,4,4,4,4),
                   c(5,5,6,6,6,6,6,7,7,7,7,7),
                   c(5,5,6,6,6,6,6,7,7,7,7,7),
                   c(5,5,6,6,6,6,6,7,7,7,7,7),
                   c(5,5,6,6,6,6,6,7,7,7,7,7))

    load(coip_western_rdata)
    load(set2_metagene_rdata)
    set2_metagene = metagene
    load(spt6_metagene_rdata)
    spt6_metagene = metagene
    load(set2_maplot_rdata)
    set2_maplot = maplot
    load(spt6_maplot_rdata)
    spt6_maplot = maplot
    load(set2_abundance_chipseq_barplot_rdata)
    set2_abundance_chipseq_barplot = chipseq_abundance_barplot
    load(spt6_abundance_chipseq_barplot_rdata)
    spt6_abundance_chipseq_barplot = chipseq_abundance_barplot

    figure_set2_spt6 = arrangeGrob(coip_western,
                                   spt6_abundance_chipseq_barplot,
                                   spt6_maplot,
                                   spt6_metagene,
                                   set2_abundance_chipseq_barplot,
                                   set2_maplot,
                                   set2_metagene,
                                   layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_set2_spt6,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(coip_western = snakemake@input[["coip_western"]],
     set2_metagene_rdata = snakemake@input[["set2_metagene"]],
     spt6_metagene_rdata = snakemake@input[["spt6_metagene"]],
     set2_maplot_rdata = snakemake@input[["set2_maplot"]],
     spt6_maplot_rdata = snakemake@input[["spt6_maplot"]],
     set2_abundance_chipseq_barplot_rdata = snakemake@input[["set2_abundance_chipseq_barplot"]],
     spt6_abundance_chipseq_barplot_rdata = snakemake@input[["spt6_abundance_chipseq_barplot"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

