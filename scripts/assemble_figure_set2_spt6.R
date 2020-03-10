library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(set2_metagene_rdata,
                spt6_metagene_rdata,
                set2_abundance_chipseq_barplot_rdata,
                spt6_abundance_chipseq_barplot_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(1,1,1,1,1,1,1,1,1,1,1,1),
                   c(2,2,2,2,2,2,2,2,3,3,3,3),
                   c(2,2,2,2,2,2,2,2,3,3,3,3),
                   c(2,2,2,2,2,2,2,2,3,3,3,3),
                   c(2,2,2,2,2,2,2,2,3,3,3,3),
                   c(4,4,4,4,4,4,4,4,5,5,5,5),
                   c(4,4,4,4,4,4,4,4,5,5,5,5),
                   c(4,4,4,4,4,4,4,4,5,5,5,5),
                   c(4,4,4,4,4,4,4,4,5,5,5,5))

    load(set2_metagene_rdata)
    set2_metagene = metagene
    load(spt6_metagene_rdata)
    spt6_metagene = metagene
    load(set2_abundance_chipseq_barplot_rdata)
    set2_abundance_chipseq_barplot = chipseq_abundance_barplot
    load(spt6_abundance_chipseq_barplot_rdata)
    spt6_abundance_chipseq_barplot = chipseq_abundance_barplot

    temp_panel = ggplot() +
        annotate(geom="text",
                 x=0,
                 y=0,
                 label="Co-IP",
                 size=2) +
        labs(tag="a") +
        theme_void() +
        theme(plot.tag=element_text(size=9,
                                    face="bold",
                                    family="FreeSans"),
              plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"))

    figure_set2_spt6 = arrangeGrob(temp_panel,
                                   spt6_metagene,
                                   spt6_abundance_chipseq_barplot,
                                   set2_metagene,
                                   set2_abundance_chipseq_barplot,
                                   layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_set2_spt6,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(set2_metagene_rdata = snakemake@input[["set2_metagene"]],
     spt6_metagene_rdata = snakemake@input[["spt6_metagene"]],
     set2_abundance_chipseq_barplot_rdata = snakemake@input[["set2_abundance_chipseq_barplot"]],
     spt6_abundance_chipseq_barplot_rdata = snakemake@input[["spt6_abundance_chipseq_barplot"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

