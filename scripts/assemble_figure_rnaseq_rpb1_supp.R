library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(extrafont)

main = function(differential_expression_rtqpcr_rdata,
                rna_single_v_custom_rdata,
                esr_fgsea_rdata,
                slow_growth_signature_rdata,
                antisense_single_locus_datavis_rdata,
                chipseq_abundance_barplots_rpb1_rdata,
                fig_width=8.5,
                fig_height=9/16 * 8.5 * 2,
                pdf_out="test.pdf"){
    layout = rbind(c(2,2,2,2,2,2,1,1,1,1,1,1),
                   c(2,2,2,2,2,2,1,1,1,1,1,1),
                   c(2,2,2,2,2,2,1,1,1,1,1,1),
                   c(2,2,2,2,2,2,1,1,1,1,1,1),
                   c(3,3,3,3,3,3,4,4,4,4,4,4),
                   c(3,3,3,3,3,3,4,4,4,4,4,4),
                   c(3,3,3,3,3,3,4,4,4,4,4,4),
                   c(3,3,3,3,3,3,4,4,4,4,4,4),
                   # c(5,5,5,5,5,5,5,5,6,6,6,6),
                   c(5,5,5,5,5,5,5,5,6,6,6,6),
                   c(5,5,5,5,5,5,5,5,6,6,6,6),
                   c(5,5,5,5,5,5,5,5,6,6,6,6))

    load(differential_expression_rtqpcr_rdata)
    load(rna_single_v_custom_rdata)
    load(esr_fgsea_rdata)
    load(slow_growth_signature_rdata)
    load(antisense_single_locus_datavis_rdata)
    load(chipseq_abundance_barplots_rpb1_rdata)

    figure_rnaseq_rpb1_supp = arrangeGrob(diffexp_rtqpcr_scatter,
                                          rna_single_vs_custom,
                                          esr_fgsea,
                                          slow_growth_signature,
                                          rnaseq_single_locus_datavis,
                                          chipseq_abundance_barplot,
                                          layout_matrix=layout)

    ggsave(pdf_out,
           plot=figure_rnaseq_rpb1_supp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
}

main(differential_expression_rtqpcr_rdata = snakemake@input[["differential_expression_rtqpcr"]],
     rna_single_v_custom_rdata = snakemake@input[["rna_single_v_custom"]],
     esr_fgsea_rdata = snakemake@input[["esr_fgsea"]],
     slow_growth_signature_rdata = snakemake@input[["slow_growth_signature"]],
     antisense_single_locus_datavis_rdata = snakemake@input[["antisense_single_locus_datavis"]],
     chipseq_abundance_barplots_rpb1_rdata = snakemake@input[["chipseq_abundance_barplots_rpb1"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]])

