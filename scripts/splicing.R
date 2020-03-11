main = function(data_path="Spn1-IAA-v-all-controls_intron_retention_results.tsv",
                aliases_path = "Scer-sys-to-common-name.tsv",
                rp_genes_path = "Scer_RPgene_transcripts.bed",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "b",
                fig_width=8.5,
                fig_height=6,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){
    source(theme_path)
    library(ggrepel)

    rp_genes = read_tsv(rp_genes_path,
             col_names=c("chrom", "start", "end",
                         "name", "score", "strand"))

    df = read_tsv(data_path) %>%
        separate(name,
                 into=c("systematic", "junk", "number"),
                 sep="_") %>%
        left_join(read_tsv(aliases_path,
                           col_names=c("systematic", "short")),
             by="systematic") %>%
        mutate(name = if_else(is.na(short),
                              systematic,
                              short),
               gene_class = if_else(name %in% rp_genes[["name"]],
                                    "ribosomal protein genes",
                                    "other genes"),
               gene_class = ordered(gene_class,
                                    levels=c("ribosomal protein genes",
                                             "other genes")))

    splicing_volcano = ggplot() +
        geom_vline(xintercept=0,
                   color="grey70",
                   size=0.2) +
        geom_text_repel(data=df,
                        aes(x=ir_est_condition - ir_est_controls,
                            y=ifelse(ir_est_condition >= ir_est_controls,
                                     -log10(qvalue_increase),
                                     -log10(qvalue_decrease)),
                            label=if_else(name %in% c("RPS27A",
                                                      "RPL34A",
                                                      "RPL14B",
                                                      "RPS21B",
                                                      "RPL43B",
                                                      "SFT1"),
                                          name,
                                          "")),
                        size=5/72*25.4,
                        family="FreeSans",
                        fontface="italic",
                        force=0.5,
                        box.padding=unit(3, "pt"),
                        point.padding=unit(1, "pt"),
                        min.segment.length=unit(0, "pt"),
                        segment.size=0.3,
                        segment.color="gray50",
                        hjust=0.90) +
        geom_point(data=df,
                   aes(x=ir_est_condition - ir_est_controls,
                       y=ifelse(ir_est_condition >= ir_est_controls,
                                -log10(qvalue_increase),
                                -log10(qvalue_decrease)),
                       color=gene_class,
                       alpha=(category=="not significant")),
                   shape=16,
                   size=0.8) +
        xlab(expression("(intron retention)"["Spn1-depleted"] -
                            "(intron retention)"["controls"])) +
        scale_y_continuous(name = expression("-log"[10]("FDR")),
                           limits = function(x) c(0, x[2] * 1.05),
                           expand=c(0,0)) +
        scale_alpha_manual(values=c(1, 0.55),
                           guide=FALSE) +
        scale_color_brewer(palette="Set2",
                           name=NULL,
                           guide=guide_legend(override.aes=list(alpha=1),
                                              keywidth=unit(8, "pt"),
                                              label.vjust=0.67)) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              legend.justification=c(0,1),
              legend.position=c(0.003,0.99),
              legend.spacing.x=unit(-1, "pt"),
              legend.background=element_blank(),
              legend.text=element_text(size=5.5))

    ggsave(pdf_out,
           splicing_volcano,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(splicing_volcano,
         file=grob_out)
}

main(data_path=snakemake@input[["data"]],
     aliases_path=snakemake@input[["aliases"]],
     rp_genes_path=snakemake@input[["rp_genes"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

