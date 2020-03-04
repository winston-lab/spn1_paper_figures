main = function(rnaseq_path="Spn1-IAA-v-Spn1-DMSO_rnaseq-spikenorm-verified-coding-genes-diffexp-results-all.tsv",
                rpb1_path="depleted-v-non-depleted_Rpb1-chipseq-libsizenorm-verified-coding-genes-diffbind-results-all.tsv",
                # nonoverlapping_path="Scer_transcripts_w_verifiedORFs-nonoverlapping.bed",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "x",
                fig_width=8.5,
                fig_height=9/16*8.5,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = read_tsv(rnaseq_path) %>%
        left_join(read_tsv(rpb1_path),
                  by=c("chrom", "start", "end", "name", "strand"),
                  suffix=c("_rna", "_rpb1")) #%>%
        # semi_join(read_tsv(nonoverlapping_path,
        #                    col_names=c("chrom", "start", "end",
        #                                "name", "score", "strand")),
        #           by=c("chrom", "start", "end", "name", "strand"))

    pearson = cor(df$log2_foldchange,
                  df$log2FC_enrichment)

    rnaseq_v_rpb1 = ggplot(data=df,
           aes(y=log2_foldchange,
               x=log2FC_enrichment)) +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_vline(xintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_abline(slope=1,
                    intercept=0,
                    size=0.2,
                    color="gray70")+
        # geom_smooth(method="lm",
                    # size=0.5) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=150,
                    shape=16,
                    size=0.25,
                    alpha=0.8,
                    position=position_jitter(width=0.004,
                                             height=0.004)) +
        annotate(geom="text",
                 x=min(df[["log2FC_enrichment"]]),
                 y=max(df[["log2_foldchange"]]),
                 label=glue::glue("R={signif(pearson, 2)}"),
                 vjust=1,
                 hjust=0,
                 size=5/72*25.4,
                 family="FreeSans") +
        scale_color_viridis_c() +
        scale_x_continuous(breaks=scales::pretty_breaks(4),
                           name=expression("Rpb1 ChIP-seq enrichment: log"[2] ~
                                           textstyle(frac("Spn1-depleted",
                                                          "non-depleted")))) +
        scale_y_continuous(breaks=scales::pretty_breaks(4),
                           name=expression(atop("RNA-seq:",
                                                "log"[2] ~
                                                    textstyle(frac("Spn1-depleted",
                                                                   "non-depleted"))))) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.title.y=element_text(angle=0,
                                        hjust=1,
                                        vjust=0.5))

    ggsave(pdf_out,
           rnaseq_v_rpb1,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(rnaseq_v_rpb1,
         file=grob_out)
}

main(rnaseq_path=snakemake@input[["rnaseq_path"]],
     rpb1_path=snakemake@input[["rpb1_path"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

