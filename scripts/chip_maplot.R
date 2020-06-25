main = function(theme_path = "spn1_2020_theme.R",
                data_path = "depleted-v-non-depleted_Set2-over-Rpb1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                rpb1_path = "depleted-v-non-depleted_Rpb1-chipseq-libsizenorm-verified-coding-genes-diffbind-results-all.tsv",
                rpb1_norm = FALSE,
                factor_id = "Set2",
                panel_letter = "x",
                fig_width = 8.5,
                fig_height = 8.5 * 9 / 16,
                pdf_out = "test.pdf",
                grob_out="test.Rdata"){
    source(theme_path)

    df = read_tsv(data_path) %>%
        left_join(read_tsv(rpb1_path),
                  by=c("chrom", "start", "end", "name", "strand"),
                  suffix=c("_factor", "_rpb1"))

    y_text = ifelse(rpb1_norm, "/ Rpb1:", "ChIP-seq:")
    y_text = paste(factor_id, y_text)

    maplot = ggplot(data=df,
           aes(x=control_enrichment_rpb1,
               y=log2FC_enrichment_factor)) +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=175,
                    shape=16,
                    size=0.3,
                    alpha=0.8) +
        scale_color_viridis_c(option="cividis") +
        scale_x_continuous(name="non-depleted Rpb1 enrichment",
                           breaks=scales::pretty_breaks(3)) +
        scale_y_continuous(name=bquote(atop(.(y_text),
                                            "log"[2] ~ textstyle(frac("Spn1-depleted",
                                                                      "non-depleted")))),
                           breaks=scales::pretty_breaks(4)) +
        labs(tag=panel_letter) +
        theme_default +
        theme(legend.position="none",
              panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5))

    ggsave(pdf_out,
           plot=maplot,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(maplot,
         file=grob_out)
}

main(theme_path = snakemake@input[["theme"]],
     data_path = snakemake@input[["data"]],
     rpb1_path = snakemake@input[["rpb1"]],
     rpb1_norm = snakemake@params[["rpb1_norm"]],
     factor_id = snakemake@params[["factor_id"]],
     panel_letter = snakemake@params[["panel_letter"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]],
     grob_out = snakemake@output[["grob"]])

