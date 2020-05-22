main = function(single_path="Spn1-IAA-v-Spn1-DMSO_rnaseq-spikenorm-verified-coding-genes-diffexp-results-all.tsv",
                custom_path="custom-diffexp-spikenorm-verified-coding-genes-results-all.tsv",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "x",
                fig_width=8.5,
                fig_height=9/16*8.5,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = read_tsv(single_path) %>%
        inner_join(read_tsv(custom_path),
                  by=c("chrom", "start", "end", "name", "strand"),
                  suffix=c("_single", "_custom"))

    df_cor = df %>%
        summarize(x=min(log2_foldchange_single, na.rm=TRUE),
                  y=max(log2_foldchange_custom, na.rm=TRUE),
                  pearson=cor(log2_foldchange_single,
                              log2_foldchange_custom,
                              use="complete.obs"))

    rna_single_vs_custom = ggplot(data=df,
           aes(y=log2_foldchange_single,
               x=log2_foldchange_custom)) +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_vline(xintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_text(data=df_cor,
                  aes(x=x,
                      y=y,
                      label=glue::glue("R={round(pearson, 2)}")),
                  hjust=0,
                  vjust=1,
                  size=5/72*25.4,
                  family="FreeSans")+
        # geom_smooth(size=0.5) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=150,
                    shape=16,
                    size=0.25,
                    alpha=0.8,
                    # position=position_jitter(width=0.004,
                                             # height=0.004)
                    ) +
        geom_abline(slope=1,
                    intercept=0,
                    size=0.2,
                    color="gray70",
                    alpha=0.5)+
        scale_color_viridis_c(option="cividis") +
        scale_x_continuous(breaks=scales::pretty_breaks(4),
                           name=quote("log"[2] ~ textstyle(frac("Spn1-depleted",
                                                                         "non-depleted")))) +
        scale_y_continuous(breaks=scales::pretty_breaks(4),
                           name=quote(atop( "log"[2] ~ textstyle(frac("Spn1-depleted",
                                                                      "non-depleted")), "(adjusted for all controls)"))) +
        labs(tag=panel_letter,
             title="RNA-seq") +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.title.y=element_text(angle=0,
                                        hjust=1,
                                        vjust=0.5))

    ggsave(pdf_out,
           rna_single_vs_custom,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(rna_single_vs_custom,
         file=grob_out)
}

main(single_path=snakemake@input[["single"]],
     custom_path=snakemake@input[["custom"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])
