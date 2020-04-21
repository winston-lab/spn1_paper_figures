import = function(h3_path,
                  rpb1_path){
    read_tsv(h3_path) %>%
        left_join(read_tsv(rpb1_path),
                  by=c("chrom", "name", "strand"),
                  suffix=c("_h3", "_rpb1")) %>%
        return()
}

plot_ma = function(df,
                   panel_letter,
                   y_title=expression(atop("H3 ChIP-seq\n(first 500bp):",
                                                "log"[2] ~
                                                    textstyle(frac("Spn1-depleted",
                                                                   "non-depleted"))))){
    ma_plot = ggplot(data=df,
                     aes(y=log2FC_enrichment_h3,
                         x=control_enrichment_rpb1)) +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_point(aes(color=reduced_h3,
                       alpha=reduced_h3),
                   size=0.5,
                   shape=1,
                   stroke=0.1) +
        geom_smooth(size=0.2,
                    # color=viridisLite::viridis(2, end=0.8)[2]) +
                    color="#feb24c") +
        scale_x_continuous(breaks=scales::pretty_breaks(4),
                           name=expression("non-depleted Rpb1 enrichment")) +
        scale_y_continuous(breaks=scales::pretty_breaks(4),
                           name=y_title) +
        scale_color_manual(values=c("black", "#e41a1c")) +
        scale_alpha_manual(values=c(0.4, 0.8)) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        hjust=1,
                                        vjust=0.5),
              legend.position="none")
    return(ma_plot)
}

main = function(h3_path="depleted-v-non-depleted_H3-chipseq-libsizenorm-5p500bp-diffbind-results-all.tsv",
                rpb1_path="depleted-v-non-depleted_Rpb1-chipseq-libsizenorm-verified-coding-genes-diffbind-results-all.tsv",
                theme_path = "spn1_2020_theme.R",
                fig_width=8.5,
                fig_height=9/16*8.5,
                panel_letter="a",
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = import(h3_path,
                rpb1_path)

    df %<>%
        mutate(reduced_h3 = (log10_padj_h3 > -log10(0.1)) & log2FC_enrichment_h3 < 0) %>%
        arrange(reduced_h3)

    h3_vs_rpb1_ma_5p500bp = plot_ma(df,
                                    panel_letter)

    ggsave(pdf_out,
           h3_vs_rpb1_ma_5p500bp,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(h3_vs_rpb1_ma_5p500bp,
         file=grob_out)
}

main(h3_path=snakemake@input[["h3_path"]],
     rpb1_path=snakemake@input[["rpb1_path"]],
     theme_path=snakemake@input[["theme"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

