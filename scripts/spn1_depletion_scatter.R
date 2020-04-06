main = function(data_path="depleted-v-non-depleted_Spn1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "a",
                fig_width=8.5,
                fig_height=9/16*8.5,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = read_tsv(data_path)

    axis_limits=range(c(df[["control_enrichment"]],
                   df[["condition_enrichment"]]),
                 na.rm=TRUE)

    spn1_depletion_scatter = ggplot(data=df,
           aes(x=control_enrichment,
               y=condition_enrichment)) +
        geom_abline(slope=1,
                    intercept=0,
                    size=0.2,
                    color="gray70") +
        stat_binhex(aes(color=..count..),
                    bins=150,
                    geom="point",
                    shape=16,
                    size=0.3,
                    alpha=0.8) +
        scale_x_continuous(limits=axis_limits,
                           name="non-depleted",
                           breaks=scales::pretty_breaks(3)) +
        scale_y_continuous(limits=axis_limits,
                           name="Spn1-depleted",
                           breaks=scales::pretty_breaks(3)) +
        scale_color_viridis_c(option="cividis") +
        labs(tag=panel_letter,
             title="Spn1 ChIP enrichment") +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              axis.text=element_text(size=5))

    ggsave(pdf_out,
           spn1_depletion_scatter,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(spn1_depletion_scatter,
         file=grob_out)
}

main(data_path=snakemake@input[["data"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

