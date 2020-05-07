main = function(spn1_path="depleted-v-non-depleted_Spn1-over-Rpb1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                rpb1_path="depleted-v-non-depleted_Rpb1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "b",
                fig_width=8.5,
                fig_height=9/16*8.5,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = read_tsv(spn1_path) %>%
        left_join(read_tsv(rpb1_path),
                  by=c("chrom", "start", "end", "name", "strand"),
                  suffix=c("_spn1", "_rpb1"))

    spn1_rpb1norm_v_rpb1_nondepleted = ggplot(data=df,
           aes(x=control_enrichment_rpb1,
               y=control_enrichment_spn1)) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=150,
                    shape=16,
                    size=0.3,
                    alpha=0.8) +
        geom_smooth(method="lm",
                    size=0.2,
                    color=viridisLite::viridis(2, end=0.8)[2]) +
        scale_color_viridis_c(option="cividis") +
        scale_x_continuous(name="non-depleted Rpb1 enrichment") +
        scale_y_continuous(name=expression(atop("log"[2] ~ textstyle(frac("Spn1", "Rpb1") * ","), "non-depleted"))) +
        labs(tag=panel_letter) +
        theme_default +
        theme(legend.position="none",
              panel.grid=element_blank(),
              axis.text=element_text(size=5),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5))


    ggsave(pdf_out,
           spn1_rpb1norm_v_rpb1_nondepleted,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(spn1_rpb1norm_v_rpb1_nondepleted,
         file=grob_out)
}

main(spn1_path=snakemake@input[["spn1"]],
     rpb1_path=snakemake@input[["rpb1"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])
