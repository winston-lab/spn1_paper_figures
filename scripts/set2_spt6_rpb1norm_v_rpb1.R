main = function(spt6_path="depleted-v-non-depleted_Spt6-over-Rpb1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                set2_path="depleted-v-non-depleted_Set2-over-Rpb1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                rpb1_path="depleted-v-non-depleted_Rpb1-chipseq-libsizenorm-verified-coding-genes-diffbind-results-all.tsv",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "b",
                fig_width=8.5,
                fig_height=7,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = read_tsv(spt6_path) %>%
        # mutate(chip_factor="atop(\"log\"[2] ~ textstyle(frac(\"Spt6\", \"Rpb1\")) * \",\", \"non-depleted\")") %>%
        mutate(chip_factor="atop(\"log\"[2] ~ frac(\"Spt6\", \"Rpb1\") * \",\", \"non-depleted\")") %>%
        bind_rows(read_tsv(set2_path) %>%
                      # mutate(chip_factor="atop(\"log\"[2] ~ textstyle(frac(\"Set2\", \"Rpb1\")) * \",\", \"non-depleted\")")) %>%
                      mutate(chip_factor="atop(\"log\"[2] ~ frac(\"Set2\", \"Rpb1\") * \",\", \"non-depleted\")")) %>%
        mutate(chip_factor = fct_inorder(chip_factor, ordered=TRUE)) %>%
        left_join(read_tsv(rpb1_path),
                  by=c("chrom", "start", "end", "name", "strand"),
                  suffix=c("_factor", "_rpb1"))

    # df_cor = df %>%
    #     group_by(chip_factor) %>%
    #     summarize(pearson=cor(control_enrichment_rpb1,
    #                           control_enrichment_factor,
    #                           use="complete.obs"),
    #               x=min(control_enrichment_rpb1,
    #                     na.rm=TRUE),
    #               y=max(control_enrichment_factor,
    #                     na.rm=TRUE))

    set2_spt6_rpb1norm_v_rpb1 = ggplot(data=df,
           aes(x=control_enrichment_rpb1,
               y=control_enrichment_factor)) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=150,
                    shape=16,
                    size=0.3,
                    alpha=0.8) +
        geom_smooth(method="lm",
                    size=0.2,
                    color=viridisLite::viridis(2, end=0.8)[2]) +
        # geom_text(data=df_cor,
        #           aes(x=x,
        #               y=y,
        #               label=glue::glue("R={round(pearson, 2)}")),
        #           family="FreeSans",
        #           size=5/72*25.4,
        #           hjust=0,
        #           vjust=1) +
        scale_color_viridis_c(option="cividis") +
        scale_x_continuous(name="non-depleted Rpb1 enrichment") +
        facet_grid(chip_factor ~ .,
                   scales="free_y",
                   switch="y",
                   labeller=label_parsed) +
        labs(tag=panel_letter) +
        theme_default +
        theme(legend.position="none",
              panel.grid=element_blank(),
              axis.text=element_text(size=5),
              axis.title.y=element_blank(),
              strip.placement="outside",
              strip.text.y=element_text(angle=-180))

    ggsave(pdf_out,
           set2_spt6_rpb1norm_v_rpb1,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(set2_spt6_rpb1norm_v_rpb1,
         file=grob_out)
}

main(spt6_path=snakemake@input[["spt6"]],
     set2_path=snakemake@input[["set2"]],
     rpb1_path=snakemake@input[["rpb1"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

