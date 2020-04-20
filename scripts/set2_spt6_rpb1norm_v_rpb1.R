import = function(data_path,
                  chip_factor_id){
    read_tsv(data_path) %>%
        select(chrom, start, end, name, strand, condition_enrichment, control_enrichment) %>%
        pivot_longer(c(condition_enrichment,
                       control_enrichment),
                     names_to="condition",
                     values_to="enrichment") %>%
        mutate(chip_factor=chip_factor_id) %>%
        return()
}

main = function(spt6_path="depleted-v-non-depleted_Spt6-over-Rpb1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                set2_path="depleted-v-non-depleted_Set2-over-Rpb1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                rpb1_path="depleted-v-non-depleted_Rpb1-chipseq-libsizenorm-verified-coding-genes-diffbind-results-all.tsv",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "b",
                fig_width=17.4/2,
                fig_height=7,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = import(spt6_path,
           "\"log\"[2] ~ frac(\"Spt6\", \"Rpb1\")") %>%
        bind_rows(import(set2_path,
                         "\"log\"[2] ~ frac(\"Set2\", \"Rpb1\")")) %>%
        left_join(import(rpb1_path,
                         "Rpb1 enrichment"),
                  by=c("chrom", "start", "end", "name", "strand", "condition"),
                  suffix=c("_y",
                           "_x")) %>%
        mutate(chip_factor_y=fct_inorder(chip_factor_y,
                                         ordered=TRUE),
               condition=ordered(condition,
                                 levels=c("control_enrichment",
                                          "condition_enrichment"),
                                 labels=c("\"non-depleted\"",
                                          "\"Spn1-depleted\"")))

    # df_cor = df %>%
    #     group_by(chip_factor_y) %>%
    #     mutate(x=min(enrichment_x,
    #                  na.rm=TRUE),
    #            y=max(enrichment_y,
    #                  na.rm=TRUE)) %>%
    #     group_by(chip_factor_y,
    #              condition) %>%
    #     summarize(pearson=cor(enrichment_x,
    #                           enrichment_y,
    #                           use="complete.obs"),
    #               x=first(x),
    #               y=first(y))

    set2_spt6_rpb1norm_v_rpb1 = ggplot(data=df %>%
                                           filter(enrichment_y > -6),
                                       aes(x=enrichment_x,
                                           y=enrichment_y)) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=125,
                    shape=16,
                    size=0.2,
                    alpha=0.8) +
        # geom_text(data=df_cor,
        #           aes(x=x,
        #               y=y,
        #               label=glue::glue("R={round(pearson, 2)}")),
        #           family="FreeSans",
        #           size=5/72*25.4,
        #           hjust=0,
        #           vjust=1) +
        # geom_smooth(method="lm",
        #             size=0.2,
        #             color=viridisLite::viridis(2, end=0.8)[2]) +
        scale_color_viridis_c(option="cividis") +
        scale_x_continuous(name="Rpb1 enrichment") +
        facet_grid(chip_factor_y ~ condition,
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

