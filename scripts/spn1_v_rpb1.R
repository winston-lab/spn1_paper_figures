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

main = function(spn1_path="depleted-v-non-depleted_Spn1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                rpb1_path="depleted-v-non-depleted_Rpb1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "b",
                fig_width=8.5,
                fig_height=9/16*8.5,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = import(spn1_path,
                "Spn1") %>%
        left_join(import(rpb1_path,
                         "Rpb1"),
                  by=c("chrom", "start", "end", "name", "strand", "condition"),
                  suffix=c("_y", "_x")) %>%
        mutate(condition=ordered(condition,
                                 levels=c("control_enrichment",
                                          "condition_enrichment"),
                                 labels=c("non-depleted",
                                          "Spn1-depleted")))

    df_cor = df %>%
        mutate(x=min(enrichment_x, na.rm=TRUE),
               y=max(enrichment_y, na.rm=TRUE)) %>%
        group_by(condition) %>%
        summarize(x=first(x),
                  y=first(y),
                  pearson=cor(enrichment_x,
                              enrichment_y,
                              use="complete.obs"))

    spn1_v_rpb1 = ggplot(data=df,
                         aes(x=enrichment_x,
                             y=enrichment_y)) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=125,
                    shape=16,
                    size=0.2,
                    alpha=0.8) +
        geom_text(data=df_cor,
                  aes(x=x,
                      y=y,
                      label=glue::glue("R={round(pearson, 2)}")),
                  family="FreeSans",
                  size=5/72*25.4,
                  hjust=0,
                  vjust=1) +
        facet_grid(.~condition) +
        scale_color_viridis_c(option="cividis") +
        scale_x_continuous(name="Rpb1 enrichment",
                           breaks=scales::pretty_breaks(3)) +
        scale_y_continuous(name="Spn1\nenrichment",
                           breaks=scales::pretty_breaks(3)) +
        labs(tag=panel_letter) +
        theme_default +
        theme(legend.position="none",
              panel.grid=element_blank(),
              axis.text=element_text(size=5),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5))

    ggsave(pdf_out,
           spn1_v_rpb1,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(spn1_v_rpb1,
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
