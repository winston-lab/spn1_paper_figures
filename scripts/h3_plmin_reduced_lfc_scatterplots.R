import = function(df,
                  data_path,
                  chip_factor_id){
    read_tsv(data_path) %>%
        mutate(chip_factor=chip_factor_id) %>%
        bind_rows(df, .) %>%
        return()
}

plot_gene_group_lfc = function(df,
                               x_var=log2FC_enrichment_h3,
                               x_title=quote(atop("H3 enrichment:",
                                                  "log"[2] ~ textstyle(frac("Spn1-depleted",
                                                                            "non-depleted")))),
                               color_var=high_expression,
                               highlight_color="green",
                               legend_title="top 5% non-depleted Rpb1",
                               panel_letter="a"){
    color_var = enquo(color_var)
    x_var = enquo(x_var)

    df %<>%
        filter(log2FC_enrichment > -6 &
                   !is.na(log2FC_enrichment) &
                   !is.na(log2FC_enrichment_h3) &
                   !is.na(h3_affected))

    df_cor = df %>%
        group_by(chip_factor) %>%
        mutate(x_max = max(log2FC_enrichment, na.rm=TRUE),
               y_min = min(!!x_var, na.rm=TRUE)) %>%
        group_by(chip_factor,
                 !!color_var) %>%
        summarize(pearson=cor(log2FC_enrichment,
                              !!x_var,
                              use="complete.obs"),
                  x_max = first(x_max),
                  y_min = first(y_min))

    plot = ggplot(data = df %>%
                      arrange(chip_factor,
                              !!color_var),
                  aes(y=!!x_var,
                      x=log2FC_enrichment)) +
        geom_hline(yintercept=0,
                   color="gray70",
                   size=0.2) +
        geom_vline(xintercept=0,
                   color="gray70",
                   size=0.2) +
        geom_smooth(aes(group=!!color_var,
                        color=!!color_var),
                    method="lm",
                    size=0.2,
                    show.legend=FALSE) +
        geom_point(aes(color=!!color_var),
                   size=0.3,
                   shape=1,
                   stroke=0.1,
                   alpha=0.5) +
        geom_text(data=df_cor,
                  aes(x=x_max,
                      y=y_min+(0.25*!(!!color_var)),
                      label=round(pearson, 2),
                      color=!!color_var),
                  hjust=1,
                  vjust=0,
                  family="FreeSans",
                  size=5/72*25.4,
                  show.legend=FALSE) +
        facet_wrap(~chip_factor,
                   scales="free_x",
                   nrow=1,
                   strip.position="bottom") +
        scale_y_continuous(breaks=scales::pretty_breaks(4),
                           name=x_title) +
        scale_x_continuous(breaks=scales::pretty_breaks(4),
                           name=quote("log"[2] ~ textstyle(frac("Spn1-depleted",
                                                                "non-depleted")))) +
        scale_color_manual(values=c("black", highlight_color),
                           na.value="black",
                           name=legend_title,
                           labels=function(x) str_to_title(x),
                           guide=guide_legend(#label.position="top",
                                              override.aes=list(size=2, stroke=0.4, alpha=1),
                                              keyheight=unit(10, "pt"),
                                              keywidth=unit(12, "pt"))) +
        labs(tag=panel_letter) +
        theme_default +
        theme(axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              panel.grid=element_blank(),
              legend.position="right",
              legend.spacing.y=unit(0, "pt"),
              legend.spacing.x=unit(0, "pt"),
              strip.placement="outside",
              strip.text.x=element_text(margin=margin(t=0, b=4, unit="pt")),
              axis.text.x=element_text(margin=margin(b=-3, unit="pt")),
              plot.margin=margin(11/2, 11/2, 0, 11/2, "pt"))
    return(plot)
}

main = function(data_paths=c("depleted-v-non-depleted_Rpb1-chipseq-libsizenorm-verified-coding-genes-diffbind-results-all.tsv",
                             "depleted-v-non-depleted_Spt6-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv",
                             "depleted-v-non-depleted_Spt6-over-Rpb1-chipseq-spikenorm-verified-coding-genes-diffbind-results-all.tsv"),
                chip_factors=c("Rpb1",
                               "Spt6",
                               "Spt6/Rpb1"),
                h3_path="depleted-v-non-depleted_H3-chipseq-libsizenorm-verified-coding-genes-diffbind-results-all.tsv",
                theme_path = "spn1_2020_theme.R",
                fig_width=17.4,
                fig_height=12/3,
                panel_letter="x",
                pdf_out="test.pdf",
                grob_out="test.Rdata"){
    source(theme_path)

    df = tibble()
    for (i in 1:length(data_paths)){
        df %<>%
            import(data_paths[i],
                   chip_factors[i])
    }

    df %<>%
        left_join(read_tsv(h3_path),
                  by=c("chrom","start","end","name","strand"),
                  suffix=c("", "_h3")) %>%
        mutate(h3_affected = log10_padj_h3 > -log10(0.1) & log2FC_enrichment_h3 < 0,
               chip_factor=ordered(chip_factor,
                                   levels=c("Rpb1",
                                            "Spn1",
                                            "Spn1/Rpb1",
                                            "Spt6",
                                            "Spt6/Rpb1",
                                            "H3K36me3",
                                            "H3K36me3/H3")))

    h3_plmin_reduced_lfc_scatterplots = plot_gene_group_lfc(df=df,
                                                            color_var=h3_affected,
                                                            highlight_color="red",
                                                            legend_title="H3 reduced",
                                                            panel_letter=panel_letter)

    ggsave(pdf_out,
           h3_plmin_reduced_lfc_scatterplots,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(h3_plmin_reduced_lfc_scatterplots,
         file=grob_out)
}

main(data_paths = snakemake@input[["data_paths"]],
     chip_factors = snakemake@params[["chip_factors"]],
     h3_path=snakemake@input[["h3_path"]],
     theme_path=snakemake@input[["theme"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

