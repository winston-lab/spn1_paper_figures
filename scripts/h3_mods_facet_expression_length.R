import = function(df,
                  data_path){
    read_tsv(data_path,
             col_names=c("group",
                         "sample",
                         "annotation",
                         "assay",
                         "index",
                         "position",
                         "signal")) %>%
        group_by(group, annotation, assay, index, position) %>%
        summarize(signal=mean(signal, na.rm=TRUE)) %>%
        group_by(group, annotation, assay) %>%
        mutate(n = n_distinct(index)) %>%
        group_by(group, annotation, assay, position) %>%
        summarize(low=quantile(signal, 0.25, na.rm=TRUE),
                  mid=median(signal, na.rm=TRUE),
                  high=quantile(signal, 0.75, na.rm=TRUE),
                  n=first(n),
                  n_position=n_distinct(index)) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(theme_path = "spn1_2020_theme.R",
                data_paths = c("verified-transcripts-nonoverlapping-TSS-facet-expression-length_ChIPseq-H3K4me3-H3norm.tsv.gz"),
                modification_id="H3K4me3",
                panel_letter="b",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=17.4/2,
                fig_height=4){
    source(theme_path)

    df = tibble()
    for (path in data_paths){
        df %<>% import(path)
    }

    df %<>%
        ungroup() %>%
        separate(annotation,
                 into=c("expression_class",
                        "length_class"),
                 sep="_") %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               expression_class=ordered(expression_class,
                                        levels=c("low",
                                                 "mid",
                                                 "high")),
               length_class=ordered(length_class,
                                    levels=c("short",
                                             "medium",
                                             "long")))

    h3_mods_facet_expression_length = ggplot(data=df %>%
                      filter(n_position > 20),
           aes(x=position,
               y=mid,
               ymin=low,
               ymax=high,
               color=group,
               fill=group)) +
        geom_vline(xintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_ribbon(linetype="blank",
                    alpha=0.1) +
        geom_line(size=0.25,
                  alpha=0.8) +
        facet_grid(expression_class ~ length_class) +
        scale_x_continuous(expand=c(0,0),
                           limits=c(NA, NA),
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==2 ~ paste(x, "kb"),
                                                        x>2 ~ "",
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           name=bquote("log"[2] ~ textstyle(frac(.(modification_id), "H3")))) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              panel.spacing=unit(1, "pt"),
              axis.ticks.length=unit(1, "pt"),
              strip.text.y=element_text(angle=0,
                                        hjust=0,
                                        margin=margin(l=1, r=1, unit="pt")),
              strip.text.x=element_text(margin=margin(b=1, unit="pt")),
              legend.title=element_blank(),
              legend.justification=c(0.5,0.5),
              legend.key.size=unit(6, "pt"),
              legend.position=c(-0.16, 0.85),
              legend.background=element_blank(),
              legend.text=element_text(size=5),
              legend.spacing.x=unit(1, "pt"),
              axis.text=element_text(size=5),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              plot.title=element_text(hjust=0.5),
              plot.margin=margin(0, 0, 0, 11/2, "pt"))

    ggsave(pdf_out,
           plot=h3_mods_facet_expression_length,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(h3_mods_facet_expression_length,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_paths=snakemake@input[["data"]],
     modification_id=snakemake@wildcards[["mod"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])

