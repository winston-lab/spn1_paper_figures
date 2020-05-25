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
                  n=first(n)) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(theme_path = "spn1_2020_theme.R",
                data_paths = c("reduced-H3-expression-length-match_ChIPseq-H3.tsv.gz",
                               "reduced-H3-expression-length-match_ChIPseq-Rpb1.tsv.gz",
                               "reduced-H3-expression-length-match_ChIPseq-Spt6.tsv.gz"),
                panel_letter="b",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=11.4,
                fig_height=11.4 * 9 / 16){
    source(theme_path)

    df = tibble()
    for (path in data_paths){
        df %<>% import(path)
    }

    df %<>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               # annotation_numbered = paste0("\"", n, " ", annotation, "\""),
               annotation_numbered = ordered(annotation,
                                             levels=c("genes with reduced H3",
                                                      "matched genes without reduced H3"),
                                             labels=c("\"77 genes with reduced H3\"",
                                                      "\"80 matched genes\nwithout reduced H3\"")),
               assay=ordered(assay,
                             levels=c("ChIPseq-H3",
                                      # "ChIPseq-Rpb1",
                                      "ChIPseq-Rpb1-spikenorm",
                                      "ChIPseq-Spt6",
                                      # "ChIPseq-Spt6-Rpb1norm"),
                                      "ChIPseq-Spt6-Rpb1norm-batchAB"),
                             labels=c("\"H3: log\"[2] ~ textstyle(frac(\"IP\", \"input\"))",
                                      "\"Rpb1: log\"[2] ~ textstyle(frac(\"IP\", \"input\"))",
                                      "\"Spt6: log\"[2] ~ textstyle(frac(\"IP\", \"input\"))",
                                      # "\"log\"[2] ~ frac(\"Spt6\", \"Rpb1\")")))
                                      "\"log\"[2] ~ textstyle(frac(\"Spt6\", \"Rpb1\"))")))

    reduced_h3_matched_metagenes = ggplot(data=df,
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
        geom_line(size=0.5,
                  alpha=0.8) +
        facet_grid(assay ~ annotation_numbered,
                   scales="free_y",
                   switch="y",
                   labeller=label_parsed) +
        scale_x_continuous(expand=c(0,0),
                           limits=c(NA, NA),
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==2 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           name=NULL) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              panel.spacing=unit(2, "pt"),
              strip.text.x=element_text(margin=margin(t=8, b=1, unit="pt")),
              strip.text.y=element_text(angle=-180,
                                        hjust=1,
                                        margin=margin(l=6, unit="pt")),
              strip.placement="outside",
              legend.key.width=unit(10, "pt"),
              legend.title=element_blank(),
              legend.justification=c(0.5,0.5),
              legend.position=c(0.33, 0.86),
              legend.background=element_blank(),
              legend.spacing.x=unit(1, "pt"),
              axis.text.y=element_text(size=5),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              plot.title=element_text(hjust=0.5),
              plot.margin=margin(11/2, 11/2 + 1, 11/2, 11/2, "pt"))

    ggsave(pdf_out,
           plot=reduced_h3_matched_metagenes,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(reduced_h3_matched_metagenes,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_paths=snakemake@input[["data"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])

