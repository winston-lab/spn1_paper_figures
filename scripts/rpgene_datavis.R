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
        group_by(group, assay, annotation, index, position) %>%
        summarize(signal=mean(signal, na.rm=TRUE)) %>%
        group_by(group, assay, annotation, position) %>%
        summarize(low=quantile(signal, 0.25, na.rm=TRUE),
                  mid=median(signal, na.rm=TRUE),
                  high=quantile(signal, 0.75, na.rm=TRUE)) %>%
        bind_rows(df, .) %>%
        return()
}


main = function(theme_path = "spn1_2020_theme.R",
                data_paths = c("RP-genes-plmin-introns-TSS_RNAseq-sense.tsv.gz",
                               "RP-genes-plmin-introns-TSS_GC-pct.tsv.gz",
                               "RP-genes-plmin-introns-TSS_ChIPseq-H3K36me2-H3norm.tsv.gz",
                               "RP-genes-plmin-introns-TSS_ChIPseq-H3K36me3-H3norm.tsv.gz"),
                panel_letter="a",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=17.4/2,
                fig_height=8){
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
               assay=ordered(assay,
                             levels=c("RNAseq-sense",
                                      "GC-pct",
                                      "ChIPseq-Rpb1",
                                      "ChIPseq-Ser5P-Rpb1norm",
                                      "ChIPseq-Ser2P-Rpb1norm",
                                      "ChIPseq-Spn1-Rpb1norm",
                                      "ChIPseq-Spt6-Rpb1norm",
                                      "ChIPseq-Set2-Rpb1norm",
                                      "ChIPseq-H3",
                                      "ChIPseq-H3K4me3-H3norm",
                                      "ChIPseq-H3K36me2-H3norm",
                                      "ChIPseq-H3K36me3-H3norm"),
                             labels=c("\"RNA-seq\"",
                                      "\"GC%\"",
                                      "\"Rpb1: log\"[2] ~ textstyle(frac(\"IP\", \"input\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"Ser5-P\", \"Rpb1\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"Ser2-P\", \"Rpb1\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"Spn1\", \"Rpb1\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"Spt6\", \"Rpb1\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"Set2\", \"Rpb1\"))",
                                      "\"H3: log\"[2] ~ textstyle(frac(\"IP\", \"input\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"H3K4me3\", \"H3\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"H3K36me2\", \"H3\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"H3K36me3\", \"H3\"))")),
               annotation=ordered(annotation,
                                  levels=c("RP genes with introns",
                                           "RP genes without introns"),
                                  labels=c("\"RP genes with introns\"",
                                           "\"RP genes without introns\"")))
    rpgene_datavis = ggplot(data=df %>%
               filter(between(position, -0.1, 1.2)),
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
        facet_grid(assay~annotation,
                   scales="free_y",
                   labeller=label_parsed,
                   switch="y") +
        scale_x_continuous(expand=c(0,0),
                           limits=c(NA, NA),
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==1 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           name=NULL) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              panel.spacing.y=unit(2, "pt"),
              panel.spacing.x=unit(4, "pt"),
              axis.ticks.length=unit(1, "pt"),
              strip.text.x=element_text(margin=margin(b=0, unit="pt")),
              strip.text.y=element_text(angle=-180,
                                        vjust=0.5,
                                        hjust=1,
                                        margin=margin(l=8, r=0, unit="pt")),
              strip.placement="outside",
              legend.title=element_blank(),
              legend.text=element_text(size=5),
              legend.justification=c(0,1),
              # legend.position=c(0.55,0.21),
              legend.position=c(0.04,1.01),
              legend.background=element_blank(),
              legend.key.width=unit(8, "pt"),
              legend.key.height=unit(7, "pt"),
              legend.spacing.x=unit(1, "pt"),
              axis.text.y=element_text(size=5),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"))

    ggsave(pdf_out,
           plot=rpgene_datavis,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(rpgene_datavis,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_paths=snakemake@input[["data"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])

