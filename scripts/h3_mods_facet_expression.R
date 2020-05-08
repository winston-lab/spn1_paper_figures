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

plot_metagene = function(df,
                         filter_assay,
                         leftmost=FALSE){
    df_temp = df %>%
        filter(assay==filter_assay &
                   n_position > 20)
    if (leftmost){
        legend_position = c(0.63, 0.76)
    } else {
        legend_position = "none"
    }

    metagene = ggplot(data=df_temp,
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
        facet_grid(annotation ~ assay,
                   labeller=label_parsed) +
        scale_x_continuous(expand=c(0,0),
                           limits=c(NA, NA),
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==3 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           name=NULL) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        theme_default +
        theme(panel.grid=element_blank(),
              panel.spacing.y=unit(2, "pt"),
              strip.text.y=element_blank(),
              strip.text.x=element_text(margin=margin(t=2, b=0, unit="pt")),
              legend.key.width=unit(10, "pt"),
              legend.key.height=unit(8, "pt"),
              legend.title=element_blank(),
              legend.justification=c(0.5,0.5),
              legend.position=legend_position,
              legend.background=element_blank(),
              legend.spacing.x=unit(1, "pt"),
              axis.text.y=element_text(size=5),
              plot.margin=margin(0,6,0,0,"pt"))
    return(metagene)
}

main = function(theme_path = "spn1_2020_theme.R",
                data_paths = c("verified-transcripts-nonoverlapping-TSS-facet-expression_ChIPseq-H3K36me2-H3norm.tsv.gz",
                               "verified-transcripts-nonoverlapping-TSS-facet-expression_ChIPseq-H3K36me3-H3norm.tsv.gz",
                               "verified-transcripts-nonoverlapping-TSS-facet-expression_ChIPseq-H3K4me3-H3norm.tsv.gz"),
                panel_letter="b",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=17.4 * 9 / 12,
                fig_height=12 * 6 / 12){
    source(theme_path)
    library(cowplot)
    library(ggplotify)

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
               annotation=ordered(annotation,
                                  levels=c("high",
                                           "mid",
                                           "low"),
                                  labels=c("top 1/3\nexpression",
                                           "middle 1/3\nexpression",
                                           "bottom 1/3\nexpression")),
               assay=ordered(assay,
                             levels=c("ChIPseq-H3K36me2-H3norm",
                                      "ChIPseq-H3K36me3-H3norm",
                                      "ChIPseq-H3K4me3-H3norm"),
                             labels=c("\"log\"[2] ~ textstyle(frac(\"H3K36me2\",\"H3\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"H3K36me3\",\"H3\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"H3K4me3\",\"H3\"))")))

    annotation_plot = ggplot(data=distinct(df, annotation)) +
        geom_text(aes(label=annotation,
                      y=0),
                  x=0,
                  hjust=0,
                  size=8/72*25.4,
                  family="FreeSans") +
        facet_grid(annotation ~ .) +
        theme_void() +
        theme(strip.text=element_blank())

    h3_mods_facet_expression = plot_grid(plot_metagene(df,
                                                       "\"log\"[2] ~ textstyle(frac(\"H3K36me2\",\"H3\"))",
                                                       leftmost=TRUE),
                                         plot_metagene(df,
                                                       "\"log\"[2] ~ textstyle(frac(\"H3K36me3\",\"H3\"))"),
                                         plot_metagene(df,
                                                       "\"log\"[2] ~ textstyle(frac(\"H3K4me3\",\"H3\"))"),
                                         annotation_plot,
                                         nrow=1,
                                         rel_widths=c(1,1,1,0.5),
                                         align="h",
                                         axis="tb") %>%
        as.ggplot()

    h3_mods_facet_expression = h3_mods_facet_expression +
        labs(tag=panel_letter) +
        theme(plot.margin=margin(11/2, 0, 0, 11/2, "pt"),
              plot.tag=element_text(size=9,
                                    face="bold"))

    ggsave(pdf_out,
           plot=h3_mods_facet_expression,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(h3_mods_facet_expression,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_paths=snakemake@input[["data"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])

