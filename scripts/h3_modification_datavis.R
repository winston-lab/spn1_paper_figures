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
        group_by(group, assay, index, position) %>%
        summarize(signal=mean(signal, na.rm=TRUE)) %>%
        bind_rows(df, .) %>%
        return()
}

plot_heatmap = function(df,
                        filter_assay,
                        quantile_low=0.05,
                        quantile_high=0.95,
                        colorbar_title=expression("relative enrichment," ~
                                                      textstyle(frac("H3K4me3",
                                                                     "H3"))),
                        leftmost=TRUE){
    df %<>%
        filter(assay==filter_assay)

    heatmap = ggplot() +
        geom_raster(data=df,
                     aes(x=position,
                         y=sorted_index,
                         fill=signal)) +
        geom_text(data=df %>%
                      group_by(group) %>%
                      summarize(sorted_index = max(sorted_index)),
                  aes(y=sorted_index,
                      label=group),
                  x=max(df[["position"]]) - 0.5,
                  size=ifelse(leftmost,
                              8/72*25.4,
                              0),
                  hjust=1,
                  nudge_y=-500,
                  family="FreeSans") +
        facet_grid(group ~ .) +
        scale_fill_viridis_c(limits=quantile(df[["signal"]],
                                             probs=c(quantile_low,
                                                     quantile_high),
                                             na.rm=TRUE),
                             oob=scales::squish,
                             option="viridis",
                             na.value=viridisLite::viridis(1),
                             name=colorbar_title,
                             guide=guide_colorbar(barwidth=unit(4, "cm"),
                                                  barheight=unit(0.2, "cm"),
                                                  title.position="top",
                                                  title.hjust=0.5,
                                                  breaks=scales::pretty_breaks(3))) +
        scale_x_continuous(expand=c(0,0),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==3 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x))) +
        scale_y_continuous(expand=c(0,25),
                           name="3087 non-overlapping coding genes") +
        theme_default +
        theme(axis.text.y=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_text(color=ifelse(leftmost,
                                                     "black",
                                                     "white"),
                                        margin=margin(l=0, r=4, unit="pt")),
              panel.grid=element_blank(),
              panel.grid.major.x=element_line(size=0.1,
                                              color="gray85"),
              panel.border=element_blank(),
              legend.position="top",
              legend.margin=margin(b=-8, unit="pt"),
              legend.title=element_text(margin=margin(b=-4, unit="pt")),
              strip.text=element_blank(),
              panel.spacing.y=unit(2, "pt"),
              axis.ticks.length.y=unit(1, "pt"),
              plot.margin=margin(11/2, 6, 0, -3, "pt"))

    return(heatmap)

    # ggsave("test.pdf",
    #        plot=heatmap,
    #        width=17.4/3,
    #        height=9/16 * 17.4,
    #        units="cm",
    #        device=cairo_pdf)
}

plot_metagene = function(df_metagene,
                         filter_assay,
                         plot_title=expression(textstyle(frac("H3K36me3",
                                                              "H3")) ~
                                                   "ChIP-seq enrichment"),
                         # yaxis_label = expression("log"[2] ~
                         #                              textstyle(frac("H3K4me3",
                         #                                             "H3"))),
                         legend_position="none",
                         leftmost=FALSE){
    df_metagene %<>%
        filter(assay==filter_assay)

    metagene = ggplot(data=df_metagene,
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
        scale_x_continuous(expand=c(0,0),
                           limits=c(NA, NA),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==3 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           # name=yaxis_label) +
                           name=expression("log"[2] ~ "ratio")) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        ggtitle(plot_title) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.title=element_blank(),
              legend.justification=c(1,1),
              legend.position=legend_position,
              legend.background=element_blank(),
              legend.spacing.x=unit(1, "pt"),
              axis.text.y=element_text(size=5),
              axis.title.y=element_text(color=ifelse(leftmost,
                                                     "black",
                                                     "white"),
                                        margin=margin(l=0, r=-3, unit="pt")),
              plot.margin=margin(11/2, 6, 11/2, -3, "pt"),
              plot.title=element_text(hjust=0.5))

    return(metagene)
}


main = function(theme_path = "spn1_2020_theme.R",
                data_paths = c("verified-transcripts-nonoverlapping-TSS_ChIPseq-H3K4me3-H3norm.tsv.gz",
                               "verified-transcripts-nonoverlapping-TSS_ChIPseq-H3K36me2-H3norm.tsv.gz",
                               "verified-transcripts-nonoverlapping-TSS_ChIPseq-H3K36me3-H3norm.tsv.gz"),
                annotation_path = "Scer_transcripts_w_verifiedORFs-nonoverlapping.bed",
                panel_letter="c",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=17.4,
                fig_height=17.4 * 9 / 16){
    source(theme_path)
    library(cowplot)
    library(ggplotify)

    df = tibble()
    for (path in data_paths){
        df %<>% import(path)
    }

    annotation = read_tsv(annotation_path,
                          col_names=c("chrom", "start", "end",
                                      "name", "score", "strand")) %>%
        mutate(index = row_number()) %>%
        arrange(desc(end - start)) %>%
        mutate(sorted_index = row_number())

    df %<>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               assay=ordered(assay,
                             levels=c("ChIPseq-H3K4me3-H3norm",
                                      "ChIPseq-H3K36me2-H3norm",
                                      "ChIPseq-H3K36me3-H3norm"),
                             labels=c("textstyle(frac(\"H3K4me3\",\"H3\"))",
                                      "textstyle(frac(\"H3K36me2\",\"H3\"))",
                                      "textstyle(frac(\"H3K36me3\",\"H3\"))"))) %>%
        left_join(annotation,
                  by="index")

    df_metagene = df %>%
        group_by(group, assay, position) %>%
        summarize(low=quantile(signal, 0.25, na.rm=TRUE),
                  mid=median(signal, na.rm=TRUE),
                  high=quantile(signal, 0.75, na.rm=TRUE))


    assays = levels(df[["assay"]])

    plots = list(plot_metagene(df=df_metagene,
                               filter_assay=assays[1],
                               # yaxis_label=expression("log"[2] ~ textstyle(frac("H3K4me3",
                               #                                                  "H3"))),
                               # plot_title=expression(textstyle(frac("H3K4me3",
                               #                                      "H3"))),
                               plot_title="H3K4me3 / H3",
                               legend_position=c(0.95, 0.95),
                               leftmost=TRUE),
                 plot_metagene(df=df_metagene,
                               filter_assay=assays[2],
                               # yaxis_label=expression("log"[2] ~ textstyle(frac("H3K36me2", "H3")))),
                               # plot_title=expression(textstyle(frac("H3K36me2",
                               #                                      "H3")))),
                               plot_title="H3K36me2 / H3"),
                 plot_metagene(df=df_metagene,
                               filter_assay=assays[3],
                               # yaxis_label=expression("log"[2] ~ textstyle(frac("H3K36me3", "H3")))))
                               # plot_title=expression(textstyle(frac("H3K36me3",
                               #                                      "H3")))),
                               plot_title="H3K36me3 / H3"),
                 plot_heatmap(df=df,
                              filter_assay=assays[1],
                              quantile_low=0.10,
                              quantile_high=0.90,
                              colorbar_title=expression("log"[2] ~ textstyle(frac("H3K4me3", "H3")))),
                 plot_heatmap(df=df,
                              filter_assay=assays[2],
                              quantile_low=0.10,
                              quantile_high=0.90,
                              colorbar_title=expression("log"[2] ~ textstyle(frac("H3K36me2", "H3"))),
                              leftmost=FALSE),
                 plot_heatmap(df=df,
                              filter_assay=assays[3],
                              quantile_low=0.10,
                              quantile_high=0.90,
                              colorbar_title=expression("log"[2] ~ textstyle(frac("H3K36me3", "H3"))),
                              leftmost=FALSE))

    h3_modification_datavis = plot_grid(plotlist=plots,
                                        align="v",
                                        axis="rl",
                                        nrow=2,
                                        ncol=3,
                                        rel_heights=c(0.4, 1)) %>%
        as.ggplot()

    h3_modification_datavis = h3_modification_datavis +
        labs(tag=panel_letter) +
        theme(plot.tag=element_text(family="FreeSans",
                                    size=9,
                                    face="bold"),
              plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"))

    ggsave(pdf_out,
           plot=h3_modification_datavis,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(h3_modification_datavis,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_paths=snakemake@input[["data"]],
     annotation_path=snakemake@input[["annotation"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])

