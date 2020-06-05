import_for_heatmap = function(df,
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

import_for_metagene = function(df,
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
        group_by(group, assay, position) %>%
        summarize(low=quantile(signal, 0.25, na.rm=TRUE),
                  mid=median(signal, na.rm=TRUE),
                  high=quantile(signal, 0.75, na.rm=TRUE)) %>%
        bind_rows(df, .) %>%
        return()
}

plot_heatmap = function(df,
                        filter_assay,
                        colorbar_title=expression("relative enrichment," ~
                                                      textstyle(frac("H3K4me3",
                                                                     "H3"))),
                        leftmost=FALSE){
    df %<>%
        filter(assay==filter_assay) %>%
        pivot_wider(names_from=group,
                    values_from=signal) %>%
        mutate(signal = `Spn1-depleted` - `non-depleted`)

    df_anno = distinct(df,
                       chrom, start, end, name, score, strand, sorted_index) %>%
        mutate(cps = (end - start - 300) / 1000) %>%
        arrange(desc(sorted_index))

    heatmap = ggplot() +
        geom_raster(data=df,
                     aes(x=position,
                         y=sorted_index,
                         fill=signal)) +
        geom_path(data=filter(df_anno,
                              cps <= 3),
                  aes(x=cps,
                      y=sorted_index),
                  size=0.3,
                  linetype="dotted",
                  alpha=1,
                  color="white") +
        scale_fill_gradientn(colors=rev(ocean.curl(200)),
                             limits=c(-1.2,1.2),
                             oob=scales::squish,
                             name=colorbar_title,
                             guide=guide_colorbar(barwidth=unit(0.2, "cm"),
                                                  barheight=unit(2.5, "cm"),
                                                  title.position="top",
                                                  title.hjust=1,
                                                  breaks=scales::pretty_breaks(3))) +
        scale_x_continuous(expand=c(0,0),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==3 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x))) +
        scale_y_continuous(limits=function(x) c(x[1], x[2] + 25),
                           expand=c(0,0),
                           name="3,087 non-overlapping coding genes") +
        theme_default +
        theme(axis.text.y=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_text(color=ifelse(leftmost,
                                                     "black",
                                                     "white"),
                                        margin=margin(l=0, r=4, unit="pt")),
              panel.grid=element_blank(),
              panel.border=element_blank(),
              legend.justification=c(1,1),
              legend.direction="vertical",
              legend.title=element_text(margin=margin(l=-40, b=-3, unit="pt")),
              strip.text=element_blank(),
              panel.spacing.y=unit(2, "pt"),
              axis.ticks.length.y=unit(1, "pt"),
              plot.margin=margin(0, 6, 0, 0, "pt"))
    if (leftmost){
        heatmap = heatmap +
            theme(legend.position=c(0.97, 0.97))
    } else {
        heatmap = heatmap +
            theme(legend.position="none",
                  panel.grid.major.x=element_line(size=0.1,
                                                  color="gray85"))
    }

    return(heatmap)
}

plot_metagene = function(df_metagene,
                         df_quant,
                         filter_assay,
                         plot_title=expression(textstyle(frac("H3K36me3",
                                                              "H3")) ~
                                                   "ChIP-seq enrichment"),
                         legend_position="none",
                         panel_letter="x",
                         leftmost=FALSE){
    df_metagene %<>%
        filter(assay==filter_assay)
    df_quant %<>%
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
        geom_vline(data=df_quant,
                   aes(xintercept=position,
                       color=group),
                   size=0.1,
                   alpha=0.5,
                   linetype="dotted",
                   show.legend=FALSE) +
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
                           name=expression("log"[2] ~ "ratio")) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        labs(title=plot_title,
             tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.title=element_blank(),
              legend.justification=c(1,0),
              legend.position=legend_position,
              legend.background=element_blank(),
              legend.spacing.x=unit(1, "pt"),
              axis.text.y=element_text(size=5),
              axis.title.y=element_text(color=ifelse(leftmost,
                                                     "black",
                                                     "white"),
                                        margin=margin(l=0, r=-3, unit="pt")),
              plot.margin=margin(0, 6, 11/2, 0, "pt"),
              plot.title=element_text(hjust=0.5))

    return(metagene)
}


main = function(theme_path = "spn1_2020_theme.R",
                metagene_data_paths = c("verified-transcripts-nonoverlapping-TSS_ChIPseq-H3K36me3-H3norm.tsv.gz",
                                       "verified-transcripts-nonoverlapping-TSS_ChIPseq-H3K36me2-H3norm.tsv.gz",
                                       "verified-transcripts-nonoverlapping-TSS_ChIPseq-H3K4me3-H3norm.tsv.gz"),
                heatmap_data_paths = c("verified-transcripts-nonoverlapping-TSS-slopR300_ChIPseq-H3K36me3-H3norm.tsv.gz",
                                       "verified-transcripts-nonoverlapping-TSS-slopR300_ChIPseq-H3K36me2-H3norm.tsv.gz",
                                       "verified-transcripts-nonoverlapping-TSS-slopR300_ChIPseq-H3K4me3-H3norm.tsv.gz"),
                annotation_path = "Scer_transcripts_w_verifiedORFs-nonoverlapping_slopR300.bed",
                # panel_letter="c",
                df_quant_out="df_quant.tsv",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=17.4,
                fig_height=10){
    source(theme_path)
    library(cowplot)
    library(pals)

    df_heatmap = tibble()
    for (path in heatmap_data_paths){
        df_heatmap %<>% import_for_heatmap(path)
    }
    df_metagene = tibble()
    for (path in metagene_data_paths){
        df_metagene %<>% import_for_metagene(path)
    }

    annotation = read_tsv(annotation_path,
                          col_names=c("chrom", "start", "end",
                                      "name", "score", "strand")) %>%
        mutate(index = row_number()) %>%
        arrange(desc(end - start)) %>%
        mutate(sorted_index = row_number())

    df_heatmap %<>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               assay=ordered(assay,
                             levels=c("ChIPseq-H3K36me2-H3norm",
                                      "ChIPseq-H3K36me3-H3norm",
                                      "ChIPseq-H3K4me3-H3norm"),
                             labels=c("textstyle(frac(\"H3K36me2\",\"H3\"))",
                                      "textstyle(frac(\"H3K36me3\",\"H3\"))",
                                      "textstyle(frac(\"H3K4me3\",\"H3\"))"))) %>%
        left_join(annotation,
                  by="index")

    df_metagene %<>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               assay=ordered(assay,
                             levels=c("ChIPseq-H3K36me2-H3norm",
                                      "ChIPseq-H3K36me3-H3norm",
                                      "ChIPseq-H3K4me3-H3norm"),
                             labels=c("textstyle(frac(\"H3K36me2\",\"H3\"))",
                                      "textstyle(frac(\"H3K36me3\",\"H3\"))",
                                      "textstyle(frac(\"H3K4me3\",\"H3\"))")))

    df_quant_increase = df_metagene %>%
        filter(position >= 0) %>%
        group_by(group, assay) %>%
        mutate(scaled = scales::rescale(mid)) %>%
        filter(scaled >= 0.9) %>%
        arrange(position) %>%
        slice(1) %>%
        select(assay, group,
               position_90_increase=position, signal_90_increase=mid, scaled_signal_90_increase=scaled) %>%
        arrange(assay, group)
    df_quant_decrease = df_metagene %>%
        filter(position >= 0) %>%
        group_by(group, assay) %>%
        mutate(scaled = scales::rescale(mid)) %>%
        filter(scaled <= 0.1) %>%
        arrange(position) %>%
        slice(1) %>%
        select(assay, group,
               position_10_decrease=position, signal_10_decrease=mid, scaled_signal_10_decrease=scaled) %>%
        arrange(assay, group)
    df_quant = inner_join(df_quant_increase, df_quant_decrease,
                          by=c("assay", "group")) %>%
        mutate_at(vars(signal_90_increase,
                       scaled_signal_90_increase,
                       signal_10_decrease,
                       scaled_signal_10_decrease),
                  ~signif(., 4)) %>%
        write_tsv(df_quant_out)

    assays = levels(df_heatmap[["assay"]])

    plots = list(plot_metagene(df=df_metagene,
                               df_quant=select(df_quant, assay, group, position=position_90_increase),
                               filter_assay=assays[1],
                               plot_title="H3K36me2 / H3",
                               legend_position=c(0.90, 0.11),
                               panel_letter="a",
                               leftmost=TRUE),
                 plot_metagene(df=df_metagene,
                               df_quant=select(df_quant, assay, group, position=position_90_increase),
                               filter_assay=assays[2],
                               plot_title="H3K36me3 / H3",
                               panel_letter="b"),
                 plot_metagene(df=df_metagene,
                               df_quant=select(df_quant, assay, group, position=position_10_decrease),
                               filter_assay=assays[3],
                               plot_title="H3K4me3 / H3",
                               panel_letter="c"),
                 plot_heatmap(df=df_heatmap,
                              filter_assay=assays[1],
                              colorbar_title=expression("log"[2] ~
                                                            textstyle(frac("Spn1-depleted",
                                                                           "non-depleted"))),
                              leftmost=TRUE),
                 plot_heatmap(df=df_heatmap,
                              filter_assay=assays[2],
                              colorbar_title=expression("log"[2] ~
                                                            textstyle(frac("Spn1-depleted",
                                                                           "non-depleted")))),
                 plot_heatmap(df=df_heatmap,
                              filter_assay=assays[3],
                              colorbar_title=expression("log"[2] ~
                                                            textstyle(frac("Spn1-depleted",
                                                                           "non-depleted")))))

    h3_modification_datavis = plot_grid(plotlist=plots,
                                        align="v",
                                        axis="rl",
                                        nrow=2,
                                        ncol=3,
                                        rel_heights=c(0.55, 1))

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
     metagene_data_paths=snakemake@input[["metagene_data"]],
     heatmap_data_paths=snakemake@input[["heatmap_data"]],
     annotation_path=snakemake@input[["annotation"]],
     # panel_letter=snakemake@params[["panel_letter"]],
     df_quant_out=snakemake@output[["quantification"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])
