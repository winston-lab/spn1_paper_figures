import = function(df,
                  data_path,
                  control_group){
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


main = function(theme_path = "spn1_2020_theme.R",
                data_paths = c("verified-transcripts-TSS_ChIPseq-H3K36me2.tsv.gz",
                               "verified-transcripts-TSS_ChIPseq-H3K36me3.tsv.gz",
                               "verified-transcripts-TSS_ChIPseq-H3K4me3.tsv.gz"),
                df_quant_out = "df_quant.tsv",
                panel_letter="a",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=17.4 * 7 / 12,
                fig_height=12 * 3 / 12){
    source(theme_path)

    df = tibble()
    for (path in data_paths){
        df %<>% import(path, control_group="non-depleted")
    }

    df %<>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               assay=ordered(assay,
                             levels=c("ChIPseq-H3K36me2",
                                      "ChIPseq-H3K36me3",
                                      "ChIPseq-H3K4me3"),
                             labels=c("H3K36me2",
                                      "H3K36me3",
                                      "H3K4me3")))
    
    df_quant_increase = df %>%
        filter(position >= 0) %>%
        group_by(group, assay) %>%
        mutate(scaled = scales::rescale(mid)) %>%
        filter(scaled >= 0.9) %>%
        arrange(position) %>%
        slice(1) %>%
        select(assay, group,
               position_90_increase=position, signal_90_increase=mid, scaled_signal_90_increase=scaled) %>%
        arrange(assay, group)
    df_quant_decrease = df %>%
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

    h3_mods_non_h3_norm = ggplot(data=df,
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
        facet_wrap(~assay,
                   scales="free_y",
                   nrow=1) +
        scale_x_continuous(expand=c(0,0),
                           limits=c(NA, NA),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==3 ~ paste(x, "kb"),
                                                        # x>2 ~ "",
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           name=expression("log"[2] ~ textstyle(frac("IP", "input")))) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              # panel.spacing.x=unit(2, "pt"),
              # axis.ticks.length=unit(1, "pt"),
              strip.text.x=element_text(margin=margin(b=1, unit="pt")),
              legend.title=element_blank(),
              legend.justification=c(0.5,0.5),
              legend.position=c(0.20,0.38),
              legend.background=element_blank(),
              legend.key.width=unit(10, "pt"),
              legend.key.height=unit(10, "pt"),
              legend.spacing.x=unit(1, "pt"),
              axis.text.y=element_text(size=5),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              plot.margin=margin(11/2, 11/2, 0, 11/2, "pt"))

    ggsave(pdf_out,
           plot=h3_mods_non_h3_norm,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(h3_mods_non_h3_norm,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_paths=snakemake@input[["data"]],
     df_quant_out=snakemake@output[["quantification"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])

