import = function(df,
                  data_path,
                  control_group){
    df_temp = read_tsv(data_path,
             col_names=c("group",
                         "sample",
                         "annotation",
                         "assay",
                         "index",
                         "position",
                         "signal"))

    df_mean_sd = df_temp %>%
        filter(group == control_group) %>%
        group_by(index) %>%
        summarize(control_mean = mean(signal, na.rm=TRUE),
                  control_sd  = sd(signal, na.rm=TRUE))

    df_temp %>%
        group_by(group, assay, index, position) %>%
        summarize(signal=mean(signal, na.rm=TRUE)) %>%
        left_join(df_mean_sd,
                  by="index") %>%
        mutate(standard_score = (signal - control_mean) / control_sd) %>%
        group_by(group, assay, position) %>%
        summarize(low=quantile(standard_score, 0.25, na.rm=TRUE),
                  mid=median(standard_score, na.rm=TRUE),
                  high=quantile(standard_score, 0.75, na.rm=TRUE)) %>%
        bind_rows(df, .) %>%
        return()
}


main = function(theme_path = "spn1_2020_theme.R",
                data_paths = c("verified-transcripts-TSS_ChIPseq-H3.tsv.gz",
                               "verified-transcripts-nonoverlapping-TSS_ChIPseq-H3K4me3-H3norm.tsv.gz",
                               "verified-transcripts-nonoverlapping-TSS_ChIPseq-H3K36me2-H3norm.tsv.gz",
                               "verified-transcripts-nonoverlapping-TSS_ChIPseq-H3K36me3-H3norm.tsv.gz"),
                panel_letter="A",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=8.5,
                fig_height=8){
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
                             levels=c("ChIPseq-H3",
                                      "ChIPseq-H3K4me3-H3norm",
                                      "ChIPseq-H3K36me2-H3norm",
                                      "ChIPseq-H3K36me3-H3norm"),
                             labels=c("\"H3\"",
                                      "textstyle(frac(\"H3K4me3\",\"H3\"))",
                                      "textstyle(frac(\"H3K36me2\",\"H3\"))",
                                      "textstyle(frac(\"H3K36me3\",\"H3\"))")))


    histone_metagenes = ggplot(data=df,
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
                           name="standard score") +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        facet_grid(assay~.,
                   scales="free_y",
                   labeller=label_parsed) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.title=element_blank(),
              legend.justification=c(1,0.5),
              legend.position=c(0.95,0.86),
              legend.background=element_blank(),
              legend.spacing.x=unit(1, "pt"),
              strip.text.y=element_text(angle=0, hjust=0),
              panel.spacing.y=unit(2, "pt"),
              axis.text.y=element_text(size=5))

    ggsave(pdf_out,
           plot=histone_metagenes,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(histone_metagenes,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_paths=snakemake@input[["data"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])

