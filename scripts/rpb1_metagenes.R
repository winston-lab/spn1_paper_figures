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
                rpb1_path = "verified-transcripts-nonoverlapping-TSS_ChIPseq-Rpb1.tsv.gz",
                modification_paths= c("verified-transcripts-nonoverlapping-TSS_ChIPseq-Ser5P-Rpb1norm.tsv.gz",
                                      "verified-transcripts-nonoverlapping-TSS_ChIPseq-Ser2P-Rpb1norm.tsv.gz"),
                panel_letter="x",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=8.5,
                fig_height=8){
    source(theme_path)
    library(cowplot)
    library(ggplotify)

    df_rpb1 = read_tsv(rpb1_path,
                       col_names=c("group",
                                   "sample",
                                   "annotation",
                                   "assay",
                                   "index",
                                   "position",
                                   "signal"))

    df_rpb1_mean_sd = df_rpb1 %>%
        filter(group == "non-depleted") %>%
        group_by(index) %>%
        summarize(control_mean = mean(signal, na.rm=TRUE),
                  control_sd  = sd(signal, na.rm=TRUE))

    df_rpb1 %<>%
        group_by(group, assay, index, position) %>%
        summarize(signal=mean(signal, na.rm=TRUE)) %>%
        left_join(df_rpb1_mean_sd,
                  by="index") %>%
        mutate(standard_score = (signal - control_mean) / control_sd) %>%
        group_by(group, assay, position) %>%
        summarize(low=quantile(standard_score, 0.25, na.rm=TRUE),
                  mid=median(standard_score, na.rm=TRUE),
                  high=quantile(standard_score, 0.75, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               assay=ordered(assay,
                             levels="ChIPseq-Rpb1-spikenorm-batchB",
                             labels="\"Rpb1\""))

    rpb1_metagene = ggplot(data=df_rpb1,
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
        facet_grid(assay~.,
                   scales="free_y",
                   labeller=label_parsed) +
        scale_x_continuous(expand=c(0,0),
                           limits=c(NA, NA),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==3 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           name="(IP/input) z-score") +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        theme_default +
        theme(panel.grid=element_blank(),
              strip.text.y=element_text(angle=0, hjust=0),
              legend.title=element_blank(),
              legend.justification=c(0.5,0.5),
              legend.position=c(0.55,0.3),
              legend.background=element_blank(),
              legend.spacing.x=unit(1, "pt"),
              axis.text.y=element_text(size=5),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin=margin(b=-0.5, unit="pt"))

    df_mods = tibble()
    for (path in modification_paths){
        df_mods %<>% import(path)
    }

    df_mods %<>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               assay=ordered(assay,
                             levels=c("ChIPseq-Ser5P-Rpb1norm",
                                      "ChIPseq-Ser2P-Rpb1norm"),
                             # labels=c("textstyle(frac(\"Rpb1-Ser5P\",\"Rpb1\"))",
                                      # "textstyle(frac(\"Rpb1-Ser2P\",\"Rpb1\"))")))
                             labels=c("frac(\"Ser5-P\",\"Rpb1\")",
                                      "frac(\"Ser2-P\",\"Rpb1\")")))

    modification_metagenes = ggplot(data=df_mods,
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
                           name=expression("log"[2] ~ "ratio")) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        facet_grid(assay~.,
                   scales="free_y",
                   labeller=label_parsed) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.title=element_blank(),
              legend.position="none",
              strip.text.y=element_text(angle=0, hjust=0),
              panel.spacing.y=unit(2, "pt"),
              axis.text.y=element_text(size=5),
              axis.title.y=element_text(margin=margin(l=0, r=-3, unit="pt")),
              plot.margin=margin(t=-0.5, unit="pt"))

    rpb1_metagenes = plot_grid(rpb1_metagene,
                               modification_metagenes,
                               ncol=1,
                               rel_heights=c(1, 2),
                               align="v",
                               axis="lr") %>%
        as.ggplot()

    rpb1_metagenes = rpb1_metagenes +
        labs(tag=panel_letter) +
        theme(plot.tag=element_text(family="FreeSans",
                                    size=9,
                                    face="bold"),
              plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"))

    ggsave(pdf_out,
           plot=rpb1_metagenes,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(rpb1_metagenes,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     rpb1_path=snakemake@input[["rpb1"]],
     modification_paths=snakemake@input[["modifications"]] ,
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])
