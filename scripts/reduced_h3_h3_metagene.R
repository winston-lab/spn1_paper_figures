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
        group_by(group, annotation, assay, position) %>%
        summarize(low=quantile(signal, 0.25, na.rm=TRUE),
                  mid=median(signal, na.rm=TRUE),
                  high=quantile(signal, 0.75, na.rm=TRUE)) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(theme_path = "custom_theme.R",
                data_paths = c("plmin-reduced-H3_ChIPseq-H3.tsv.gz"),
                panel_letter="b",
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=17.4/2,
                fig_height=(17.4/2) * 9 / 16){
    source(theme_path)

    df = tibble()
    for (path in data_paths){
        df %<>% import(path)
    }

    df %<>%
        ungroup() %>%
        filter(annotation=="genes with reduced H3") %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")))

    reduced_h3_h3_metagene = ggplot(data=df,
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
                           breaks=scales::pretty_breaks(4),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        x==2 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           name=expression("log"[2] ~ textstyle(frac("IP", "input")))) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        labs(tag=panel_letter,
             title="H3 ChIP-seq",
             subtitle="77 genes with reduced H3 enrichment") +
        theme_default +
        theme(panel.grid=element_blank(),
              # panel.spacing.x=unit(20, "pt"),
              legend.title=element_blank(),
              legend.justification=c(0.5,0.5),
              legend.position=c(0.6,0.55),
              # legend.position="top",
              legend.background=element_blank(),
              legend.spacing.x=unit(1, "pt"),
              axis.text.y=element_text(size=5),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              # plot.margin=margin(0, 11/2, 0, 11/2, "pt")
              )

    ggsave(pdf_out,
           plot=reduced_h3_h3_metagene,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(reduced_h3_h3_metagene,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_paths=snakemake@input[["data"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]])

