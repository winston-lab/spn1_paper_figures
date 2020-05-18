main = function(data_path="verified-transcripts-nonoverlapping-TSS_ChIPseq-Set2-Rpb1norm-batchAB.tsv.gz",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "b",
                fig_width=17.4 * 5 / 12,
                fig_height=13 * 4 / 12,
                numerator_factor="Spt6",
                denominator_factor="Rpb1",
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = read_tsv(data_path,
                  col_names=c("group", "sample", "annotation", "assay",
                              "index", "position", "signal")) %>%
        group_by(group, index, position) %>%
        summarize(signal = mean(signal, na.rm=TRUE)) %>%
        group_by(group, position) %>%
        summarize(low=quantile(signal, 0.25, na.rm=TRUE),
                  mid=median(signal, na.rm=TRUE),
                  high=quantile(signal, 0.75, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")))

    metagene = ggplot(data=df,
           aes(x=position,
               y=mid,
               ymin=low,
               ymax=high,
               fill=group,
               color=group)) +
        geom_vline(xintercept=0,
                   color="gray70",
                   size=0.2) +
        geom_ribbon(linetype="blank",
                    alpha=0.15) +
        geom_line(alpha=0.9,
                  size=0.5) +
        scale_x_continuous(expand=c(0,0),
                           labels=function(x){case_when(x==0 ~ "TSS",
                                                        x==3 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x))},
                           name=NULL) +
        scale_y_continuous(name=bquote("log"[2] ~ frac(.(numerator_factor),
                                                       .(denominator_factor)))) +
        scale_color_viridis_d(end=0.6) +
        scale_fill_viridis_d(end=0.6) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.title=element_blank(),
              legend.justification=c(1,0),
              legend.position=c(0.99, 0.01),
              legend.key.width=unit(14, "pt"),
              legend.background=element_blank(),
              legend.spacing.x=unit(1, "pt"),
              axis.text.y=element_text(size=5),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5,
                                        margin=margin(l=-4, unit="pt")),
              plot.margin=margin(11/2, 11/2+1, 11/2+4, 11/2, "pt"))

    ggsave(pdf_out,
           metagene,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(metagene,
         file=grob_out)
}

main(data_path=snakemake@input[["data"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     numerator_factor=snakemake@params[["numerator_factor"]],
     denominator_factor=snakemake@params[["denominator_factor"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])
