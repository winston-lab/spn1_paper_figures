main = function(data_path="verified-transcripts-nonoverlapping-TSS_ChIPseq-Spn1.tsv.gz",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "C",
                fig_width=8.5,
                fig_height=9/16*8.5,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = read_tsv(data_path,
                  col_names=c("group", "sample", "annotation", "assay",
                              "index", "position", "signal")) %>%
        group_by(group, index, position) %>%
        summarize(signal = mean(signal, na.rm=TRUE)) %>%
        group_by(index) %>%
        mutate(standard_score = (signal - mean(signal, na.rm=TRUE)) / sd(signal, na.rm=TRUE)) %>%
        group_by(group, position) %>%
        summarize(low=quantile(standard_score, 0.25, na.rm=TRUE),
                  mid=median(standard_score, na.rm=TRUE),
                  high=quantile(standard_score, 0.75, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")))

    spn1_depletion_metagene = ggplot(data=df,
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
        annotate(geom="text",
                 x=2.9,
                 y=c(0.66, -0.14),
                 label=c("non-depleted",
                         "Spn1-depleted"),
                 hjust=1,
                 family="FreeSans",
                 size=7/72*25.4) +
        scale_x_continuous(expand=c(0,0),
                           labels=function(x){case_when(x==0 ~ "TSS",
                                                        x==3 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x))},
                           name=NULL) +
        scale_y_continuous(name="standard score") +
        scale_color_viridis_d(end=0.6,
                              guide=FALSE) +
        scale_fill_viridis_d(end=0.6,
                             guide=FALSE) +
        ggtitle(label="Spn1 ChIP-seq enrichment",
                subtitle="3087 non-overlapping coding genes") +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank())

    ggsave(pdf_out,
           spn1_depletion_metagene,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(spn1_depletion_metagene,
         file=grob_out)
}

main(data_path=snakemake@input[["data"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

