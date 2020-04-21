
main = function(theme_path = "spn1_2020_theme.R",
                data_path = "promoter_swap_rtqpcr_data.txt",
                rnaseq_path="Spn1-IAA-v-Spn1-DMSO_rnaseq-spikenorm-transcripts-diffexp-results-genic-all.tsv",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=8.5,
                fig_height=9/16*8.5,
                panel_letter="x"){
    source(theme_path)
    library(broom)
    library(mratios)

    df = read_tsv(data_path) %>%
        mutate(strain=ordered(strain,
                              levels=c("parent",
                                       "noDEL",
                                       "pUBI4",
                                       "pGCV3-1",
                                       "pGCV3-2"),
                              labels=c("\"wild type\"[\"unmarked\"]",
                                       "\"wild type\"[\"marked\"]",
                                       "\"p\" * italic(\"UBI4\")",
                                       "\"p\" * italic(\"GCV3\")[\"[-232, -1]\"]",
                                       "\"p\" * italic(\"GCV3\")[\"[-232, +90]\"]")),
               condition=fct_inorder(condition,
                                     ordered=TRUE))

    df_ratio = df %>%
        pivot_wider(id_cols=c(strain, replicate),
                    names_from=condition,
                    values_from=signal) %>%
        group_by(strain) %>%
        do(ttestratio(.$IAA,
                      .$DMSO,
                      var.equal=TRUE) %>%
               tidy()) %>%
        rename(mean_IAA = estimate1,
               mean_DMSO = estimate2,
               mean_ratio = estimate3) %>%
        mutate(gene = case_when(str_detect(strain, "YLR454W") | str_detect(strain, "wild type")  ~ "FMP27",
                                str_detect(strain, "UBI4") ~ "UBI4",
                                str_detect(strain, "GCV3") ~ "GCV3")) %>%
        left_join(read_tsv(rnaseq_path),
                  by=c("gene"="name")) %>%
        mutate(conf_low_rnaseq = log2_foldchange + qnorm(0.025) * lfc_SE,
               conf_high_rnaseq = log2_foldchange + qnorm(0.975) * lfc_SE,
               gene = ordered(gene,
                              levels=c("GCV3",
                                       "FMP27",
                                       "UBI4"))) %>%
        group_by(gene) %>%
        arrange(mean_ratio) %>%
        mutate(label_hjust = if_else(n() == 1,
                                     0,
                                     1),
               label_hjust = label_hjust * row_number() %% 2,
               jitter_direction = scales::rescale(label_hjust, to=c(-1, 1)),
               x_jitter = jitter_direction * 0.02)

    promoter_swap_scatter = ggplot(data=df_ratio) +
        geom_hline(yintercept=0,
                   color="gray70",
                   size=0.2) +
        geom_vline(xintercept=0,
                   color="gray70",
                   size=0.2) +
        geom_abline(slope=1,
                    intercept=0,
                    color="gray70",
                    size=0.2) +
        geom_segment(aes(x=log2_foldchange + x_jitter,
                         xend=log2_foldchange + x_jitter,
                         y=log2(conf.low),
                         yend=log2(conf.high),
                         color=gene),
                     size=0.2,
                     alpha=0.5) +
        geom_label(aes(x=log2_foldchange + x_jitter * 3 - if_else(x_jitter==0, 0.05, 0),
                       y=log2(mean_ratio) - x_jitter * 3,
                       label=strain,
                       hjust=(1 - label_hjust),
                       vjust=scales::rescale(label_hjust, to=c(0.15, 1.10))),
                   size=5/72*25.4,
                   label.size=NA,
                   label.r=unit(0, "pt"),
                   label.padding=unit(0, "pt"),
                   alpha=0.8,
                   parse=TRUE,
                   family="FreeSans") +
        geom_segment(aes(x=conf_low_rnaseq + x_jitter,
                         xend=conf_high_rnaseq + x_jitter,
                         y=log2(mean_ratio),
                         yend=log2(mean_ratio),
                         color=gene),
                     size=0.2,
                     alpha=0.5) +
        geom_point(aes(x=log2_foldchange + x_jitter,
                       y=log2(mean_ratio),
                       color=gene),
                   size=0.5) +
        scale_x_continuous(limits=c(-4.15, NA),
                           name=expression("RNA-seq: log"[2] ~ textstyle(frac("Spn1-depleted",
                                                                              "non-depleted")) ~
                                               ", native gene"),
                           breaks=scales::pretty_breaks(4)) +
        scale_y_continuous(breaks=scales::pretty_breaks(4),
                           name=expression(atop("RT-qPCR:",
                                                displaystyle(atop("log"[2] ~ textstyle(frac("Spn1-depleted",
                                                                          "non-depleted") ~ ","),
                                                                  "gene fusion"))))) +
        scale_color_brewer(palette="Set1") +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.title.y=element_text(angle=0,
                                        vjust=0.5))

    ggsave(pdf_out,
           plot=promoter_swap_scatter,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(promoter_swap_scatter,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_path=snakemake@input[["data"]],
     rnaseq_path=snakemake@input[["rnaseq"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])

