
main = function(theme_path = "spn1_2020_theme.R",
                data_path = "DE_rtqpcr_data.txt",
                rnaseq_path = "Spn1-IAA-v-Spn1-DMSO_rnaseq-spikenorm-transcripts-diffexp-results-genic-all.tsv",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=8.5,
                fig_height=9/16*8.5,
                panel_letter="x"){
    source(theme_path)
    library(broom)
    library(mratios)

    df = read_tsv(data_path)

    df_ratio = df %>%
        group_by(gene, condition) %>%
        mutate(replicate = row_number()) %>%
        pivot_wider(id_cols=c(gene, replicate),
                    names_from=condition,
                    values_from=signal) %>%
        group_by(gene) %>%
        do(ttestratio(.$IAA,
                      .$DMSO,
                      var.equal=FALSE) %>%
               tidy()) %>%
        rename(mean_IAA = estimate1,
               mean_DMSO = estimate2,
               mean_ratio = estimate3) %>%
        left_join(read_tsv(rnaseq_path),
                  by=c("gene"="name")) %>%
        mutate(conf_low_rnaseq = log2_foldchange + qnorm(0.025) * lfc_SE,
               conf_high_rnaseq = log2_foldchange + qnorm(0.975) * lfc_SE,
               h_just=if_else(gene %in% c("GCV3",
                                         "HSP82",
                                         "SER3",
                                         "GTR2"),
                             1, 0),
               nudge_y=case_when(gene=="GCV3" ~ -1.4,
                                 gene=="KAP120" ~ -0.8,
                                 gene=="SRB4" ~ 0.65,
                                 gene=="HSP82" ~ 0.65,
                                 gene=="UBI4" ~ -0.75,
                                 gene=="SER3" ~ 0.6,
                                 gene=="GTR2" ~ 0.65),
               nudge_x=case_when(gene=="GCV3" ~ -0.01,
                                 gene=="KAP120" ~ 0.32,
                                 gene=="SRB4" ~ 0.05,
                                 gene=="HSP82" ~ -0.05,
                                 gene=="UBI4" ~ 0.02,
                                 gene=="SER3" ~ -0.07,
                                 gene=="GTR2" ~ -0.01))

    diffexp_rtqpcr_scatter = ggplot(data=df_ratio) +
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
        geom_segment(aes(x=log2_foldchange,
                         xend=log2_foldchange,
                         y=log2(conf.low),
                         yend=log2(conf.high),
                         color=gene),
                     size=0.2,
                     alpha=0.5) +
        geom_segment(aes(x=conf_low_rnaseq,
                         xend=conf_high_rnaseq,
                         y=log2(mean_ratio),
                         yend=log2(mean_ratio),
                         color=gene),
                     size=0.2,
                     alpha=0.5) +
        geom_point(aes(x=log2_foldchange,
                       y=log2(mean_ratio),
                       color=gene),
                   size=0.7,
                   shape=16,
                   alpha=1) +
        geom_text(aes(x=log2_foldchange + nudge_x,
                       y=log2(mean_ratio) + nudge_y,
                       label=gene,
                       color=gene,
                       hjust=h_just),
                  size=5/72*25.4,
                  # label.size=NA,
                  # label.r=unit(0, "pt"),
                  # label.padding=unit(0.5, "pt"),
                  # alpha=0.8,
                  family="FreeSans",
                  fontface="italic") +
        scale_x_continuous(name=expression("RNA-seq: log"[2] ~ textstyle(frac("Spn1-depleted",
                                                                              "non-depleted"))),
                           breaks=scales::pretty_breaks(4)) +
        scale_y_continuous(breaks=scales::pretty_breaks(4),
                           name=expression(atop("RT-qPCR:",
                                                "log"[2] ~ textstyle(frac("Spn1-depleted",
                                                                          "non-depleted"))))) +
        scale_color_brewer(palette="Dark2") +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.title.y=element_text(angle=0,
                                        vjust=0.5))

    ggsave(pdf_out,
           plot=diffexp_rtqpcr_scatter,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(diffexp_rtqpcr_scatter,
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

