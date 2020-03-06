
main = function(theme_path = "spn1_2020_theme.R",
                data_path = "Spn1-chipseq_spikein-counts.tsv",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=4,
                fig_height=6,
                panel_letter="X",
                plot_title="Spn1 ChIP-seq"){
    source(theme_path)

    df = read_tsv(data_path) %>%
        mutate(abundance = (experimental_counts_IP / spikein_counts_IP) *
                   (spikein_counts_input / experimental_counts_input),
               group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")))

    baseline_mean = df %>%
        filter(group == "non-depleted") %>%
        summarize(mean=mean(abundance)) %>%
        pull(mean)

    df %<>%
        mutate(scaled_abundance=scales::rescale(abundance,
                                                from=c(0, baseline_mean)))

    df_summary = df %>%
        group_by(group) %>%
        summarize(mean_scaled_abundance = mean(scaled_abundance),
               sd_scaled_abundance = sd(scaled_abundance)) %>%
        mutate(label_direction = if_else(mean_scaled_abundance > 0.5,
                                         -1,
                                         1))

    chipseq_abundance_barplot = ggplot() +
        geom_col(data=df_summary,
                 aes(x=group,
                     y=mean_scaled_abundance,
                     fill=group),
                 alpha=0.5,
                 width=0.8) +
        geom_errorbar(data=df_summary,
                      aes(x=group,
                          ymin=mean_scaled_abundance - sd_scaled_abundance,
                          ymax=mean_scaled_abundance + sd_scaled_abundance),
                      width=0.2,
                      size=0.3,
                      alpha=0.8) +
        geom_jitter(data=df,
                    aes(x=group,
                        y=scaled_abundance),
                    width=0.2,
                    shape=16,
                    size=0.7,
                    alpha=0.8) +
        geom_text(data=df_summary,
                  aes(x=group,
                      y=mean_scaled_abundance + (sd_scaled_abundance + 0.18) * label_direction,
                      label=paste0(sprintf("%.2f", mean_scaled_abundance),
                                   "\n±",
                                  sprintf("%.2f", sd_scaled_abundance))),
                  size=5/72*25.4,
                  family="FreeSans") +
        scale_x_discrete(name=NULL,
                         expand=c(0,0.5)) +
        scale_y_continuous(limits=c(0,
                                    (max(df_summary[["mean_scaled_abundance"]] + df_summary["sd_scaled_abundance"],
                                        df[["scaled_abundance"]]) + 0.05) * 1.05),
                           expand=c(0,0),
                           breaks=c(0,1),
                           name="total signal") +
        scale_fill_viridis_d(end=0.6,
                             guide=FALSE) +
        labs(tag=panel_letter,
             title=plot_title) +
        theme_default +
        theme(panel.grid.major.x=element_blank(),
              panel.border=element_blank(),
              axis.text.x=element_text(angle=20,
                                       hjust=0.9),
              axis.ticks=element_blank())

    ggsave(pdf_out,
           plot=chipseq_abundance_barplot,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(chipseq_abundance_barplot,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_path=snakemake@input[["data"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]],
     plot_title=snakemake@params[["plot_title"]])

