main = function(chip_data_path = "h3k36me3_chipqpcr_chd1D.tsv",
                theme_path = "spn1_2020_theme.R",
                fig_width=8.5,
                fig_height=4.5,
                panel_letter="b",
                pdf_out="test.pdf",
                grob_out="test.Rdata"){
    source(theme_path)

    df_qpcr = read_tsv(chip_data_path) %>%
        mutate(strain = ordered(strain,
                                levels=c("Spn1-AID",
                                         "Spn1-AID_chd1D")),
               condition=ordered(condition,
                                 levels=c("DMSO",
                                          "IAA"),
                                 labels=c("non-depleted",
                                          "Spn1-depleted")))

    df_summary = df_qpcr %>%
        group_by(strain, condition, gene) %>%
        summarize(mean = mean(signal, na.rm=TRUE),
                  sd = sd(signal, na.rm=TRUE))

    chd1_qpcr = ggplot() +
        geom_col(data=df_summary,
                 aes(x=strain,
                     group=condition,
                     y=mean),
                 position=position_dodge(width=0.5),
                 size=0.2,
                 color="black",
                 fill="gray80",
                 width=0.5,
                 alpha=0.7) +
        geom_errorbar(data=df_summary,
                      aes(x=strain,
                          group=condition,
                          ymin=mean - sd,
                          ymax=mean + sd),
                      position=position_dodge(width=0.5),
                      width=0.2,
                      size=0.3,
                      alpha=0.8) +
        geom_point(data=df_qpcr,
                   aes(x=strain,
                       color=condition,
                      y=signal),
                   position=position_jitterdodge(dodge.width=0.5),
                   shape=16,
                   size=0.8,
                   alpha=0.8) +
        facet_wrap(~gene) +
        scale_color_viridis_d(end=0.6,
                              name=NULL,
                              guide=guide_legend(keywidth=unit(4, "pt"),
                                                 override.aes = list(size=1.5))) +
        scale_y_continuous(limits=function(x) c(0, x[2] * 1.05),
                           expand=c(0,0),
                           breaks=scales::pretty_breaks(3),
                           name=quote(atop("H3K36me3:", textstyle(frac("IP", "input"))))) +
        scale_x_discrete(name=NULL,
                         labels=c(quote(textstyle(atop("Spn1-AID", phantom(".")))),
                                  quote(textstyle(atop("Spn1-AID", italic("chd1Δ"))))),
                         expand=c(0.5,0)) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              panel.spacing.x=unit(2, "pt"),
              legend.spacing.x=unit(1, "pt"),
              legend.background=element_blank(),
              legend.box.background=element_blank(),
              legend.text=element_text(size=7),
              legend.position=c(0.66, 0.6),
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(size=10),
              strip.text.x=element_text(face="italic",
                                        margin=margin(b=2, unit="pt")),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5,
                                        margin=margin(r=0, unit="pt")))

    ggsave(pdf_out,
           chd1_qpcr,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(chd1_qpcr,
         file=grob_out)
}

main(chip_data_path=snakemake@input[["chip_data"]],
     theme_path=snakemake@input[["theme"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

