
main = function(theme_path = "spn1_2020_theme.R",
                data_path = "IR_rtqpcr_data.txt",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=8.5,
                fig_height=9/16*8.5,
                panel_letter="x"){
    source(theme_path)

    df = read_tsv(data_path) %>%
        mutate(genotype = ordered(genotype,
                                  levels=c("WT", "Spn1-AID"),
                                  labels=c("wild type", "Spn1-AID")),
               group = glue::glue("{genotype} + {condition}"),
               group = fct_inorder(group,
                                   ordered=TRUE),
               gene = fct_inorder(gene,
                                  ordered=TRUE))

    df_summary = df %>%
        group_by(gene, group) %>%
        summarize(mean = mean(signal),
                  sd = sd(signal))

    splicing_rtqpcr_barplot = ggplot() +
        geom_col(data=df_summary,
                 aes(x=gene,
                     group=group,
                     y=mean),
                 position=position_dodge(width=0.8),
                 size=0.2,
                 color="black",
                 fill="gray80",
                 width=0.8,
                 alpha=0.7) +
        geom_errorbar(data=df_summary,
                      aes(x=gene,
                          group=group,
                          ymin=mean - sd,
                          ymax=mean + sd),
                      position=position_dodge(width=0.8),
                      width=0.2,
                      size=0.3,
                      alpha=0.8) +
        geom_point(data=df,
                   aes(x=gene,
                       color=group,
                       y=signal),
                   position=position_jitterdodge(dodge.width=0.8),
                   shape=16,
                   size=0.7,
                   alpha=0.8) +
        scale_color_manual(values=c("#377eb8",
                                    "#BB5566"),
                           name=NULL,
                           guide=guide_legend(keywidth=unit(4, "pt"),
                                              override.aes = list(size=1.5))) +
        scale_y_continuous(limits=function(x) c(0, x[2] * 1.05),
                           expand=c(0,0),
                           breaks=scales::pretty_breaks(3),
                           name=expression(atop("RT-qPCR:",
                                                textstyle(frac("unspliced",
                                                               "spliced"))))) +
        scale_x_discrete(name=NULL) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.text.x=element_text(face="italic",
                                       size=6),
              axis.ticks.length.x=unit(1, "pt"),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5,
                                        margin=margin(r=0, unit="pt")),
              legend.spacing.x=unit(2, "pt"),
              legend.justification=c(0.5, 1),
              legend.position=c(4/6, 0.95))

    ggsave(pdf_out,
           plot=splicing_rtqpcr_barplot,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(splicing_rtqpcr_barplot,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_path=snakemake@input[["data"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])
