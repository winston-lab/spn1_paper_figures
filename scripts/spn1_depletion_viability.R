
main = function(theme_path = "spn1_2020_theme.R",
                data_path = "viability_data.txt",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=6,
                fig_height=9/16*6,
                panel_letter="x"){
    source(theme_path)

    df = read_tsv(data_path,
                  col_types="ccidiiidd",
                  na=c("N/A")) %>%
        mutate(strain = ordered(strain,
                                levels=c("WT",
                                         "noAID",
                                         "Spn1-AID"),
                                labels=c("wild type",
                                         # "Spn1-untagged",
                                         "Spn1-noAID",
                                         "Spn1-AID")),
               condition = ordered(condition,
                                   levels=c("T0",
                                            "DMSO",
                                            "IAA"),
                                   labels=c("T0",
                                            "DMSO",
                                            "IAA")),
               cfu_per_od_unit = `CFU/mL_normal` / OD600) %>%
        filter(condition != "T0")

    df_summary = df %>%
        group_by(strain, condition) %>%
        summarize(mean = mean(cfu_per_od_unit, na.rm=TRUE),
                  sd = sd(cfu_per_od_unit, na.rm=TRUE))

    viability_barplot = ggplot() +
        geom_col(data=df_summary,
                 aes(x=strain,
                     group=condition,
                     y=mean),
                 position=position_dodge(width=0.8),
                 size=0.2,
                 color="black",
                 fill="gray80",
                 width=0.8,
                 alpha=0.7) +
        geom_errorbar(data=df_summary,
                      aes(x=strain,
                          group=condition,
                          ymin=mean - sd,
                          ymax=mean + sd),
                      position=position_dodge(width=0.8),
                      width=0.2,
                      size=0.3,
                      alpha=0.8) +
        geom_point(data=df,
                   aes(x=strain,
                       color=condition,
                      y=cfu_per_od_unit),
                   position=position_jitterdodge(dodge.width=0.8),
                   shape=16,
                   size=0.5,
                   alpha=0.8) +
        scale_color_brewer(palette="Set1",
                          direction = -1,
                          name=NULL,
                          guide=guide_legend(keywidth=unit(4, "pt"),
                                             override.aes = list(size=1.5))) +
        # scale_color_manual(values=rep("black", 2),
                           # guide=FALSE) +
        scale_y_continuous(limits=c(0,
                                    max(c(df_summary[["mean"]] + df_summary[["sd"]],
                                          df[["cfu_per_od_unit"]]),
                                        na.rm=TRUE) * 1.05),
                           expand=c(0,0),
                           name=expression(textstyle(frac("CFU", "OD"[600] ~ "units"))),
                           breaks=c(0, 1e7, 2e7, 3e7),
                           labels=c(expression(scriptstyle(0)),
                                    expression(scriptstyle(1 %*% 10^7)),
                                    expression(scriptstyle(2 %*% 10^7)),
                                    expression(scriptstyle(3 %*% 10^7)))) +
        scale_x_discrete(name=NULL) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.spacing.x=unit(2, "pt"),
              legend.background=element_blank(),
              legend.box.background=element_blank(),
              legend.margin=margin(l=-6, unit="pt"),
              axis.text.x=element_text(angle=20,
                                       hjust=0.8),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5,
                                        margin=margin(r=4, unit="pt")))

    ggsave(pdf_out,
           plot=viability_barplot,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(viability_barplot,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_path=snakemake@input[["data"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])

