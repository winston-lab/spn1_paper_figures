main = function(data_path="Spn1-IAA-v-Spn1-DMSO_rnaseq-spikenorm-verified-coding-genes-diffexp-results-all.tsv",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "a",
                fig_width=8.5,
                fig_height=9/16*8.5,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)
    library(cowplot)

    df = read_tsv(data_path) %>%
        filter(control_expr > 0)

    lfc_density = ggplot(data=df,
           aes(x=log2_foldchange)) +
        geom_vline(xintercept=0,
                   color="gray70",
                   size=0.2) +
        geom_density(aes(y=..scaled..),
                     size=0.2,
                     fill="#114477") +
        scale_x_continuous(limits=range(df[["log2_foldchange"]]),
                           breaks=seq(-4,4,2)) +
        scale_y_continuous(limits=c(-0.02,1.1),
                           expand=c(0,0),
                           name="density") +
        coord_flip() +
        theme_default +
        theme(axis.title.y=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              panel.grid=element_blank(),
              plot.margin=margin(l=-11, r=11, unit="pt"))

    maplot = ggplot(data=df,
           aes(x=control_expr,
               y=log2_foldchange)) +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=153,
                    shape=16,
                    size=0.3,
                    alpha=0.8) +
        # geom_rug(sides="l",
        #          size=0.1,
        #          alpha=0.1) +
        scale_x_log10(breaks=c(1e1, 1e3, 1e5),
                      labels=c(expression(10^1),
                               expression(10^3),
                               expression(10^5)),
                      name="non-depleted transcript levels") +
        scale_y_continuous(limits=range(df[["log2_foldchange"]]),
                           breaks=seq(-4,4,2),
                           name=expression(atop("RNA-seq:",
                                                "log"[2] ~
                                                    textstyle(frac("Spn1-depleted",
                                                                   "non-depleted"))))) +
        scale_color_viridis_c(option="inferno") +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.title.y=element_text(angle=0,
                                        hjust=1,
                                        vjust=0.5),
              axis.text.x=element_text(margin=margin(t=-0.5, unit="pt")))

    rnaseq_maplot = plot_grid(maplot,
                              lfc_density,
                              align="h",
                              axis="tb",
                              nrow=1,
                              rel_widths=c(0.9,0.1))

    ggsave(pdf_out,
           rnaseq_maplot,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(rnaseq_maplot,
         file=grob_out)
}

main(data_path=snakemake@input[["data"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

