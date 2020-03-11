
main = function(theme_path = "spn1_2020_theme.R",
                data_path = "promoter_swap_diagram.tsv",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=8.5,
                fig_height=9/16*8.5,
                panel_letter="a"){
    source(theme_path)

    x_limits = c(-1050, 500)
    df = read_tsv(data_path) %>%
        mutate(strain=ordered(strain,
                              levels=c("parent",
                                       "noDEL",
                                       "pUBI4",
                                       "pGCV3-1",
                                       "pGCV3-2"),
                              labels=c("\"wild type\"[\"unmarked\"]",
                                       "\"wild type\"[\"marked\"]",
                                       "\"pUBI4\"",
                                       "\"pGCV3\"[\"[-232, -1]\"]",
                                       "\"pGCV3\"[\"[-232, +90]\"]")),
               promoter=factor(promoter,
                               levels=c("UBI4",
                                        "FMP27",
                                        "GCV3-1",
                                        "GCV3-2"),
                               labels=c("\"pUBI4\"",
                                        "\"pYLR454W\"",
                                        "\"pGCV3\"[scriptscriptstyle(\"[-232, -1]\")]",
                                        "\"pGCV3\"[scriptscriptstyle(\"[-232, +90]\")]")),
               orf = ordered(orf,
                             levels="FMP27",
                             labels="YLR454W"),
               orf_label_x = if_else(orf_end > x_limits[2],
                                     (orf_start + x_limits[2]) / 2,
                                     (orf_start + orf_end) / 2),
               promoter_gene = ordered(promoter_gene,
                                       levels=c("GCV3",
                                                "FMP27",
                                                "UBI4")))

    df_natmx = df %>%
        filter(! is.na(natmx_start)) %>%
        mutate(notch = natmx_start + 0.85 * (natmx_end - natmx_start)) %>%
        expand(nesting(strain,
                       natmx_start,
                       notch,
                       natmx_end),
               natmx_y=c(-1,1)) %>%
        pivot_longer(c(natmx_start, notch, natmx_end),
                     names_to="type",
                     values_to="natmx_x") %>%
        mutate(natmx_y = ifelse(type=="natmx_end",
                                0,
                                natmx_y)) %>%
        distinct() %>%
        group_by(strain) %>%
        arrange(type,
                natmx_y,
                .by_group=TRUE) %>%
        mutate(order = c(1,3,4,2,5)) %>%
        arrange(order,
                .by_group=TRUE) %>%
        mutate(natmx_x = pmax(x_limits[1],
                              pmin(x_limits[2], natmx_x, na.rm=TRUE),
                              na.rm=TRUE))

    promoter_swap_diagram = ggplot(data=df) +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_polygon(data=df_natmx,
                     aes(x=natmx_x,
                         y=natmx_y,
                         group=strain),
                     fill="gray80") +
        geom_text(aes(x=(x_limits[1] + 0.85 * natmx_end + 0.15 * natmx_start) / 2,
                      y=0),
                  label="natMX6",
                  size=6/72*25.4,
                  family="FreeSans",
                  fontface="italic") +
        geom_rect(aes(xmin=orf_start,
                      xmax=pmin(orf_end, x_limits[2]),
                      ymin=-1,
                      ymax=1),
                  fill="gray80") +
        geom_text(aes(x=orf_label_x,
                      y=0,
                      label=orf),
                  size=6/72*25.4,
                  family="FreeSans",
                  fontface="italic") +
        geom_segment(aes(x=promoter_start,
                         xend=promoter_end,
                         color=promoter_gene),
                     y=0,
                     yend=0,
                     size=1) +
        geom_text(aes(x=(promoter_start + promoter_end) / 2,
                      y=0.1,
                      label=promoter),
                  size=5/72*25.4,
                  nudge_x=-10,
                  vjust=0,
                  parse=TRUE,
                  family="FreeSans") +
        geom_vline(xintercept=x_limits,
                   size=0.2,
                   color="gray50") +
        facet_grid(strain~.,
                   switch="y",
                   labeller=label_parsed) +
        scale_x_continuous(expand=c(0,0),
                           limits=x_limits) +
        scale_y_continuous(limits=c(-1.5,1.5)) +
        scale_color_brewer(palette="Set1") +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              panel.border=element_blank(),
              panel.spacing.y=unit(1, "pt"),
              axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              strip.text.y=element_text(angle=-180,
                                        hjust=1),
              legend.position="none")

    ggsave(pdf_out,
           plot=promoter_swap_diagram,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(promoter_swap_diagram,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_path=snakemake@input[["data"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])

