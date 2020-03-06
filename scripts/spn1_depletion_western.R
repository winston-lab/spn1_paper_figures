
main = function(theme_path = "spn1_2020_theme.R",
                spn1_blot_path="Spn1_depletion_Spn1_9.5.tif",
                pgk1_blot_path="Spn1_depletion_Pgk1_for_Spn1_9.tif",
                # quant_data_path = "protein_level_quants_norm_Pgk1_H3.txt",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=4,
                fig_height=6,
                panel_letter="X",
                units="cm"){
    source(theme_path)
    library(grid)
    library(gridExtra)
    library(tiff)
    library(ggplotify)

    # df = read_tsv(quant_data_path) %>%
    #     rename(levels_norm_pgk1_h3 = `levels norm. to Pgk1 or H3`) %>%
    #     mutate_at(vars(strain, condition, protein),
    #               ~fct_inorder(., ordered=TRUE))
    #
    # df_rescale = df %>%
    #     filter(strain=="Spn1-AID",
    #            condition=="DMSO") %>%
    #     group_by(protein) %>%
    #     summarize(control_mean = mean(levels_norm_pgk1_h3))
    #
    # df %<>%
    #     left_join(df_rescale,
    #               by="protein") %>%
    #     group_by(protein) %>%
    #     mutate(scaled_levels = scales::rescale(levels_norm_pgk1_h3,
    #                                            from=c(0, first(control_mean))))
    #
    # df_summary = df %>%
    #     group_by(strain, condition, protein) %>%
    #     summarize(mean = mean(scaled_levels,
    #                           na.rm=TRUE),
    #               sd = sd(scaled_levels,
    #                       na.rm=TRUE))
    #
    #
    # ggplot() +
    #     geom_col(data = df_summary,
    #        aes(x=strain,
    #            fill=condition,
    #            y=mean),
    #        position=position_dodge(width=0.8),
    #        width=0.8) +
    #     geom_errorbar(data = df_summary,
    #                   aes(x=strain,
    #                       group=condition,
    #                       ymin=mean-sd,
    #                       ymax=mean+sd),
    #                   position=position_dodge(width=0.8),
    #                   width=0.2) +
    #     geom_point(data=df,
    #                 aes(x=strain,
    #                     group=condition,
    #                     y=scaled_levels),
    #                position=position_dodge(width=0.8)) +
    #     facet_wrap(~protein,
    #                scales="free_y")

    antigen_label_edge = 0.55
    condition_labels_y = 0.7
    lane_centers = c(0.15, 0.40)
    outline_height = 0.12
    antigen_labels_y = c(condition_labels_y - 0.02 - outline_height / 2,
                         condition_labels_y - 0.04 - outline_height / 2 * 3)

    antigen_labels = textGrob(x=antigen_label_edge,
                              y=antigen_labels_y,
                              label=c("Spn1", "Pgk1"),
                              hjust=0,
                              gp=gpar(fontsize=7,
                                      fontfamily="FreeSans"))
    blot_outlines = rectGrob(width = antigen_label_edge - 0.04,
                             height=outline_height,
                             x=antigen_label_edge / 2,
                             y=antigen_labels_y,
                             gp=gpar(lwd=1,
                                     fill=NA))
    spn1_image = readTIFF(spn1_blot_path)
    spn1_raster = rasterGrob(spn1_image[139:157, 337:395, 1:3],
                             width=antigen_label_edge - 0.04,
                             height=outline_height,
                             x=antigen_label_edge / 2,
                             y=antigen_labels_y[1])

    pgk1_image = readTIFF(pgk1_blot_path)
    pgk1_raster = rasterGrob(pgk1_image[217:235, 341:399, 1:3],
                             width=antigen_label_edge - 0.04,
                             height=outline_height,
                             x=antigen_label_edge / 2,
                             y=antigen_labels_y[2])

    condition_labels = textGrob(x=lane_centers,
                                y=condition_labels_y,
                                label=c("non-depleted",
                                        "Spn1-depleted"),
                                gp=gpar(fontsize=7,
                                        fontfamily="FreeSans"),
                                rot=41,
                                hjust=0.08,
                                vjust=-0.2)
    quantification_labels = textGrob(x=lane_centers,
                                     y=min(antigen_labels_y - outline_height / 2 - 0.02),
                                     label=c("0.00\n±0.00",
                                             "0.00\n±0.00"),
                                     #                                   phantom(.) %+-% "0.00" ~
                                     #                                       phantom(.)))),
                                     #         expression(textstyle(atop("0.00",
                                     #                                   phantom(.) %+-% "0.00" ~
                                     #                                       phantom(.))))),
                                     gp=gpar(fontsize=5,
                                             fontfamily="FreeSans"),
                                     hjust=0.5,
                                     vjust=1)
    lane_alignment = segmentsGrob(x0=lane_centers,
                                  x1=lane_centers,
                                  y=0,
                                  y1=1)

    western = gTree(children=gList(
        antigen_labels,
        condition_labels,
        pgk1_raster,
        spn1_raster,
        blot_outlines,
        # lane_alignment,
        quantification_labels))

    spn1_depletion_western= as.ggplot(western) +
        labs(tag=panel_letter) +
        theme(plot.margin=margin(11/2, 0, 0, 0, "pt"),
              plot.tag=element_text(family="FreeSans",
                                    size=9,
                                    face="bold"))

    ggsave(pdf_out,
           plot=spn1_depletion_western,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(spn1_depletion_western,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     spn1_blot_path=snakemake@input[["spn1_blot"]],
     pgk1_blot_path=snakemake@input[["pgk1_blot"]],
     # quant_data_path=snakemake@input[["quant_data"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])

