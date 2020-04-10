
main = function(theme_path = "spn1_2020_theme.R",
                spn1_blot_path="Spn1_depletion_Spn1_9.5.tif",
                pgk1_blot_path="Spn1_depletion_Pgk1_for_Spn1_9_rot_2.2.tif",
                quant_data_path = "spn1_depletion_spn1_western_quantification.tsv",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=11.4/6,
                fig_height=9/16*11.4*8/12,
                panel_letter="X",
                units="cm"){
    source(theme_path)
    library(grid)
    library(gridExtra)
    library(tiff)
    library(ggplotify)

    df = read_tsv(quant_data_path) %>%
        mutate(area = bounding_box_width * bounding_box_height) %>%
        select(strain, condition, measurement_type, antigen, replicate,
               raw_integrated_density, area) %>%
        pivot_wider(names_from=measurement_type,
                    values_from=c(raw_integrated_density,
                                  area)) %>%
        mutate(norm_integrated_density_background = raw_integrated_density_background / area_background * area_band) %>%
        select(-c(raw_integrated_density_background,
                  area_band,
                  area_background)) %>%
        mutate(background_subtracted_signal = raw_integrated_density_band - norm_integrated_density_background,
               strain = fct_inorder(strain, ordered=TRUE),
               condition = fct_inorder(condition, ordered=TRUE))

    df_rescale = df %>%
        filter(strain=="Spn1-AID",
               condition=="DMSO") %>%
        select(antigen,
               replicate,
               rescale_value = background_subtracted_signal)

    df %<>%
        left_join(df_rescale,
                  by=c("antigen", "replicate")) %>%
        group_by(antigen, replicate) %>%
        mutate(rescaled = scales::rescale(pmax(background_subtracted_signal, 0),
                                          from=c(0, first(rescale_value)))) %>%
        select(strain, condition, antigen, replicate, rescaled) %>%
        pivot_wider(names_from=antigen,
                    values_from=rescaled) %>%
        mutate(test=Spn1/Pgk1) %>%
        arrange(strain, condition)

    # replicate 3 is excluded because the signal to noise ratio
    # for that replicate is much lower in the Spn1 blot
    df_summary = df %>%
        filter(replicate != 3) %>%
        group_by(strain, condition) %>%
        summarize(mean=mean(test),
                  sd=sd(test))

    antigen_label_edge = 0.55
    condition_labels_y = 0.62
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
    spn1_raster = rasterGrob(spn1_image[139:(139+18), 337:(337+58), 1:3],
                             width=antigen_label_edge - 0.04,
                             height=outline_height,
                             x=antigen_label_edge / 2,
                             y=antigen_labels_y[1])

    pgk1_image = readTIFF(pgk1_blot_path)
    pgk1_raster = rasterGrob(pgk1_image[220:(220+18), 341:(341+58), 1:3],
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
                                rot=49,
                                hjust=0.02,
                                vjust=-0.2)
    quantification_labels = textGrob(x=lane_centers,
                                     y=min(antigen_labels_y - outline_height / 2 - 0.02),
                                     label=c("1",
                                             "0.08\n±0.02"),
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
        theme(plot.margin=margin(11/2, 0, 0, 11/2, "pt"),
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
     quant_data_path=snakemake@input[["quant_data"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])
