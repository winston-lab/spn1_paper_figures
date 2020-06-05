
main = function(theme_path = "spn1_2020_theme.R",
                spn1_blot_path="All replicates - raw western blots/Rep1/Spn1_Rep1.tif",
                rpb1_blot_path="All replicates - raw western blots/Rep3/Rpb1-Pgk1_Rep3.tif",
                spt6_blot_path="All replicates - raw western blots/Rep3/Spt6_Rep3.tif",
                set2_blot_path="All replicates - raw western blots/Rep2/Set2_Rep2.tif",
                pgk1_blot_path="All replicates - raw western blots/Rep3/Pgk1_for_Spn1_Rep3.tif",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=6,
                fig_height=6,
                panel_letter="X"){
    source(theme_path)
    library(grid)
    library(gridExtra)
    library(tiff)
    library(ggplotify)

    y_margin_size = 0.0225
    blot_max_y = 0.88
    blot_min_y = 0.04
    blot_height = ((blot_max_y - blot_min_y) - 4 * y_margin_size) / 6.2
    blot_pixel_height = 30
    blot_pixel_width = 169
    antigen_labels_x = 0.75
    antigen_labels_y = blot_max_y - c((2.45 / 2) * blot_height,
                                      (2.45 + 0.5) * blot_height + y_margin_size,
                                      (2.45 + 1.5) * blot_height + 2 * y_margin_size,
                                      (2.45 + 2.5) * blot_height + 3 * y_margin_size,
                                      (2.45 + 3.5) * blot_height + 4 * y_margin_size)
    blot_width = antigen_labels_x - 0.02

    antigen_labels = textGrob(x=antigen_labels_x,
                              y=antigen_labels_y,
                              label=c("", "Rpb1", "Spt6", "Set2", "Pgk1"),
                              hjust=0,
                              gp=gpar(fontsize=7,
                                      fontfamily="FreeSans"))
    spn1_labels = textGrob(x=antigen_labels_x,
                           y=c(antigen_labels_y[1] + 0.09,
                               antigen_labels_y[1] - 0.09),
                           label=c("Spn1-AID", "Spn1"),
                           hjust=0,
                           gp=gpar(fontsize=7,
                                   fontfamily="FreeSans"))
    lane_centers = scales::rescale(seq(0.5/6, 5.5/6, length.out=6),
                                   from=c(0,1),
                                   to=c(0.01, 0.01 + blot_width))

    lane_alignment = segmentsGrob(x0=lane_centers,
                                  x1=lane_centers,
                                  y=0,
                                  y1=1)

    condition_labels = textGrob(x=lane_centers,
                                y=blot_max_y + (y_margin_size / 2),
                                label=rep(c("DMSO", "IAA"), 3),
                                vjust=0,
                                gp=gpar(fontsize=5,
                                        fontfamily="FreeSans"))

    strain_centers = lane_centers %>%
        matrix(ncol=2, byrow=TRUE) %>%
        rowMeans()

    strain_lines = segmentsGrob(x0=strain_centers - blot_width / 6 + 0.008,
                                x1=strain_centers + blot_width / 6 - 0.008,
                                y0=0.935,
                                y1=0.935)
    strain_labels = textGrob(x=strain_centers,
                             y=1,
                             label=c("wild type",
                                     "Spn1-noAID",
                                     "Spn1-AID"),
                             vjust=1,
                             gp=gpar(fontsize=7,
                                     fontfamily="FreeSans"))

    blot_outlines = rectGrob(width = blot_width,
                             height = c(blot_height * 2.45, rep(blot_height, 4)) ,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y,
                             gp=gpar(lwd=1,
                                     fill=NA))

    spn1_image = readTIFF(spn1_blot_path)
    spn1_raster = rasterGrob(spn1_image[95:(95 + 2.45 * blot_pixel_height),
                                        53:(53 + blot_pixel_width),
                                        1:3],
                             width=blot_width,
                             height=blot_height * 2.45,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[1])

    rpb1_image = readTIFF(rpb1_blot_path)
    rpb1_raster = rasterGrob(rpb1_image[50:(50 + blot_pixel_height),
                                        55:(55 + blot_pixel_width),
                                        1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[2])

    spt6_image = readTIFF(spt6_blot_path)
    spt6_raster = rasterGrob(spt6_image[50:(50 + blot_pixel_height),
                                        24:(24 + blot_pixel_width - 5),
                                        1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[3])

    set2_image = readTIFF(set2_blot_path)
    set2_raster = rasterGrob(set2_image[104:(104+blot_pixel_height),
                                        17:(17+blot_pixel_width),
                                        1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[4])

    pgk1_image = readTIFF(pgk1_blot_path)
    pgk1_raster = rasterGrob(pgk1_image[175:(175 + blot_pixel_height),
                                        47:(47 + blot_pixel_width),
                                        1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[5])

    western = gTree(children=gList(
        antigen_labels,
        spn1_labels,
        strain_lines,
        strain_labels,
        condition_labels,
        spn1_raster,
        rpb1_raster,
        spt6_raster,
        set2_raster,
        pgk1_raster,
        blot_outlines#,
        # lane_alignment
        ))

    elongation_spn1_depletion_westerns = as.ggplot(western) +
        labs(tag=panel_letter) +
        theme(plot.margin=margin(11/2, 0, 0, 11/2, "pt"),
              plot.tag=element_text(family="FreeSans",
                                    size=9,
                                    face="bold"))

    ggsave(pdf_out,
           plot=elongation_spn1_depletion_westerns,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(elongation_spn1_depletion_westerns,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     spn1_blot_path=snakemake@input[["spn1_blot"]],
     rpb1_blot_path=snakemake@input[["rpb1_blot"]],
     spt6_blot_path=snakemake@input[["spt6_blot"]],
     set2_blot_path=snakemake@input[["set2_blot"]],
     pgk1_blot_path=snakemake@input[["pgk1_blot"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])
