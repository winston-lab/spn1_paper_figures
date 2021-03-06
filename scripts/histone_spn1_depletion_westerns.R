
main = function(theme_path = "spn1_2020_theme.R",
                spn1_blot_path="All replicates - raw western blots/Rep1/Spn1_Rep1.tif",
                pgk1_blot_path="All replicates - raw western blots/Rep3/Pgk1_for_Spn1_Rep3.tif",
                h3k36me3_blot_path="All replicates - raw western blots/Rep3/H3K36me3_Rep3.tif",
                h3k36me2_blot_path="All replicates - raw western blots/Rep3/H3K36me2_Rep3.tif",
                h3_blot_path="All replicates - raw western blots/Rep3/H3_for_H3K36me2_Rep3.tif",
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
                              label=c("", "H3K36me2", "H3K36me3", "H3", "Pgk1"),
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

    pgk1_image = readTIFF(pgk1_blot_path)
    pgk1_raster = rasterGrob(pgk1_image[175:(175 + blot_pixel_height),
                                        47:(47 + blot_pixel_width),
                                        1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[5])

    h3_image = readTIFF(h3_blot_path)
    h3_raster = rasterGrob(h3_image[200:(200 + blot_pixel_height),
                                    73:(73 + blot_pixel_width * 1.5 - 10),
                                    1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[4])

    h3k36me3_image = readTIFF(h3k36me3_blot_path)
    h3k36me3_raster = rasterGrob(h3k36me3_image[199:(199 + blot_pixel_height),
                                                52:(52 + blot_pixel_width * 1.5 - 10),
                                                1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[3])

    h3k36me2_image = readTIFF(h3k36me2_blot_path)
    h3k36me2_raster = rasterGrob(h3k36me2_image[197:(197 + blot_pixel_height),
                                                75:(75 + blot_pixel_width * 1.5 - 10),
                                                1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[2])


    western = gTree(children=gList(
        antigen_labels,
        spn1_labels,
        strain_lines,
        strain_labels,
        condition_labels,
        spn1_raster,
        h3_raster,
        h3k36me3_raster,
        h3k36me2_raster,
        pgk1_raster,
        blot_outlines#,
        # lane_alignment
        ))

    histone_spn1_depletion_westerns = as.ggplot(western) +
        labs(tag=panel_letter) +
        theme(plot.margin=margin(11/2, 0, 0, 11/2, "pt"),
              plot.tag=element_text(family="FreeSans",
                                    size=9,
                                    face="bold"))

    ggsave(pdf_out,
           plot=histone_spn1_depletion_westerns,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(histone_spn1_depletion_westerns,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     spn1_blot_path=snakemake@input[["spn1_blot"]],
     pgk1_blot_path=snakemake@input[["pgk1_blot"]],
     h3k36me3_blot_path=snakemake@input[["h3k36me3_blot"]],
     h3k36me2_blot_path=snakemake@input[["h3k36me2_blot"]],
     h3_blot_path=snakemake@input[["h3_blot"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])
