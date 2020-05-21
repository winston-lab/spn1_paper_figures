
main = function(theme_path = "spn1_2020_theme.R",
                spn1_blot_path="All replicates - raw western blots/Rep1/Spn1_Rep1.tif",
                rpb1_blot_path="All replicates - raw western blots/Rep3/Rpb1-Pgk1_Rep3.tif",
                spt6_blot_path="All replicates - raw western blots/Rep3/Spt6_Rep3.tif",
                set2_blot_path="All replicates - raw western blots/Rep2/Set2_Rep2.tif",
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

    y_margin_size = 0.015
    blot_max_y = 0.88
    blot_min_y = 0.055
    blot_height = ((blot_max_y - blot_min_y) - 4 * y_margin_size) / 9.2
    blot_pixel_height = 30
    blot_pixel_width = 169
    antigen_labels_x = 0.75
    antigen_labels_y = blot_max_y - c((2.2 / 2) * blot_height,
                                      (2.2 + 0.5) * blot_height + y_margin_size,
                                      (2.2 + 1.5) * blot_height + 2 * y_margin_size,
                                      (2.2 + 2.5) * blot_height + 3 * y_margin_size,
                                      (2.2 + 3.5) * blot_height + 4 * y_margin_size,
                                      (2.2 + 4.5) * blot_height + 5 * y_margin_size,
                                      (2.2 + 5.5) * blot_height + 6 * y_margin_size,
                                      (2.2 + 6.5) * blot_height + 7 * y_margin_size)
    blot_width = antigen_labels_x - 0.02

    antigen_labels = textGrob(x=antigen_labels_x,
                              y=antigen_labels_y,
                              label=c("Spn1", "Rpb1", "Spt6", "Set2", "H3", "H3K36me3", "H3K36me2", "Pgk1"),
                              hjust=0,
                              gp=gpar(fontsize=7,
                                      fontfamily="FreeSans"))
    lane_centers = scales::rescale(seq(0.5/6, 5.5/6, length.out=6),
                                   from=c(0,1),
                                   to=c(0.02, 0.02 + blot_width))

    lane_alignment = segmentsGrob(x0=lane_centers,
                                  x1=lane_centers,
                                  y=0,
                                  y1=1)

    condition_labels = textGrob(x=lane_centers,
                                y=blot_max_y + y_margin_size,
                                label=rep(c("DMSO", "IAA"), 3),
                                vjust=0,
                                gp=gpar(fontsize=5,
                                        fontfamily="FreeSans"))

    strain_centers = lane_centers %>%
        matrix(ncol=2, byrow=TRUE) %>%
        rowMeans()

    strain_lines = segmentsGrob(x0=strain_centers - blot_width / 6 + 0.02,
                                x1=strain_centers + blot_width / 6 - 0.02,
                                y0=0.94,
                                y1=0.94)
    strain_labels = textGrob(x=strain_centers,
                             y=1,
                             label=c("wild type",
                                     "Spn1-noAID",
                                     "Spn1-AID"),
                             vjust=1,
                             gp=gpar(fontsize=7,
                                     fontfamily="FreeSans"))

    blot_outlines = rectGrob(width = blot_width,
                             height = c(blot_height * 2.2, rep(blot_height, 7)) ,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y,
                             gp=gpar(lwd=1,
                                     fill=NA))

    spn1_image = readTIFF(spn1_blot_path)
    spn1_raster = rasterGrob(spn1_image[100:(100 + 2.2 * blot_pixel_height),
                                        53:(53 + blot_pixel_width),
                                        1:3],
                             width=blot_width,
                             height=blot_height * 2.2,
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
                             y=antigen_labels_y[8])

    h3_image = readTIFF(h3_blot_path)
    h3_raster = rasterGrob(h3_image[200:(200 + blot_pixel_height),
                                    73:(73 + blot_pixel_width * 1.5 - 10),
                                    1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[5])

    h3k36me3_image = readTIFF(h3k36me3_blot_path)
    h3k36me3_raster = rasterGrob(h3k36me3_image[202:(202 + blot_pixel_height),
                                                52:(52 + blot_pixel_width * 1.5 - 10),
                                                1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[6])

    h3k36me2_image = readTIFF(h3k36me2_blot_path)
    h3k36me2_raster = rasterGrob(h3k36me2_image[194:(194 + blot_pixel_height),
                                                75:(75 + blot_pixel_width * 1.5 - 10),
                                                1:3],
                             width=blot_width,
                             height=blot_height,
                             x=antigen_labels_x / 2,
                             y=antigen_labels_y[7])


    western = gTree(children=gList(
        antigen_labels,
        strain_lines,
        strain_labels,
        condition_labels,
        spn1_raster,
        rpb1_raster,
        spt6_raster,
        set2_raster,
        h3_raster,
        h3k36me3_raster,
        h3k36me2_raster,
        pgk1_raster,
        blot_outlines#,
        # lane_alignment
        ))

    general_spn1_depletion_westerns = as.ggplot(western) +
        labs(tag=panel_letter) +
        theme(plot.margin=margin(11/2, 0, 0, 11/2, "pt"),
              plot.tag=element_text(family="FreeSans",
                                    size=9,
                                    face="bold"))

    ggsave(pdf_out,
           plot=general_spn1_depletion_westerns,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(general_spn1_depletion_westerns,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     spn1_blot_path=snakemake@input[["spn1_blot"]],
     rpb1_blot_path=snakemake@input[["rpb1_blot"]],
     spt6_blot_path=snakemake@input[["spt6_blot"]],
     set2_blot_path=snakemake@input[["set2_blot"]],
     pgk1_blot_path=snakemake@input[["pgk1_blot"]],
     h3k36me3_blot_path=snakemake@input[["h3k36me3_blot"]],
     h3k36me2_blot_path=snakemake@input[["h3k36me2_blot"]],
     h3_blot_path=snakemake@input[["h3_blot"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])

