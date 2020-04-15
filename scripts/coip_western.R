
main = function(theme_path = "spn1_2020_theme.R",
                input_rpb3_blot_path = "Rpb3-Flag WB Rep1 inputs_10.tif",
                input_spt6_blot_path = "Spt6 WB Rep1 inputs_8.tif",
                input_spn1_blot_path = "Spn1 WB Rep1 inputs_8.5_2.tif",
                input_set2_blot_path = "Set2 WB Rep3 inputs_8-4.tif",
                ip_rpb3_blot_path = "Rpb3-Flag WB Rep1 IPs_9.tif",
                ip_spt6_blot_path = "Spt6 WB Rep1 IPs_7.5.tif",
                ip_spn1_blot_path = "Spn1 WB Rep1 IPs_8_3.tif",
                ip_set2_blot_path = "Set2 WB Rep2 IPs_8-7.tif",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=11.4,
                fig_height=13/3,
                panel_letter="a"){
    source(theme_path)
    library(grid)
    library(gridExtra)
    library(tiff)
    library(ggplotify)

    antigen_label_edge = 0.08
    y_margin_size = 0.015
    blot_max_y = 0.67
    blot_min_y = 0.045
    blot_height= ((blot_max_y - blot_min_y) - 3 * y_margin_size) / 6
    antigen_labels_y = c(blot_max_y - (1/2) * blot_height,
                         blot_max_y - (3/2) * blot_height - y_margin_size,
                         blot_max_y - (7/2) * blot_height - 2 * y_margin_size,
                         blot_max_y - (11/2) * blot_height - 3 * y_margin_size)

    blot_width = ((1 - antigen_label_edge) - (3 * 0.02)) / 2
    blot_centers = c(antigen_label_edge + 0.02 + blot_width / 2,
                     1 - 0.02 - blot_width / 2)
    input_lane_centers = c(0.125, 0.175, 0.223, 0.27,
                           0.315, 0.362, 0.408, 0.455, 0.505)
    ip_lane_centers = c(0.578, 0.63, 0.673, 0.720, 0.765,
                               0.81, 0.858, 0.905, 0.955)

    rpb3_centers = c(input_lane_centers[-1],
      ip_lane_centers[-1]) %>%
        matrix(ncol=2, byrow=TRUE) %>%
        rowMeans()

    strain_centers = rpb3_centers %>%
        matrix(ncol=2, byrow=TRUE) %>%
        rowMeans()

    antigen_labels = textGrob(x=antigen_label_edge,
                              y=antigen_labels_y,
                              label=c("Rpb3",
                                      "Spt6",
                                      "Spn1",
                                      "Set2"),
                              hjust=0.8,
                              gp=gpar(fontsize=7,
                                      fontfamily="FreeSans"))
    input_ip_labels = textGrob(x=blot_centers,
                               y=1,
                               label=c("input", "IP"),
                               vjust=1,
                               gp=gpar(fontsize=7,
                                       fontfamily="FreeSans"))
    input_ip_lines = segmentsGrob(x0=blot_centers - blot_width / 2,
                                  x1=blot_centers + blot_width / 2,
                                  y0=0.93, y1=0.93)
    strain_labels = textGrob(x=strain_centers,
                             y=0.91,
                             label=rep(c("Spn1-AID", "Spt6-AID"), 2),
                             vjust=1,
                             gp=gpar(fontsize=7,
                                     fontfamily="FreeSans"))
    strain_lines = segmentsGrob(x0=strain_centers - 0.08,
                                x1=strain_centers + 0.08,
                                y0=0.84, y1=0.84)
    rpb3_labels = textGrob(x=rpb3_centers,
                           y=0.81,
                           label=rep(c("untagged",
                                       "Rpb3-FLAG"),
                                     4),
                             vjust=1,
                             gp=gpar(fontsize=5,
                                     fontfamily="FreeSans"))
    rpb3_lines = segmentsGrob(x0=rpb3_centers - 0.035,
                              x1=rpb3_centers + 0.035,
                              y0=0.75, y1=0.75)
    rpb3_flag_labels = textGrob(x=c(input_lane_centers[1],
                                    ip_lane_centers[1]),
                                y=0.675,
                                label="Rpb3-FLAG",
                                rot=90,
                                hjust=0,
                                gp=gpar(fontsize=5,
                                        fontfamily="FreeSans"))
    dmso_iaa_labels = textGrob(x=c(input_lane_centers[-1],
                                   ip_lane_centers[-1]),
                               y=0.73,
                               label=c("D", "I"),
                               vjust=1,
                               gp=gpar(fontsize=7,
                                       fontfamily="FreeSans"))
    blot_outlines = rectGrob(width = blot_width,
                             height=rep(c(blot_height,
                                          blot_height,
                                          blot_height * 3,
                                          blot_height), 2),
                             x=c(rep(antigen_label_edge + 0.02 + blot_width / 2, 4),
                                 rep(1 - 0.02 - blot_width / 2, 4)),
                             y=rep(antigen_labels_y, 2),
                             gp=gpar(lwd=1,
                                     fill=NA))

    pixel_height = 30
    input_pixel_width = 365
    input_rpb3_image = readTIFF(input_rpb3_blot_path)
    input_rpb3_blot = rasterGrob(input_rpb3_image[220:(220+pixel_height),
                                                  139:(139+input_pixel_width - 10),
                                                  1:3],
                                 width=blot_width,
                                 height=blot_height,
                                 x=blot_centers[1],
                                 y=antigen_labels_y[1])
    input_spt6_image = readTIFF(input_spt6_blot_path)
    input_spt6_blot = rasterGrob(input_spt6_image[86:(86 + pixel_height),
                                                  141:(141 + input_pixel_width - 12),
                                                  1:3],
                                 width=blot_width,
                                 height=blot_height,
                                 x=blot_centers[1],
                                 y=antigen_labels_y[2])


    input_spn1_image = readTIFF(input_spn1_blot_path)
    input_spn1_blot = rasterGrob(input_spn1_image[150:(150 + pixel_height * 3),
                                                  137:(137 + input_pixel_width - 10),
                                                  1:3],
                                 width=blot_width,
                                 height=blot_height * 3,
                                 x=blot_centers[1],
                                 y=antigen_labels_y[3])

    input_set2_image = readTIFF(input_set2_blot_path)
    input_set2_blot = rasterGrob(input_set2_image[149:(149 + pixel_height),
                                                  118:(118 + input_pixel_width),
                                                  1:3],
                                 width=blot_width,
                                 height=blot_height,
                                 x=blot_centers[1],
                                 y=antigen_labels_y[4])

    ip_pixel_width = 379
    ip_rpb3_image = readTIFF(ip_rpb3_blot_path)
    ip_rpb3_blot = rasterGrob(ip_rpb3_image[221:(221+pixel_height),
                                            121:(121+ip_pixel_width),
                                            1:3],
                                 width=blot_width,
                                 height=blot_height,
                                 x=blot_centers[2],
                                 y=antigen_labels_y[1])
    ip_spt6_image = readTIFF(ip_spt6_blot_path)
    ip_spt6_blot = rasterGrob(ip_spt6_image[88:(88+pixel_height),
                                            129:(129+ip_pixel_width-16),
                                            1:3],
                                 width=blot_width,
                                 height=blot_height,
                                 x=blot_centers[2],
                                 y=antigen_labels_y[2])

    ip_spn1_image = readTIFF(ip_spn1_blot_path)
    ip_spn1_blot = rasterGrob(ip_spn1_image[159:(159+pixel_height*3),
                                            127:(127+ip_pixel_width-5),
                                            1:3],
                                 width=blot_width,
                                 height=blot_height*3,
                                 x=blot_centers[2],
                                 y=antigen_labels_y[3])

    ip_set2_image = readTIFF(ip_set2_blot_path)
    ip_set2_blot = rasterGrob(ip_set2_image[157:(157+pixel_height),
                                            135:(135+ip_pixel_width-12),
                                            1:3],
                                 width=blot_width,
                                 height=blot_height,
                                 x=blot_centers[2],
                                 y=antigen_labels_y[4])

    lane_labels = textGrob(x=c(input_lane_centers,
                               ip_lane_centers),
                           y=0,
                           label=seq(length(c(input_lane_centers,
                                              ip_lane_centers))),
                           vjust=0,
                           gp=gpar(fontsize=5,
                                   fontfamily="FreeSans"))
    lane_alignment = segmentsGrob(x0=c(input_lane_centers, ip_lane_centers),
                                  x1=c(input_lane_centers, ip_lane_centers),
                                  y0=0, y1=1,
                                  gp=gpar(lwd=0.2))
    horizontal_alignment = segmentsGrob(x0=0,
                                        x1=1,
                                        y0=c(0.62, 0.50, 0.38, 0.19, 0.06),
                                        y1=c(0.62, 0.50, 0.38, 0.19, 0.06),
                                        gp=gpar(lwd=0.2))

    coip_western = gTree(children=gList(
        antigen_labels,
        input_ip_labels,
        input_ip_lines,
        strain_labels,
        strain_lines,
        rpb3_labels,
        rpb3_lines,
        dmso_iaa_labels,
        rpb3_flag_labels,
        input_rpb3_blot,
        input_spt6_blot,
        input_spn1_blot,
        input_set2_blot,
        ip_rpb3_blot,
        ip_spt6_blot,
        ip_spn1_blot,
        ip_set2_blot,
        blot_outlines,
        lane_labels#,
        # horizontal_alignment
        # lane_alignment
        ))

    coip_western = as.ggplot(coip_western) +
        labs(tag=panel_letter) +
        theme(plot.margin=margin(0, 11/2, 0, 11/2, "pt"),
              plot.tag=element_text(family="FreeSans",
                                    size=9,
                                    face="bold"))

    ggsave(pdf_out,
           plot=coip_western,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(coip_western,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     input_rpb3_blot_path = snakemake@input[["input_rpb3_blot"]],
     input_spt6_blot_path = snakemake@input[["input_spt6_blot"]],
     input_spn1_blot_path = snakemake@input[["input_spn1_blot"]],
     input_set2_blot_path = snakemake@input[["input_set2_blot"]],
     ip_rpb3_blot_path = snakemake@input[["ip_rpb3_blot"]],
     ip_spt6_blot_path = snakemake@input[["ip_spt6_blot"]],
     ip_spn1_blot_path = snakemake@input[["ip_spn1_blot"]],
     ip_set2_blot_path = snakemake@input[["ip_set2_blot"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])
