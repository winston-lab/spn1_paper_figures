import = function(df,
                  data_path){
    read_tsv(data_path,
             col_names=c("group",
                         "sample",
                         "annotation",
                         "assay",
                         "index",
                         "position",
                         "signal")) %>%
        group_by(group, assay, annotation, index, position) %>%
        summarize(signal=mean(signal, na.rm=TRUE)) %>%
        pivot_wider(names_from=group,
                    values_from=signal) %>%
        mutate(ratio = if_else(str_detect(assay, "RNAseq"),
                               log2((depleted + 0.01) / (`non-depleted` + 0.01)),
                               depleted - `non-depleted`)) %>%
        bind_rows(df, .) %>%
        return()
}

import_bed = function(bed_path,
                      annotation_id){
    read_tsv(bed_path,
             col_names=c("chrom", "start", "end",
                         "name", "score", "strand")) %>%
        mutate(annotation=annotation_id,
               index=row_number()) %>%
        return()
}

main = function(theme_path = "spn1_2020_theme.R",
                data_paths = c("test_RNAseq-sense.tsv.gz",
                               "test_ChIPseq-H3K36me2-H3norm.tsv.gz"),
                rp_with_intron_bed = "Scer_RPgene_with_intron_transcripts_TSS_intron_dist_sort.bed",
                nonrp_with_intron_bed = "Scer_nonRPgene_with_intron_transcripts_TSS_intron_dist_sort.bed",
                intron_bed = "Scer_introns_in_coding_genes.bed",
                name_lookup = "Scer-sys-to-common-name.tsv",
                panel_letter="x",
                fig_width=17.4/2,
                fig_height=8,
                pdf_out = "test.pdf",
                grob_out = "test.Rdata"){
    source(theme_path)

    df = tibble()
    for (data_path in data_paths){
        df %<>%
            import(data_path)
    }
    df %<>%
        mutate(assay=ordered(assay,
                             levels=c("RNAseq-sense",
                                      "ChIPseq-H3K36me2-H3norm"),
                             labels=c("sense RNA-seq",
                                      # "\"log\"[2] ~ textstyle(frac(\"H3K36me2\", \"H3\"))")),
                                      "H3K36me2 / H3")),
               annotation=ordered(annotation,
                                  levels=c("RP genes with introns",
                                           "nonRP genes with introns"),
                                  labels=c("89 RP genes with introns",
                                           "153 non-RP genes with introns")))

    df_annotation = import_bed(rp_with_intron_bed,
                    "RP genes with introns") %>%
        bind_rows(import_bed(nonrp_with_intron_bed,
                             "nonRP genes with introns")) %>%
        left_join(import_bed(intron_bed,
                             "intron") %>%
                      separate(name,
                               into=c("systematic", "junk"),
                               sep="_") %>%
                      select(-c(score,junk,annotation,index)) %>%
                      left_join(read_tsv(name_lookup,
                                         col_names=c("systematic", "short")),
                                by="systematic"),
                  by=c("chrom", "name"="short", "strand"),
                  suffix=c("_transcript",
                           "_intron")) %>%
        mutate(start = if_else(strand=="+",
                               (start_intron - start_transcript) / 1e3,
                               (end_transcript - end_intron) / 1e3),
               end = if_else(strand=="+",
                             (end_intron - start_transcript) / 1e3,
                             (end_transcript - start_intron)/ 1e3),
               annotation=ordered(annotation,
                                  levels=c("RP genes with introns",
                                           "nonRP genes with introns"),
                                  labels=c("89 RP genes with introns",
                                           "153 non-RP genes with introns")))

    intron_containing_gene_heatmaps = ggplot() +
        geom_raster(data=df,
                    aes(x=position,
                        y=index,
                        fill=ratio),
                    interpolate=FALSE) +
        geom_point(data=df_annotation,
                   aes(x=start,
                       y=index),
                   color="firebrick1",
                   size=0.3,
                   alpha=0.8,
                   shape=16) +
        geom_point(data=df_annotation,
                   aes(x=end,
                       y=index),
                   color="orange",
                   size=0.3,
                   alpha=0.6,
                   shape=16) +
        scale_fill_distiller(palette="PRGn",
                             direction=1,
                             # limits=c(-1.5, 1.5),
                             limits=c(-1.7, 1.7),
                             oob=scales::squish,
                             guide=guide_colorbar(barheight=unit(0.2, "cm"),
                                                  barwidth=unit(5, "cm")),
                             breaks=scales::pretty_breaks(5),
                             name=quote("log"[2] ~ textstyle(frac("Spn1-depleted", "non-depleted")))) +
        facet_grid(annotation~assay,
                   scales="free_y",
                   space="free_y",
                   switch="y") +
        scale_x_continuous(limits=c(-0.1, 1.2),
                           expand=c(0, 0),
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) if_else(x==0, "TSS", paste(x, "kb"))) +
        scale_y_reverse(expand=c(0,0)) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              panel.spacing.y=unit(2, "pt"),
              panel.spacing.x=unit(4, "pt"),
              strip.text.x=element_text(margin=margin(0,0,1,0,"pt")),
              strip.text.y=element_text(margin=margin(0,1,0,0,"pt")),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title=element_blank(),
              strip.placement="outside",
              legend.position="top",
              legend.margin=margin(0,0,-2,0,"pt"),
              plot.margin=margin(11/2, 11/2, 0, 11/2, "pt"))

    ggsave(pdf_out,
           plot=intron_containing_gene_heatmaps,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(intron_containing_gene_heatmaps,
         file=grob_out)
}

main(theme_path = snakemake@input[["theme"]],
     data_paths = snakemake@input[["data"]],
     rp_with_intron_bed = snakemake@input[["rp_with_intron_bed"]],
     nonrp_with_intron_bed = snakemake@input[["nonrp_with_intron_bed"]],
     intron_bed = snakemake@input[["intron_bed"]],
     name_lookup = snakemake@input[["name_lookup"]],
     panel_letter = snakemake@params[["panel_letter"]],
     fig_width = snakemake@params[["fig_width"]],
     fig_height = snakemake@params[["fig_height"]],
     pdf_out = snakemake@output[["pdf"]],
     grob_out = snakemake@output[["grob"]])

