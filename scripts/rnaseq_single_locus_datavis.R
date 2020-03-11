import = function(df,
                  data_path,
                  target_id){
    read_tsv(data_path,
                  col_names=c("group", "sample", "annotation",
                              "assay", "index", "position", "signal")) %>%
        filter(assay %in% c("RNAseq-sense",
                            "RNAseq-antisense")) %>%
        group_by(group, assay, position) %>%
        summarize(signal = mean(signal)) %>%
        ungroup %>%
        mutate(assay = str_remove(assay, "RNAseq-"),
               target = target_id) %>%
        pivot_wider(names_from=assay,
                    values_from=signal) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(data_paths=c("SER3_all-assays.tsv.gz", "UBI4_all-assays.tsv.gz", "GCV3_all-assays.tsv.gz"),
                targets = c("SER3", "UBI4", "GCV3"),
                transcript_annotation_path = "Scer_polIItranscripts-adjustedTSS.bed",
                orf_annotation_path = "Scer_nondubious_ORFs_and_blocked_reading_frames-adjustedATG.bed",
                theme_path = "spn1_2020_theme.R",
                panel_letter = "b",
                fig_width=17.4/2,
                fig_height=7/12 * 12,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = tibble()
    for (i in 1:length(data_paths)){
        df %<>% import(data_path = data_paths[i],
                       target_id = targets[i])
    }

    df %<>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                             "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               target=ordered(target,
                              levels=targets))

    x_limits = df %>%
        group_by(target) %>%
        summarize(x_min = min(position),
                  x_max = max(position))

    transcript_annotations = read_tsv(transcript_annotation_path,
             col_names=c("chrom", "left", "right", "target", "score", "strand")) %>%
        filter(target %in% targets) %>%
        select(-score) %>%
        left_join(read_tsv(transcript_annotation_path,
                           col_names=c("chrom", "transcript_left", "transcript_right", "transcript_name",
                                       "score", "transcript_strand")) %>%
                      select(-score),
                  by=c("chrom")) %>%
        mutate_at(vars(transcript_left, transcript_right),
                  ~((. - ifelse(strand=="+", left, right)) * ifelse(strand=="+", 1, -1) / 1e3)) %>%
        mutate(t_left = ifelse(strand=="+", transcript_left, transcript_right),
               t_right = ifelse(strand=="+", transcript_right, transcript_left),
               transcript_left = t_left,
               transcript_right = t_right) %>%
        select(-c(t_left, t_right)) %>%
        left_join(x_limits,
                  by="target") %>%
        filter((transcript_left >= x_min & transcript_left <= x_max) |
                  (transcript_right >= x_min & transcript_left <= x_max) |
                   (transcript_left <= x_min & transcript_right >= x_max)) %>%
        left_join(read_tsv(orf_annotation_path,
                           col_names=c("chrom", "orf_left", "orf_right", "transcript_name",
                                       "score", "transcript_strand")) %>%
                      select(-score),
                  by=c("chrom", "transcript_name", "transcript_strand")) %>%
        mutate_at(vars(orf_left, orf_right),
                  ~((. - ifelse(strand=="+", left, right)) * ifelse(strand=="+", 1, -1) / 1e3)) %>%
        mutate(o_left = ifelse(strand=="+", orf_left, orf_right),
               o_right = ifelse(strand=="+", orf_right, orf_left),
               orf_left = o_left,
               orf_right = o_right) %>%
        select(-c(o_left, o_right)) %>%
        mutate(transcript_start = if_else(transcript_strand==strand,
                                          transcript_left,
                                          transcript_right),
               transcript_end = if_else(transcript_strand==strand,
                                        transcript_right,
                                        transcript_left),
               notch = ifelse(transcript_strand==strand,
                              (orf_right - orf_left) * 0.85 + orf_left,
                              orf_right - (orf_right - orf_left) * 0.85),
               orf_start = if_else(transcript_strand==strand,
                                   orf_left,
                                   orf_right),
               orf_end = if_else(transcript_strand==strand,
                                 orf_right,
                                 orf_left),
               label_xpos = if_else(is.na(notch),
                                    (transcript_start + transcript_end) / 2,
                                    (orf_start + orf_end) / 2)) %>%
        select(target,
               transcript_name,
               transcript_start,
               transcript_end,
               orf_start,
               orf_end,
               notch,
               label_xpos,
               x_min,
               x_max) %>%
        mutate_at(vars(transcript_start, transcript_end),
                  ~pmax(x_min, pmin(x_max, .))) %>%
        left_join(df %>%
                      group_by(target) %>%
                      summarize(y_range = max(sense) + max(antisense),
                                y_min = -max(antisense)),
                  by="target") %>%
        mutate(transcript_y = y_min - 0.14 * y_range * (scales::rescale(row_number() %% 2, to=c(1,2))),
               target=ordered(target,
                              levels=targets))

    orf_annotations = transcript_annotations %>%
        filter(! is.na(notch)) %>%
        expand(nesting(target,
                       transcript_name,
                       transcript_start,
                       transcript_end,
                       orf_start,
                       notch,
                       orf_end,
                       label_xpos,
                       x_min,
                       x_max,
                       y_range),
               orf_y=c(-1, 1)) %>%
        mutate(orf_y = orf_y * 0.08 * y_range) %>%
        pivot_longer(c(orf_start, notch, orf_end),
                     names_to="type",
                     values_to="orf_x") %>%
        mutate(orf_y = ifelse(type=="orf_end",
                              0,
                              orf_y)) %>%
        distinct() %>%
        group_by(transcript_name) %>%
        arrange(type,
                orf_y,
                .by_group=TRUE) %>%
        mutate(order = c(4,2,3,5,1)) %>%
        arrange(order,
                .by_group=TRUE) %>%
        mutate(orf_x = pmax(x_min, pmin(x_max, orf_x))) %>%
        left_join(transcript_annotations %>%
                      select(transcript_name, transcript_y),
                  by="transcript_name") %>%
        mutate(orf_y = orf_y + transcript_y,
               target=ordered(target,
                              levels=targets))

    rnaseq_single_locus_datavis = ggplot() +
        geom_ribbon(data=df,
                    aes(x=position,
                        ymin=0,
                        ymax=sense,
                        fill=group),
                    linetype="blank",
                    alpha=0.35) +
        geom_ribbon(data=df,
                    aes(x=position,
                        ymin=-antisense,
                        ymax=0,
                        fill=group),
                    linetype="blank",
                    alpha=0.35) +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_segment(data=transcript_annotations,
                     aes(x=transcript_start,
                         xend=transcript_end,
                         y=transcript_y,
                         yend=transcript_y),
                     color="gray50",
                     size=0.25) +
        geom_polygon(data=orf_annotations,
                     aes(x=orf_x,
                         y=orf_y,
                         group=transcript_name),
                     fill="grey80") +
        geom_text(data=transcript_annotations %>%
                      filter(label_xpos >= x_min &
                                 label_xpos <= x_max),
                  aes(x=label_xpos,
                      y=transcript_y,
                      label=if_else(str_detect(transcript_name, "stringtie"),
                                    "",
                                    transcript_name),
                      vjust=if_else(is.na(notch),
                                    -0.1, 0.5)),
                  size=5/72*25.4,
                  family="FreeSans",
                  fontface="italic") +
        facet_wrap(~target,
                   scales="free",
                   ncol=2)  +
        scale_x_continuous(expand=c(0,0),
                           name=NULL,
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        TRUE ~ paste(x, "kb"))) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           labels=function(x) abs(x),
                           name="normalized counts") +
        # scale_color_viridis_d(end=0.6,
        #                       name=NULL,
        #                       guide=guide_legend(keywidth=unit(10, "pt"),
        #                                          keyheight=unit(8, "pt"))) +
        scale_fill_viridis_d(end=0.6,
                              name=NULL,
                              guide=guide_legend(keywidth=unit(10, "pt"),
                                                 keyheight=unit(8, "pt"))) +
        labs(tag=panel_letter,
             title="RNA-seq") +
        theme_default +
        theme(strip.text=element_blank(),
              panel.spacing.y=unit(2, "pt"),
              panel.spacing.x=unit(3, "pt"),
              panel.grid=element_blank(),
              legend.position=c(0.01, 0.99),
              legend.justification = c(0, 1),
              legend.spacing.x=unit(1, "pt"),
              legend.spacing.y=unit(1, "pt"),
              legend.background=element_blank(),
              axis.text=element_text(size=5),
              axis.ticks.length=unit(1, "pt"))

    ggsave(pdf_out,
           rnaseq_single_locus_datavis,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(rnaseq_single_locus_datavis,
         file=grob_out)
}

main(data_paths=snakemake@input[["data_paths"]],
     targets=snakemake@params[["targets"]],
     transcript_annotation_path=snakemake@input[["transcript_annotations"]],
     orf_annotation_path=snakemake@input[["orf_annotations"]],
     theme_path=snakemake@input[["theme"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

