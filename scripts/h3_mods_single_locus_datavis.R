import = function(df,
                  data_path,
                  target_id){
    read_tsv(data_path,
                  col_names=c("group", "sample", "annotation",
                              "assay", "index", "position", "signal")) %>%
        group_by(group, assay, position) %>%
        summarize(signal = mean(signal)) %>%
        mutate(target=target_id) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(data_paths=c("ARO80_all-assays.tsv.gz",
                             "PHO89_all-assays.tsv.gz",
                             "UBI4_all-assays.tsv.gz",
                             "GCV3_all-assays.tsv.gz"),
                targets = c("ARO80", "PHO89", "UBI4", "GCV3"),
                transcript_annotation_path = "Scer_polIItranscripts-adjustedTSS.bed",
                orf_annotation_path = "Scer_nondubious_ORFs_and_blocked_reading_frames-adjustedATG.bed",
                theme_path = "spn1_2020_theme.R",
                fig_width=17.4,
                fig_height=7,
                panel_letter="b",
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)
    library(cowplot)

    df = tibble()
    for (i in 1:length(data_paths)){
        df %<>% import(data_path = data_paths[i],
                       target_id = targets[i])
    }

    df %<>%
        filter(assay %in% c("ChIPseq-H3",
                            "ChIPseq-H3K36me2-H3norm",
                            "ChIPseq-H3K36me3-H3norm",
                            "ChIPseq-H3K4me3-H3norm")) %>%
        mutate(group=ordered(group,
                             levels=c("non-depleted",
                                      "depleted"),
                             labels=c("non-depleted",
                                      "Spn1-depleted")),
               target=ordered(target,
                              levels=targets),
               assay=ordered(assay,
                             levels=c("ChIPseq-H3",
                                      "ChIPseq-H3K36me2-H3norm",
                                      "ChIPseq-H3K36me3-H3norm",
                                      "ChIPseq-H3K4me3-H3norm"),
                             labels=c("\"H3: log\"[2] ~ textstyle(frac(\"IP\", \"input\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"H3K36me2\", \"H3\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"H3K36me3\", \"H3\"))",
                                      "\"log\"[2] ~ textstyle(frac(\"H3K4me3\", \"H3\"))")))

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
        mutate(transcript_y = (scales::rescale(row_number() %% 2, to=c(1,2))),
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
                       x_max),
               orf_y=c(-0.4, 0.4)) %>%
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

    annotation_plot = ggplot() +
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
        geom_vline(data=x_limits %>%
                    pivot_longer(-target,
                                 names_to="type",
                                 values_to="x_intercept"),
                   aes(xintercept=x_intercept),
                   size=0.5,
                   color="gray70") +
        facet_grid(.~target,
                   scales="free",
                   space="free") +
        scale_x_continuous(expand=c(0,0),
                           name=NULL,
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        TRUE ~ paste(x, "kb"))) +
        scale_y_continuous(limits=range(orf_annotations["orf_y"])) +
        theme_default +
        theme(panel.grid=element_blank(),
              panel.spacing.x=unit(2, "pt"),
              axis.title.y=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_text(size=5,
                                       color="#ffffffff"),
              axis.ticks.length=unit(1, "pt"),
              axis.ticks=element_line(color="#ffffffff"),
              strip.text=element_blank(),
              plot.margin=margin(t=2, b=0, unit="pt"),
              panel.border=element_blank())

    histones_plot = ggplot() +
        geom_line(data=df,
                  aes(x=position,
                      y=signal,
                      color=group)) +
        facet_grid(assay ~ target,
                   scales="free",
                   space="free_x",
                   switch="y",
                   labeller=label_parsed) +
        scale_x_continuous(expand=c(0,0),
                           name=NULL,
                           breaks=scales::pretty_breaks(3),
                           labels=function(x) case_when(x==0 ~ "TSS",
                                                        TRUE ~ paste(x, "kb"))) +
        scale_y_continuous(breaks=scales::pretty_breaks(3),
                           name=NULL) +
        scale_color_viridis_d(end=0.6,
                              name=NULL,
                              guide=guide_legend(keywidth=unit(10, "pt"),
                                                 keyheight=unit(8, "pt"))) +
        theme_default +
        theme(strip.text.x=element_blank(),
              strip.text.y=element_text(angle=-180,
                                        vjust=0.5),
              strip.placement="outside",
              panel.spacing=unit(2, "pt"),
              panel.grid=element_blank(),
              legend.position=c(0.17, 0.84),
              legend.spacing.x=unit(1, "pt"),
              legend.spacing.y=unit(1, "pt"),
              legend.background=element_blank(),
              axis.text.y=element_text(size=5,
                                       hjust=1),
              axis.ticks.length=unit(1, "pt"),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5,
                                        hjust=1),
              plot.margin=margin(-11/4, 0, 0, 0, "pt"))

    h3_mods_single_locus = plot_grid(annotation_plot,
              histones_plot,
              align="v",
              axis="rl",
              ncol=1,
              rel_heights=c(0.1, 1))

    h3_mods_single_locus = h3_mods_single_locus +
        labs(tag=panel_letter) +
        theme(plot.tag=element_text(family="FreeSans",
                                    size=9,
                                    face="bold"),
              plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"))

    ggsave(pdf_out,
           h3_mods_single_locus,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(h3_mods_single_locus,
         file=grob_out)
}

main(data_paths=snakemake@input[["data_paths"]],
     targets=snakemake@params[["targets"]],
     transcript_annotation_path=snakemake@input[["transcript_annotations"]],
     orf_annotation_path=snakemake@input[["orf_annotations"]],
     theme_path=snakemake@input[["theme"]],
     fig_width = as.numeric(snakemake@params[["fig_width"]]),
     fig_height = as.numeric(snakemake@params[["fig_height"]]),
     panel_letter = snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

