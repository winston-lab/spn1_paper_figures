get_feature_stats = function(diffexp_results_path){
    diffexp_results = read_tsv(diffexp_results_path) %>%
        filter(! is.na(log2_foldchange))

    diffexp_results %>%
        pull(log2_foldchange) %>%
        set_names(diffexp_results %>%
                      pull(name)) %>%
        return()
}

get_fgsea_results = function(go_mapping_list,
                             feature_stats,
                             min_go_group_size,
                             max_go_group_size,
                             n_permutations){
    fgsea(pathways=go_mapping_list,
          stats=feature_stats,
          minSize=min_go_group_size,
          maxSize=max_go_group_size,
          nperm=n_permutations) %>%
        return()
}

main = function(theme_path="spn1_2020_theme.R",
                diffexp_results_path_single = "Spn1-IAA-v-Spn1-DMSO_rnaseq-spikenorm-verified-coding-genes-diffexp-results-all.tsv",
                go_mapping_path = "spn1_go_slim_mapping.tsv",
                min_go_group_size = 15,
                max_go_group_size = Inf,
                n_permutations = 1e4,
                pdf_out = "test.pdf",
                grob_out = "test.Rdata",
                fig_width=17.4/2,
                fig_height=4/11*12,
                panel_letter="c"){
    source(theme_path)
    library(fgsea)

    go_mapping = read_tsv(go_mapping_path)

    go_mapping_list = split(go_mapping[["feature_name"]],
                            go_mapping[["go_category"]])

    feature_stats_single = get_feature_stats(diffexp_results_path_single)

    fgsea_results_single = get_fgsea_results(go_mapping_list=go_mapping_list,
                                             feature_stats=feature_stats_single,
                                             min_go_group_size=min_go_group_size,
                                             max_go_group_size=max_go_group_size,
                                             n_permutations=n_permutations)

    padj = fgsea_results_single %>%
        as_tibble() %>%
        filter(pathway=="ESR_downregulated") %>%
        pull(padj)

    enrichment_plot_single = plotEnrichment(go_mapping_list[["ESR_downregulated"]],
                                            feature_stats_single)

    df = as_tibble(enrichment_plot_single$data)

    esr_fgsea = ggplot(data=df,
           aes(x=x,
               y=y)) +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_line(size=0.3) +
        geom_rug(sides="b",
                 size=0.1,
                 alpha=0.2,
                 length=unit(3, "pt")) +
        annotate(geom="text",
                 x=2500,
                 y=-0.18,
                 label=bquote("p"["adj"] == .(round(padj, 3))),
                 size=5/72*25.4,
                 family="FreeSans") +
        scale_x_continuous(name=quote("(more upregulated)" %<-%  "gene rank"  %->% "(more downregulated)"),
                           expand=c(0,0),
                           breaks=scales::pretty_breaks(4)) +
        scale_y_continuous(name="running\nenrichment\nscore",
                           expand=c(0.07, 0),
                           breaks=scales::pretty_breaks(4)) +
        labs(tag=panel_letter,
             title="Enrichment of genes downregulated in ESR") +
        theme_default +
        theme(panel.grid=element_blank(),
              plot.margin=margin(11/2, 11, 11/2, 11/2, "pt"),
              axis.text=element_text(size=5),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              plot.title=element_text(margin=margin(0, 0, 1, 0, "pt")))

    ggsave(pdf_out,
           plot=esr_fgsea,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(esr_fgsea,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     diffexp_results_path_single = snakemake@input[["diffexp_path_single"]],
     go_mapping_path = snakemake@input[["go_mapping_path"]],
     min_go_group_size = as.numeric(snakemake@params[["min_go_group_size"]]),
     max_go_group_size = as.numeric(snakemake@params[["max_go_group_size"]]),
     n_permutations = as.numeric(snakemake@params[["n_permutations"]]),
     pdf_out = snakemake@output[["pdf"]],
     grob_out = snakemake@output[["grob"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]])

