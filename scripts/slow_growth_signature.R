main = function(theme_path = "spn1_2020_theme.R",
                oduibhir_supp4="oduibhir2014_supp_data_4.txt",
                rnaseq_results="Spn1-IAA-v-Spn1-DMSO_rnaseq-spikenorm-verified-coding-genes-diffexp-results-all.tsv",
                sys_to_common="Scer-sys-to-common-name.tsv",
                pdf_out="test.pdf",
                grob_out="test.Rdata",
                fig_width=8.5,
                fig_height=8.5*9/16,
                panel_letter="x"){
    source(theme_path)

    df = read_tsv(oduibhir_supp4,
                  skip=1,
                  col_names=c("systematic",
                              "short",
                              "lfc_slowgrowth")) %>%
        left_join(read_tsv(rnaseq_results) %>%
                      left_join(read_tsv(sys_to_common,
                                         col_names=c("systematic", "name")),
                                by="name"),
                  by="systematic")

    # pearson = cor(df[["log2_foldchange"]],
    #               df[["lfc_slowgrowth"]],
    #               use="complete.obs")

    slow_growth_signature = ggplot(data=df,
           aes(x=log2_foldchange,
               y=lfc_slowgrowth)) +
        geom_hline(yintercept=0,
                   color="gray70",
                   size=0.2) +
        geom_vline(xintercept=0,
                   color="gray70",
                   size=0.2) +
        # geom_abline(slope=1,
        #             intercept=0,
        #             color="gray70",
        #             size=0.2) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=150,
                    shape=16,
                    size=0.3,
                    alpha=0.8) +
        # annotate(geom="text",
        #          x=min(df[["log2_foldchange"]], na.rm=TRUE),
        #          y=max(df[["lfc_slowgrowth"]], na.rm=TRUE),
        #          label=paste0("r=", round(pearson, 2)),
        #          vjust=1,
        #          hjust=0,
        #          size=8/72*25.4) +
        scale_color_viridis_c(option="cividis") +
        scale_x_continuous(breaks=scales::pretty_breaks(4),
                           name=expression("RNA-seq: log"[2] ~ textstyle(frac("Spn1-depleted", "non-depleted"))))+
        scale_y_continuous(breaks=scales::pretty_breaks(4),
                           name=expression(atop("slow growth signature:", "log"[2] ~
                                               textstyle(frac("deletion mutant", "wild type"))))) +
        labs(tag=panel_letter) +
        theme_default +
        theme(legend.position="none",
              panel.grid=element_blank(),
              axis.text=element_text(color="black"),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5,
                                        hjust=1))

    ggsave(pdf_out,
           plot=slow_growth_signature,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    save(slow_growth_signature,
         file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     oduibhir_supp4=snakemake@input[["oduibhir_supp4"]],
     rnaseq_results=snakemake@input[["rnaseq_results"]],
     sys_to_common=snakemake@input[["sys_to_common"]],
     fig_width=snakemake@params[["fig_width"]],
     fig_height=snakemake@params[["fig_height"]],
     panel_letter=snakemake@params[["panel_letter"]],
     pdf_out=snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])

