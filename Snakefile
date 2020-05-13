#!/usr/bin/env python

configfile: "config.yaml"

rule target:
    input:
        "panels/spn1_depletion_western.pdf",
        "panels/spn1_depletion_scatter.pdf",
        "panels/spn1_depletion_metagene.pdf",
        "panels/spn1_depletion_chipseq_barplot.pdf",
        "figures/figure_1_spn1_depletion.pdf",
        "panels/spn1_depletion_viability.pdf",
        "panels/spn1_v_rpb1_nondepleted.pdf",
        "panels/spn1_rpb1norm_v_rpb1_nondepleted.pdf",
        "figures/figure_S1_spn1_depletion_supplemental.pdf",
        "panels/rnaseq_maplot.pdf",
        "panels/rnaseq_maplot_alternate.pdf",
        "panels/rnaseq_single_locus_datavis.pdf",
        "panels/rnaseq_vs_rpb1_single_locus.pdf",
        "panels/rpb1_metagenes.pdf",
        "panels/rnaseq_vs_rpb1.pdf",
        "figures/figure_2_rnaseq_rpb1.pdf",
        "panels/rnaseq_single_locus_datavis_supp.pdf",
        "panels/differential_expression_rtqpcr.pdf",
        "panels/rna_single_v_custom.pdf",
        "panels/antisense_single_locus_datavis.pdf",
        "panels/chipseq_abundance_barplots_rpb1.pdf",
        "panels/esr_fgsea.pdf",
        "panels/slow_growth_signature.pdf",
        "figures/figure_S2_rnaseq_rpb1_supplemental.pdf",
        "panels/splicing.pdf",
        "panels/splicing_rtqpcr.pdf",
        "figures/figure_7_splicing.pdf",
        "panels/promoter_swap_diagram.pdf",
        "panels/promoter_swap_rtqpcr.pdf",
        "figures/figure_3_promoter_swap.pdf",
        "panels/coip_western.pdf",
        "panels/set2_metagene.pdf",
        "panels/spt6_metagene.pdf",
        "panels/set2_abundance_chipseq_barplot.pdf",
        "panels/spt6_abundance_chipseq_barplot.pdf",
        expand("panels/{factor}_chip_maplot.pdf", factor=["Set2", "Spt6", "Set2-Rpb1norm", "Spt6-Rpb1norm"]),
        "figures/figure_4_set2_spt6.pdf",
        "panels/set2_spt6_v_rpb1.pdf",
        "panels/set2_spt6_rpb1norm_v_rpb1.pdf",
        "figures/figure_S3_set2_spt6_supplemental.pdf",
        "panels/h3_vs_rpb1_ma.pdf",
        "panels/reduced_h3_h3_metagene.pdf",
        "panels/h3_single_locus_datavis.pdf",
        "figures/figure_5_h3.pdf",
        "panels/h3_vs_rpb1_ma_5p500bp.pdf",
        "panels/h3_plmin_reduced_lfc_scatterplots.pdf",
        "panels/reduced_h3_matched_metagenes.pdf",
        "figures/figure_S4_h3_supplemental.pdf",
        "panels/h3_metagene.pdf",
        "panels/h3_modification_datavis_ratio.pdf",
        "panels/h3_mods_single_locus_datavis.pdf",
        "figures/figure_6_h3_mods.pdf",
        "panels/chipseq_abundance_barplots_h3.pdf",
        "panels/h3_mods_non_h3_norm.pdf",
        "panels/h3_mods_facet_expression.pdf",
        "figures/figure_S5_h3_mods_supplemental.pdf",
        expand("panels/{mod}_facet_expression_length.pdf", mod=["H3K4me3", "H3K36me2", "H3K36me3", "H3"]),
        "panels/rpgene_datavis.pdf",
        "panels/rpgene_datavis_supp.pdf",
        # "panels/rpgene_datavis_length_filtered.pdf",
        # "panels/rpgene_datavis_supp_length_filtered.pdf",
        "figures/figure_S6_splicing_supplemental.pdf",
        "panels/intron_containing_gene_heatmaps.pdf",

rule register_fonts:
    input:
        fonts_path = config["fonts_path"],
    output:
        output_path = ".fonts_registered.txt"
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/register_fonts.R"

rule spn1_depletion_western:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        spn1_blot = config["spn1_depletion_western"]["spn1_blot"],
        pgk1_blot = config["spn1_depletion_western"]["pgk1_blot"],
        quant_data = config["spn1_depletion_western"]["quant_data"],
    output:
        pdf = "panels/spn1_depletion_western.pdf",
        grob = "panels/spn1_depletion_western.Rdata",
    params:
        fig_height = eval(str(config["spn1_depletion_western"]["fig_height"])),
        fig_width = eval(str(config["spn1_depletion_western"]["fig_width"])),
        panel_letter = config["spn1_depletion_western"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_depletion_western.R"

rule spn1_depletion_scatter:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["spn1_depletion_scatter"]["data"],
    output:
        pdf = "panels/spn1_depletion_scatter.pdf",
        grob = "panels/spn1_depletion_scatter.Rdata",
    params:
        fig_height = eval(str(config["spn1_depletion_scatter"]["fig_height"])),
        fig_width = eval(str(config["spn1_depletion_scatter"]["fig_width"])),
        panel_letter = config["spn1_depletion_scatter"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_depletion_scatter.R"

rule spn1_depletion_metagene:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["spn1_depletion_metagene"]["data"],
    output:
        pdf = "panels/spn1_depletion_metagene.pdf",
        grob = "panels/spn1_depletion_metagene.Rdata",
    params:
        fig_height = eval(str(config["spn1_depletion_metagene"]["fig_height"])),
        fig_width = eval(str(config["spn1_depletion_metagene"]["fig_width"])),
        panel_letter = config["spn1_depletion_metagene"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_depletion_metagene.R"

rule spn1_depletion_chipseq_barplot:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["spn1_depletion_chipseq_barplot"]["data"],
    output:
        pdf = "panels/spn1_depletion_chipseq_barplot.pdf",
        grob = "panels/spn1_depletion_chipseq_barplot.Rdata",
    params:
        fig_height = eval(str(config["spn1_depletion_chipseq_barplot"]["fig_height"])),
        fig_width = eval(str(config["spn1_depletion_chipseq_barplot"]["fig_width"])),
        panel_letter = config["spn1_depletion_chipseq_barplot"]["panel_letter"],
        plot_title = "Spn1 ChIP-seq"
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/chipseq_abundance_barplot.R"

rule assemble_figure_spn1_depletion:
    input:
        fonts = ".fonts_registered.txt",
        spn1_depletion_western = "panels/spn1_depletion_western.Rdata",
        spn1_depletion_chipseq_barplot = "panels/spn1_depletion_chipseq_barplot.Rdata",
        spn1_depletion_metagene = "panels/spn1_depletion_metagene.Rdata",
        spn1_depletion_scatter = "panels/spn1_depletion_scatter.Rdata",
    output:
        pdf = "figures/figure_1_spn1_depletion.pdf"
    params:
        fig_width = eval(str(config["spn1_depletion_figure"]["fig_width"])),
        fig_height = eval(str(config["spn1_depletion_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_spn1_depletion.R"

rule spn1_depletion_viability:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["spn1_depletion_viability"]["data"],
    output:
        pdf = "panels/spn1_depletion_viability.pdf",
        grob = "panels/spn1_depletion_viability.Rdata",
    params:
        fig_height = eval(str(config["spn1_depletion_viability"]["fig_height"])),
        fig_width = eval(str(config["spn1_depletion_viability"]["fig_width"])),
        panel_letter = config["spn1_depletion_viability"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_depletion_viability.R"

rule spn1_v_rpb1:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        spn1 = config["spn1_v_rpb1"]["spn1"],
        rpb1 = config["spn1_v_rpb1"]["rpb1"],
    output:
        pdf = "panels/spn1_v_rpb1.pdf",
        grob = "panels/spn1_v_rpb1.Rdata",
    params:
        fig_height = eval(str(config["spn1_v_rpb1"]["fig_height"])),
        fig_width = eval(str(config["spn1_v_rpb1"]["fig_width"])),
        panel_letter = config["spn1_v_rpb1"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_v_rpb1.R"

rule spn1_rpb1norm_v_rpb1_nondepleted:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        spn1 = config["spn1_rpb1norm_v_rpb1_nondepleted"]["spn1"],
        rpb1 = config["spn1_rpb1norm_v_rpb1_nondepleted"]["rpb1"],
    output:
        pdf = "panels/spn1_rpb1norm_v_rpb1_nondepleted.pdf",
        grob = "panels/spn1_rpb1norm_v_rpb1_nondepleted.Rdata",
    params:
        fig_height = eval(str(config["spn1_rpb1norm_v_rpb1_nondepleted"]["fig_height"])),
        fig_width = eval(str(config["spn1_rpb1norm_v_rpb1_nondepleted"]["fig_width"])),
        panel_letter = config["spn1_rpb1norm_v_rpb1_nondepleted"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_rpb1norm_v_rpb1_nondepleted.R"

rule assemble_figure_spn1_depletion_supp:
    input:
        fonts = ".fonts_registered.txt",
        spn1_depletion_viability = "panels/spn1_depletion_viability.Rdata",
        spn1_v_rpb1 = "panels/spn1_v_rpb1.Rdata",
        spn1_rpb1norm_v_rpb1_nondepleted = "panels/spn1_rpb1norm_v_rpb1_nondepleted.Rdata",
    output:
        pdf = "figures/figure_S1_spn1_depletion_supplemental.pdf"
    params:
        fig_width = eval(str(config["spn1_depletion_supplemental"]["fig_width"])),
        fig_height = eval(str(config["spn1_depletion_supplemental"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_spn1_depletion_supp.R"


rule rnaseq_maplot:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["rnaseq_maplot"]["data"]
    output:
        pdf = "panels/rnaseq_maplot.pdf",
        grob = "panels/rnaseq_maplot.Rdata",
    params:
        fig_height = eval(str(config["rnaseq_maplot"]["fig_height"])),
        fig_width = eval(str(config["rnaseq_maplot"]["fig_width"])),
        panel_letter = config["rnaseq_maplot"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rnaseq_maplot.R"

rule rnaseq_maplot_alternate:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["rnaseq_maplot"]["data"]
    output:
        pdf = "panels/rnaseq_maplot_alternate.pdf",
        grob = "panels/rnaseq_maplot_alternate.Rdata",
    params:
        fig_height = eval(str(config["rnaseq_maplot"]["fig_height"])),
        fig_width = eval(str(config["rnaseq_maplot"]["fig_width"])),
        panel_letter = config["rnaseq_maplot"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rnaseq_maplot_alternate.R"

rule rnaseq_single_locus_datavis:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data_paths = list(config["rnaseq_single_locus_datavis"]["data"].values()),
        transcript_annotations = config["rnaseq_single_locus_datavis"]["transcript_annotation"],
        orf_annotations = config["rnaseq_single_locus_datavis"]["orf_annotation"]
    output:
        pdf = "panels/rnaseq_single_locus_datavis.pdf",
        grob = "panels/rnaseq_single_locus_datavis.Rdata",
    params:
        targets = list(config["rnaseq_single_locus_datavis"]["data"].keys()),
        fig_height = eval(str(config["rnaseq_single_locus_datavis"]["fig_height"])),
        fig_width = eval(str(config["rnaseq_single_locus_datavis"]["fig_width"])),
        panel_letter = config["rnaseq_single_locus_datavis"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rnaseq_single_locus_datavis.R"

rule rnaseq_vs_rpb1_single_locus:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data_paths = list(config["rnaseq_vs_rpb1_single_locus"]["data"].values()),
        transcript_annotations = config["rnaseq_vs_rpb1_single_locus"]["transcript_annotation"],
        orf_annotations = config["rnaseq_vs_rpb1_single_locus"]["orf_annotation"]
    output:
        pdf = "panels/rnaseq_vs_rpb1_single_locus.pdf",
        grob = "panels/rnaseq_vs_rpb1_single_locus.Rdata",
    params:
        targets = list(config["rnaseq_vs_rpb1_single_locus"]["data"].keys()),
        fig_height = eval(str(config["rnaseq_vs_rpb1_single_locus"]["fig_height"])),
        fig_width = eval(str(config["rnaseq_vs_rpb1_single_locus"]["fig_width"])),
        panel_letter = config["rnaseq_vs_rpb1_single_locus"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rnaseq_vs_rpb1_single_locus.R"

rule rpb1_metagenes:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        rpb1 = config["rpb1_metagenes"]["rpb1"],
        modifications = config["rpb1_metagenes"]["modifications"],
    output:
        pdf = "panels/rpb1_metagenes.pdf",
        grob = "panels/rpb1_metagenes.Rdata",
    params:
        fig_height = eval(str(config["rpb1_metagenes"]["fig_height"])),
        fig_width = eval(str(config["rpb1_metagenes"]["fig_width"])),
        panel_letter = config["rpb1_metagenes"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rpb1_metagenes.R"

rule rnaseq_vs_rpb1:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        rnaseq_path = config["rnaseq_vs_rpb1"]["rnaseq"],
        rpb1_path = config["rnaseq_vs_rpb1"]["rpb1"]
    output:
        pdf = "panels/rnaseq_vs_rpb1.pdf",
        grob = "panels/rnaseq_vs_rpb1.Rdata",
    params:
        fig_height = eval(str(config["rnaseq_vs_rpb1"]["fig_height"])),
        fig_width = eval(str(config["rnaseq_vs_rpb1"]["fig_width"])),
        panel_letter = config["rnaseq_vs_rpb1"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rnaseq_vs_rpb1.R"

rule assemble_figure_rnaseq_rpb1:
    input:
        fonts = ".fonts_registered.txt",
        rnaseq_maplot = "panels/rnaseq_maplot.Rdata",
        # rnaseq_single_locus_datavis = "panels/rnaseq_single_locus_datavis.Rdata",
        rnaseq_vs_rpb1_single_locus = "panels/rnaseq_vs_rpb1_single_locus.Rdata",
        rpb1_metagenes = "panels/rpb1_metagenes.Rdata",
        rnaseq_vs_rpb1 = "panels/rnaseq_vs_rpb1.Rdata",
    output:
        pdf = "figures/figure_2_rnaseq_rpb1.pdf"
    params:
        fig_width = eval(str(config["rnaseq_rpb1_figure"]["fig_width"])),
        fig_height = eval(str(config["rnaseq_rpb1_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_rnaseq_rpb1.R"

rule rnaseq_single_locus_datavis_supp:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data_paths = list(config["rnaseq_single_locus_datavis_supp"]["data"].values()),
        transcript_annotations = config["rnaseq_single_locus_datavis_supp"]["transcript_annotation"],
        orf_annotations = config["rnaseq_single_locus_datavis_supp"]["orf_annotation"]
    output:
        pdf = "panels/rnaseq_single_locus_datavis_supp.pdf",
        grob = "panels/rnaseq_single_locus_datavis_supp.Rdata",
    params:
        targets = list(config["rnaseq_single_locus_datavis_supp"]["data"].keys()),
        fig_height = eval(str(config["rnaseq_single_locus_datavis_supp"]["fig_height"])),
        fig_width = eval(str(config["rnaseq_single_locus_datavis_supp"]["fig_width"])),
        panel_letter = config["rnaseq_single_locus_datavis_supp"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rnaseq_single_locus_datavis.R"

rule differential_expression_rtqpcr:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["differential_expression_rtqpcr"]["data"],
        rnaseq = config["differential_expression_rtqpcr"]["rnaseq"]
    output:
        pdf = "panels/differential_expression_rtqpcr.pdf",
        grob = "panels/differential_expression_rtqpcr.Rdata",
    params:
        fig_height = eval(str(config["differential_expression_rtqpcr"]["fig_height"])),
        fig_width = eval(str(config["differential_expression_rtqpcr"]["fig_width"])),
        panel_letter = config["differential_expression_rtqpcr"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/differential_expression_rtqpcr.R"

rule rna_single_v_custom:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        single = config["rna_single_v_custom"]["single"],
        custom = config["rna_single_v_custom"]["custom"]
    output:
        pdf = "panels/rna_single_v_custom.pdf",
        grob = "panels/rna_single_v_custom.Rdata",
    params:
        fig_height = eval(str(config["rna_single_v_custom"]["fig_height"])),
        fig_width = eval(str(config["rna_single_v_custom"]["fig_width"])),
        panel_letter = config["rna_single_v_custom"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rna_single_v_custom.R"


rule antisense_single_locus_datavis:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data_paths = list(config["antisense_single_locus_datavis"]["data"].values()),
        transcript_annotations = config["antisense_single_locus_datavis"]["transcript_annotation"],
        orf_annotations = config["antisense_single_locus_datavis"]["orf_annotation"]
    output:
        pdf = "panels/antisense_single_locus_datavis.pdf",
        grob = "panels/antisense_single_locus_datavis.Rdata",
    params:
        targets = list(config["antisense_single_locus_datavis"]["data"].keys()),
        fig_height = eval(str(config["antisense_single_locus_datavis"]["fig_height"])),
        fig_width = eval(str(config["antisense_single_locus_datavis"]["fig_width"])),
        panel_letter = config["antisense_single_locus_datavis"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rnaseq_single_locus_datavis.R"

rule chipseq_abundance_barplots_rpb1:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = list(config["chipseq_abundance_barplots_rpb1"]["data"].values()),
    output:
        pdf = "panels/chipseq_abundance_barplots_rpb1.pdf",
        grob = "panels/chipseq_abundance_barplots_rpb1.Rdata",
    params:
        factor_ids = list(config["chipseq_abundance_barplots_rpb1"]["data"].keys()),
        fig_height = eval(str(config["chipseq_abundance_barplots_rpb1"]["fig_height"])),
        fig_width = eval(str(config["chipseq_abundance_barplots_rpb1"]["fig_width"])),
        panel_letter = config["chipseq_abundance_barplots_rpb1"]["panel_letter"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/chipseq_abundance_barplots_multifactor.R"

rule esr_fgsea:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        diffexp_path_single = config["esr_fgsea"]["diffexp_path_single"],
        go_mapping_path = config["esr_fgsea"]["go_mapping_path"],
    output:
        pdf = "panels/esr_fgsea.pdf",
        grob = "panels/esr_fgsea.Rdata",
    params:
        min_go_group_size=15,
        max_go_group_size="Inf",
        n_permutations=1e5,
        fig_height = eval(str(config["esr_fgsea"]["fig_height"])),
        fig_width = eval(str(config["esr_fgsea"]["fig_width"])),
        panel_letter = config["esr_fgsea"]["panel_letter"],
    conda:
        "envs/fgsea.yaml"
    script:
        "scripts/esr_fgsea.R"

rule slow_growth_signature:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        oduibhir_supp4 = config["slow_growth_signature"]["oduibhir"],
        rnaseq_results = config["slow_growth_signature"]["rnaseq"],
        sys_to_common = config["slow_growth_signature"]["sys_to_common"]
    output:
        pdf = "panels/slow_growth_signature.pdf",
        grob = "panels/slow_growth_signature.Rdata",
    params:
        fig_height = eval(str(config["slow_growth_signature"]["fig_height"])),
        fig_width = eval(str(config["slow_growth_signature"]["fig_width"])),
        panel_letter = config["slow_growth_signature"]["panel_letter"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/slow_growth_signature.R"

rule assemble_figure_rnaseq_rpb1_supp:
    input:
        fonts = ".fonts_registered.txt",
        differential_expression_rtqpcr = "panels/differential_expression_rtqpcr.Rdata",
        rna_single_v_custom = "panels/rna_single_v_custom.Rdata",
        esr_fgsea = "panels/esr_fgsea.Rdata",
        slow_growth_signature = "panels/slow_growth_signature.Rdata",
        antisense_single_locus_datavis = "panels/antisense_single_locus_datavis.Rdata",
        chipseq_abundance_barplots_rpb1 = "panels/chipseq_abundance_barplots_rpb1.Rdata",
    output:
        pdf = "figures/figure_S2_rnaseq_rpb1_supplemental.pdf"
    params:
        fig_width = eval(str(config["rnaseq_rpb1_supplemental"]["fig_width"])),
        fig_height = eval(str(config["rnaseq_rpb1_supplemental"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_rnaseq_rpb1_supp.R"

rule splicing:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["splicing"]["data"],
        aliases = config["splicing"]["aliases"],
        rp_genes = config["splicing"]["rp_genes"],
    output:
        pdf = "panels/splicing.pdf",
        grob = "panels/splicing.Rdata",
    params:
        fig_height = eval(str(config["splicing"]["fig_height"])),
        fig_width = eval(str(config["splicing"]["fig_width"])),
        panel_letter = config["splicing"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/splicing.R"

rule splicing_rtqpcr:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["splicing_rtqpcr"]["data"],
    output:
        pdf = "panels/splicing_rtqpcr.pdf",
        grob = "panels/splicing_rtqpcr.Rdata",
    params:
        fig_height = eval(str(config["splicing_rtqpcr"]["fig_height"])),
        fig_width = eval(str(config["splicing_rtqpcr"]["fig_width"])),
        panel_letter = config["splicing_rtqpcr"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/splicing_rtqpcr.R"

rule assemble_figure_splicing:
    input:
        fonts = ".fonts_registered.txt",
        splicing = "panels/splicing.Rdata",
        splicing_rtqpcr = "panels/splicing_rtqpcr.Rdata",
        rpgene_datavis = "panels/rpgene_datavis.Rdata",
    output:
        pdf = "figures/figure_7_splicing.pdf"
    params:
        fig_width = eval(str(config["splicing_figure"]["fig_width"])),
        fig_height = eval(str(config["splicing_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_splicing.R"


rule promoter_swap_diagram:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["promoter_swap_diagram"]["data"],
    output:
        pdf = "panels/promoter_swap_diagram.pdf",
        grob = "panels/promoter_swap_diagram.Rdata",
    params:
        fig_height = eval(str(config["promoter_swap_diagram"]["fig_height"])),
        fig_width = eval(str(config["promoter_swap_diagram"]["fig_width"])),
        panel_letter = config["promoter_swap_diagram"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/promoter_swap_diagram.R"

rule promoter_swap_rtqpcr:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["promoter_swap_rtqpcr"]["data"],
        rnaseq = config["promoter_swap_rtqpcr"]["rnaseq"]
    output:
        pdf = "panels/promoter_swap_rtqpcr.pdf",
        grob = "panels/promoter_swap_rtqpcr.Rdata",
    params:
        fig_height = eval(str(config["promoter_swap_rtqpcr"]["fig_height"])),
        fig_width = eval(str(config["promoter_swap_rtqpcr"]["fig_width"])),
        panel_letter = config["promoter_swap_rtqpcr"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/promoter_swap_rtqpcr.R"

rule assemble_figure_promoter_swap:
    input:
        fonts = ".fonts_registered.txt",
        promoter_swap_diagram = "panels/promoter_swap_diagram.Rdata",
        promoter_swap_rtqpcr = "panels/promoter_swap_rtqpcr.Rdata",
    output:
        pdf = "figures/figure_3_promoter_swap.pdf"
    params:
        fig_width = eval(str(config["promoter_swap_figure"]["fig_width"])),
        fig_height = eval(str(config["promoter_swap_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_promoter_swap.R"


rule coip_western:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        input_rpb3_blot = config["coip_western"]["input_rpb3_blot"],
        input_spt6_blot = config["coip_western"]["input_spt6_blot"],
        input_spn1_blot = config["coip_western"]["input_spn1_blot"],
        input_set2_blot = config["coip_western"]["input_set2_blot"],
        ip_rpb3_blot = config["coip_western"]["ip_rpb3_blot"],
        ip_spt6_blot = config["coip_western"]["ip_spt6_blot"],
        ip_spn1_blot = config["coip_western"]["ip_spn1_blot"],
        ip_set2_blot = config["coip_western"]["ip_set2_blot"],
    output:
        pdf = "panels/coip_western.pdf",
        grob = "panels/coip_western.Rdata",
    params:
        fig_height = eval(str(config["coip_western"]["fig_height"])),
        fig_width = eval(str(config["coip_western"]["fig_width"])),
        panel_letter = config["coip_western"]["panel_letter"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/coip_western.R"

rule set2_metagene:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["set2_metagene"]["data"],
    output:
        pdf = "panels/set2_metagene.pdf",
        grob = "panels/set2_metagene.Rdata",
    params:
        fig_height = eval(str(config["set2_metagene"]["fig_height"])),
        fig_width = eval(str(config["set2_metagene"]["fig_width"])),
        panel_letter = config["set2_metagene"]["panel_letter"],
        numerator_factor = "Set2",
        denominator_factor = "Rpb1",
        control_label_y = config["set2_metagene"]["control_label_y"],
        condition_label_y = config["set2_metagene"]["condition_label_y"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/single_factor_ratio_metagene.R"

rule chip_maplot:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = lambda wc: config["chip_maplot"][wc.factor]["data"],
        rpb1 = config["chip_maplot"]["rpb1"]
    output:
        pdf = "panels/{factor}_chip_maplot.pdf",
        grob = "panels/{factor}_chip_maplot.Rdata",
    params:
        factor_id = lambda wc: config["chip_maplot"][wc.factor]["factor_id"],
        rpb1_norm = lambda wc: config["chip_maplot"][wc.factor]["rpb1_norm"],
        fig_height = lambda wc: eval(str(config["chip_maplot"][wc.factor]["fig_height"])),
        fig_width = lambda wc: eval(str(config["chip_maplot"][wc.factor]["fig_width"])),
        panel_letter = lambda wc: config["chip_maplot"][wc.factor]["panel_letter"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/chip_maplot.R"

rule spt6_metagene:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["spt6_metagene"]["data"],
    output:
        pdf = "panels/spt6_metagene.pdf",
        grob = "panels/spt6_metagene.Rdata",
    params:
        fig_height = eval(str(config["spt6_metagene"]["fig_height"])),
        fig_width = eval(str(config["spt6_metagene"]["fig_width"])),
        panel_letter = config["spt6_metagene"]["panel_letter"],
        numerator_factor = "Spt6",
        denominator_factor = "Rpb1",
        control_label_y = config["spt6_metagene"]["control_label_y"],
        condition_label_y = config["spt6_metagene"]["condition_label_y"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/single_factor_ratio_metagene.R"

rule set2_abundance_chipseq_barplot:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["set2_abundance_chipseq_barplot"]["data"],
    output:
        pdf = "panels/set2_abundance_chipseq_barplot.pdf",
        grob = "panels/set2_abundance_chipseq_barplot.Rdata",
    params:
        fig_height = eval(str(config["set2_abundance_chipseq_barplot"]["fig_height"])),
        fig_width = eval(str(config["set2_abundance_chipseq_barplot"]["fig_width"])),
        panel_letter = config["set2_abundance_chipseq_barplot"]["panel_letter"],
        plot_title = "Set2 ChIP-seq"
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/chipseq_abundance_barplot.R"

rule spt6_abundance_chipseq_barplot:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["spt6_abundance_chipseq_barplot"]["data"],
    output:
        pdf = "panels/spt6_abundance_chipseq_barplot.pdf",
        grob = "panels/spt6_abundance_chipseq_barplot.Rdata",
    params:
        fig_height = eval(str(config["spt6_abundance_chipseq_barplot"]["fig_height"])),
        fig_width = eval(str(config["spt6_abundance_chipseq_barplot"]["fig_width"])),
        panel_letter = config["spt6_abundance_chipseq_barplot"]["panel_letter"],
        plot_title = "Spt6 ChIP-seq"
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/chipseq_abundance_barplot.R"

rule assemble_figure_set2_spt6:
    input:
        fonts = ".fonts_registered.txt",
        coip_western = "panels/coip_western.Rdata",
        set2_metagene = "panels/set2_metagene.Rdata",
        spt6_metagene = "panels/spt6_metagene.Rdata",
        set2_maplot = "panels/Set2_chip_maplot.Rdata",
        spt6_maplot = "panels/Spt6_chip_maplot.Rdata",
        set2_abundance_chipseq_barplot = "panels/set2_abundance_chipseq_barplot.Rdata",
        spt6_abundance_chipseq_barplot = "panels/spt6_abundance_chipseq_barplot.Rdata",
    output:
        pdf = "figures/figure_4_set2_spt6.pdf"
    params:
        fig_width = eval(str(config["set2_spt6_figure"]["fig_width"])),
        fig_height = eval(str(config["set2_spt6_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_set2_spt6.R"

rule set2_spt6_v_rpb1:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        spt6 = config["set2_spt6_v_rpb1"]["spt6"],
        set2 = config["set2_spt6_v_rpb1"]["set2"],
        rpb1 = config["set2_spt6_v_rpb1"]["rpb1"],
    output:
        pdf = "panels/set2_spt6_v_rpb1.pdf",
        grob = "panels/set2_spt6_v_rpb1.Rdata",
    params:
        fig_height = eval(str(config["set2_spt6_v_rpb1"]["fig_height"])),
        fig_width = eval(str(config["set2_spt6_v_rpb1"]["fig_width"])),
        panel_letter = config["set2_spt6_v_rpb1"]["panel_letter"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/set2_spt6_v_rpb1.R"

rule set2_spt6_rpb1norm_v_rpb1:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        spt6 = config["set2_spt6_rpb1norm_v_rpb1"]["spt6"],
        set2 = config["set2_spt6_rpb1norm_v_rpb1"]["set2"],
        rpb1 = config["set2_spt6_rpb1norm_v_rpb1"]["rpb1"],
    output:
        pdf = "panels/set2_spt6_rpb1norm_v_rpb1.pdf",
        grob = "panels/set2_spt6_rpb1norm_v_rpb1.Rdata",
    params:
        fig_height = eval(str(config["set2_spt6_rpb1norm_v_rpb1"]["fig_height"])),
        fig_width = eval(str(config["set2_spt6_rpb1norm_v_rpb1"]["fig_width"])),
        panel_letter = config["set2_spt6_rpb1norm_v_rpb1"]["panel_letter"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/set2_spt6_rpb1norm_v_rpb1.R"

rule assemble_figure_set2_spt6_supp:
    input:
        fonts = ".fonts_registered.txt",
        set2_spt6_v_rpb1 = "panels/set2_spt6_v_rpb1.Rdata",
        set2_spt6_rpb1norm_v_rpb1 = "panels/set2_spt6_rpb1norm_v_rpb1.Rdata",
    output:
        pdf = "figures/figure_S3_set2_spt6_supplemental.pdf"
    params:
        fig_width = eval(str(config["set2_spt6_supplemental"]["fig_width"])),
        fig_height = eval(str(config["set2_spt6_supplemental"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_set2_spt6_supp.R"


rule h3_metagene:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["h3_metagene"]["data"],
    output:
        pdf = "panels/h3_metagene.pdf",
        grob = "panels/h3_metagene.Rdata",
    params:
        fig_height = eval(str(config["h3_metagene"]["fig_height"])),
        fig_width = eval(str(config["h3_metagene"]["fig_width"])),
        panel_letter = config["h3_metagene"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_metagene.R"

rule h3_vs_rpb1_ma:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        h3_path = config["h3_vs_rpb1_ma"]["h3"],
        rpb1_path = config["h3_vs_rpb1_ma"]["rpb1"],
    output:
        pdf = "panels/h3_vs_rpb1_ma.pdf",
        grob = "panels/h3_vs_rpb1_ma.Rdata",
    params:
        fig_height = eval(str(config["h3_vs_rpb1_ma"]["fig_height"])),
        fig_width = eval(str(config["h3_vs_rpb1_ma"]["fig_width"])),
        panel_letter = config["h3_vs_rpb1_ma"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_vs_rpb1_ma_paper.R"

rule reduced_h3_h3_metagene:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["reduced_h3_h3_metagene"]["data"],
    output:
        pdf = "panels/reduced_h3_h3_metagene.pdf",
        grob = "panels/reduced_h3_h3_metagene.Rdata",
    params:
        fig_height = eval(str(config["reduced_h3_h3_metagene"]["fig_height"])),
        fig_width = eval(str(config["reduced_h3_h3_metagene"]["fig_width"])),
        panel_letter = config["reduced_h3_h3_metagene"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/reduced_h3_h3_metagene.R"

rule h3_single_locus_datavis:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data_paths = list(config["h3_single_locus_datavis"]["data"].values()),
        transcript_annotations = config["h3_single_locus_datavis"]["transcript_annotation"],
        orf_annotations = config["h3_single_locus_datavis"]["orf_annotation"]
    output:
        pdf = "panels/h3_single_locus_datavis.pdf",
        grob = "panels/h3_single_locus_datavis.Rdata",
    params:
        targets = list(config["h3_single_locus_datavis"]["data"].keys()),
        fig_height = eval(str(config["h3_single_locus_datavis"]["fig_height"])),
        fig_width = eval(str(config["h3_single_locus_datavis"]["fig_width"])),
        panel_letter = config["h3_single_locus_datavis"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_single_locus_datavis_paper.R"

rule assemble_figure_h3:
    input:
        fonts = ".fonts_registered.txt",
        h3_vs_rpb1_ma = "panels/h3_vs_rpb1_ma.Rdata",
        reduced_h3_h3_metagene = "panels/reduced_h3_h3_metagene.Rdata",
        h3_single_locus_datavis = "panels/h3_single_locus_datavis.Rdata",
    output:
        pdf = "figures/figure_5_h3.pdf"
    params:
        fig_width = eval(str(config["h3_figure"]["fig_width"])),
        fig_height = eval(str(config["h3_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_h3.R"

rule h3_vs_rpb1_ma_5p500bp:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        h3_path = config["h3_vs_rpb1_ma_5p500bp"]["h3"],
        rpb1_path = config["h3_vs_rpb1_ma_5p500bp"]["rpb1"],
    output:
        pdf = "panels/h3_vs_rpb1_ma_5p500bp.pdf",
        grob = "panels/h3_vs_rpb1_ma_5p500bp.Rdata",
    params:
        fig_height = eval(str(config["h3_vs_rpb1_ma_5p500bp"]["fig_height"])),
        fig_width = eval(str(config["h3_vs_rpb1_ma_5p500bp"]["fig_width"])),
        panel_letter = config["h3_vs_rpb1_ma_5p500bp"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_vs_rpb1_ma_5p500bp.R"

rule h3_plmin_reduced_lfc_scatterplots:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data_paths = [v for k,v in config["h3_plmin_reduced_lfc_scatterplots"]["data"].items()],
        h3_path = config["h3_plmin_reduced_lfc_scatterplots"]["h3"]
    output:
        pdf = "panels/h3_plmin_reduced_lfc_scatterplots.pdf",
        grob = "panels/h3_plmin_reduced_lfc_scatterplots.Rdata",
    params:
        chip_factors = [k for k,v in config["h3_plmin_reduced_lfc_scatterplots"]["data"].items()],
        fig_height = eval(str(config["h3_plmin_reduced_lfc_scatterplots"]["fig_height"])),
        fig_width = eval(str(config["h3_plmin_reduced_lfc_scatterplots"]["fig_width"])),
        panel_letter = config["h3_plmin_reduced_lfc_scatterplots"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_plmin_reduced_lfc_scatterplots.R"

rule reduced_h3_matched_metagenes:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["reduced_h3_matched_metagenes"]["data"],
    output:
        pdf = "panels/reduced_h3_matched_metagenes.pdf",
        grob = "panels/reduced_h3_matched_metagenes.Rdata",
    params:
        fig_height = eval(str(config["reduced_h3_matched_metagenes"]["fig_height"])),
        fig_width = eval(str(config["reduced_h3_matched_metagenes"]["fig_width"])),
        panel_letter = config["reduced_h3_matched_metagenes"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/reduced_h3_matched_metagenes.R"

rule assemble_figure_h3_supp:
    input:
        fonts = ".fonts_registered.txt",
        h3_metagene = "panels/h3_metagene.Rdata",
        h3_vs_rpb1_ma_5p500bp = "panels/h3_vs_rpb1_ma_5p500bp.Rdata",
        h3_plmin_reduced_lfc_scatterplots = "panels/h3_plmin_reduced_lfc_scatterplots.Rdata",
        reduced_h3_matched_metagenes = "panels/reduced_h3_matched_metagenes.Rdata",
    output:
        pdf = "figures/figure_S4_h3_supplemental.pdf"
    params:
        fig_width = eval(str(config["h3_supplemental"]["fig_width"])),
        fig_height = eval(str(config["h3_supplemental"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_h3_supp.R"

rule h3_modification_datavis_ratio:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        metagene_data = config["h3_modification_datavis"]["metagene_data"],
        heatmap_data = config["h3_modification_datavis"]["heatmap_data"],
        annotation = config["h3_modification_datavis"]["annotation"],
    output:
        quantification =  "panels/h3_modification_datavis_ratio_quantification.tsv",
        pdf = "panels/h3_modification_datavis_ratio.pdf",
        grob = "panels/h3_modification_datavis_ratio.Rdata",
    params:
        fig_height = eval(str(config["h3_modification_datavis_ratio"]["fig_height"])),
        fig_width = eval(str(config["h3_modification_datavis_ratio"]["fig_width"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_modification_datavis_ratio.R"

rule h3_mods_single_locus_datavis:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data_paths = list(config["h3_mods_single_locus_datavis"]["data"].values()),
        transcript_annotations = config["h3_mods_single_locus_datavis"]["transcript_annotation"],
        orf_annotations = config["h3_mods_single_locus_datavis"]["orf_annotation"]
    output:
        pdf = "panels/h3_mods_single_locus_datavis.pdf",
        grob = "panels/h3_mods_single_locus_datavis.Rdata",
    params:
        targets = list(config["h3_mods_single_locus_datavis"]["data"].keys()),
        fig_height = eval(str(config["h3_mods_single_locus_datavis"]["fig_height"])),
        fig_width = eval(str(config["h3_mods_single_locus_datavis"]["fig_width"])),
        panel_letter = config["h3_mods_single_locus_datavis"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_mods_single_locus_datavis.R"

rule assemble_figure_h3_mods:
    input:
        fonts = ".fonts_registered.txt",
        h3_modification_datavis = "panels/h3_modification_datavis_ratio.Rdata",
    output:
        pdf = "figures/figure_6_h3_mods.pdf"
    params:
        fig_width = eval(str(config["h3_mods_figure"]["fig_width"])),
        fig_height = eval(str(config["h3_mods_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_h3_mods.R"

rule chipseq_abundance_barplots_h3:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = list(config["chipseq_abundance_barplots_h3"]["data"].values()),
    output:
        pdf = "panels/chipseq_abundance_barplots_h3.pdf",
        grob = "panels/chipseq_abundance_barplots_h3.Rdata",
    params:
        factor_ids = list(config["chipseq_abundance_barplots_h3"]["data"].keys()),
        fig_height = eval(str(config["chipseq_abundance_barplots_h3"]["fig_height"])),
        fig_width = eval(str(config["chipseq_abundance_barplots_h3"]["fig_width"])),
        panel_letter = config["chipseq_abundance_barplots_h3"]["panel_letter"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/chipseq_abundance_barplots_multifactor.R"

rule h3_mods_non_h3_norm:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["h3_mods_non_h3_norm"]["data"],
    output:
        quantification =  "panels/h3_mods_non_h3_norm_quantification.tsv",
        pdf = "panels/h3_mods_non_h3_norm.pdf",
        grob = "panels/h3_mods_non_h3_norm.Rdata",
    params:
        fig_height = eval(str(config["h3_mods_non_h3_norm"]["fig_height"])),
        fig_width = eval(str(config["h3_mods_non_h3_norm"]["fig_width"])),
        panel_letter = config["h3_mods_non_h3_norm"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_mods_non_h3_norm.R"

rule h3_mods_facet_expression_length:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = lambda wc: config["h3_mods_facet_expression_length"][wc.mod]["data"]
    output:
        pdf = "panels/{mod}_facet_expression_length.pdf",
        grob = "panels/{mod}_facet_expression_length.Rdata",
    params:
        fig_height = eval(str(config["h3_mods_facet_expression_length"]["fig_height"])),
        fig_width = eval(str(config["h3_mods_facet_expression_length"]["fig_width"])),
        panel_letter = lambda wc: config["h3_mods_facet_expression_length"][wc.mod]["panel_letter"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_mods_facet_expression_length.R"

rule h3_mods_facet_expression:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = lambda wc: config["h3_mods_facet_expression"]["data"]
    output:
        pdf = "panels/h3_mods_facet_expression.pdf",
        grob = "panels/h3_mods_facet_expression.Rdata",
    params:
        fig_height = eval(str(config["h3_mods_facet_expression"]["fig_height"])),
        fig_width = eval(str(config["h3_mods_facet_expression"]["fig_width"])),
        panel_letter = config["h3_mods_facet_expression"]["panel_letter"],
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_mods_facet_expression.R"

rule assemble_figure_h3_mods_supp:
    input:
        fonts = ".fonts_registered.txt",
        chipseq_abundance_barplots_h3 = "panels/chipseq_abundance_barplots_h3.Rdata",
        h3_mods_non_h3_norm = "panels/h3_mods_non_h3_norm.Rdata",
        h3_mods_facet_expression = "panels/h3_mods_facet_expression.Rdata",
        # h3k4me3 = "panels/H3K4me3_facet_expression_length.Rdata",
        # h3k36me2 = "panels/H3K36me2_facet_expression_length.Rdata",
        # h3k36me3 = "panels/H3K36me3_facet_expression_length.Rdata",
    output:
        pdf = "figures/figure_S5_h3_mods_supplemental.pdf"
    params:
        fig_width = eval(str(config["h3_mods_supplemental"]["fig_width"])),
        fig_height = eval(str(config["h3_mods_supplemental"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_h3_mods_supp.R"

rule rpgene_datavis:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["rpgene_datavis"]["data"],
    output:
        pdf = "panels/rpgene_datavis.pdf",
        # pdf = "panels/rpgene_datavis_length_filtered.pdf",
        grob = "panels/rpgene_datavis.Rdata",
        # grob = "panels/rpgene_datavis_length_filtered.Rdata",
    params:
        fig_height = eval(str(config["rpgene_datavis"]["fig_height"])),
        fig_width = eval(str(config["rpgene_datavis"]["fig_width"])),
        panel_letter = config["rpgene_datavis"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rpgene_datavis.R"

rule rpgene_datavis_supp:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["rpgene_datavis_supp"]["data"],
    output:
        pdf = "panels/rpgene_datavis_supp.pdf",
        # pdf = "panels/rpgene_datavis_supp_length_filtered.pdf",
        grob = "panels/rpgene_datavis_supp.Rdata",
        # grob = "panels/rpgene_datavis_supp_length_filtered.Rdata",
    params:
        fig_height = eval(str(config["rpgene_datavis_supp"]["fig_height"])),
        fig_width = eval(str(config["rpgene_datavis_supp"]["fig_width"])),
        panel_letter = config["rpgene_datavis_supp"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/rpgene_datavis.R"

rule assemble_figure_splicing_supp:
    input:
        fonts = ".fonts_registered.txt",
        rpgene_datavis_supp = "panels/rpgene_datavis_supp.Rdata",
    output:
        pdf = "figures/figure_S6_splicing_supplemental.pdf"
    params:
        fig_width = eval(str(config["splicing_supplemental"]["fig_width"])),
        fig_height = eval(str(config["splicing_supplemental"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_splicing_supp.R"

rule intron_containing_gene_heatmaps:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["intron_containing_gene_heatmaps"]["data"],
        rp_with_intron_bed = config["intron_containing_gene_heatmaps"]["rp_with_intron_bed"],
        nonrp_with_intron_bed = config["intron_containing_gene_heatmaps"]["nonrp_with_intron_bed"],
        intron_bed = config["intron_containing_gene_heatmaps"]["intron_bed"],
        name_lookup = config["intron_containing_gene_heatmaps"]["name_lookup"],
    output:
        pdf = "panels/intron_containing_gene_heatmaps.pdf",
        grob = "panels/intron_containing_gene_heatmaps.Rdata",
    params:
        fig_height = eval(str(config["intron_containing_gene_heatmaps"]["fig_height"])),
        fig_width = eval(str(config["intron_containing_gene_heatmaps"]["fig_width"])),
        panel_letter = config["intron_containing_gene_heatmaps"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/intron_containing_gene_heatmaps.R"

