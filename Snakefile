#!/usr/bin/env python

configfile: "config.yaml"

rule target:
    input:
        "panels/spn1_depletion_western.pdf",
        "panels/spn1_depletion_scatter.pdf",
        "panels/spn1_depletion_metagene.pdf",
        "panels/spn1_depletion_chipseq_barplot.pdf",
        "figures/figure_spn1_depletion.pdf",
        "panels/spn1_depletion_viability.pdf",
        "panels/spn1_v_rpb1_nondepleted.pdf",
        "panels/spn1_rpb1norm_v_rpb1_nondepleted.pdf",
        "figures/figure_spn1_depletion_supplemental.pdf",
        "panels/rnaseq_maplot.pdf",
        "panels/rnaseq_single_locus_datavis.pdf",
        "panels/rnaseq_vs_rpb1_single_locus.pdf",
        "panels/rpb1_metagenes.pdf",
        "panels/rnaseq_vs_rpb1.pdf",
        "figures/figure_rnaseq_rpb1.pdf",
        "panels/rnaseq_single_locus_datavis_supp.pdf",
        "panels/differential_expression_rtqpcr.pdf",
        "panels/antisense_single_locus_datavis.pdf",
        "panels/chipseq_abundance_barplots_rpb1.pdf",
        "panels/slow_growth_signature.pdf",
        "figures/figure_rnaseq_rpb1_supplemental.pdf",
        "panels/splicing.pdf",
        "panels/splicing_rtqpcr.pdf",
        "figures/figure_splicing.pdf",
        "panels/promoter_swap_diagram.pdf",
        "panels/promoter_swap_rtqpcr.pdf",
        "figures/figure_promoter_swap.pdf",
        "panels/coip_western.pdf",
        "panels/set2_metagene.pdf",
        "panels/spt6_metagene.pdf",
        "panels/set2_abundance_chipseq_barplot.pdf",
        "panels/spt6_abundance_chipseq_barplot.pdf",
        "figures/figure_set2_spt6.pdf",
        "panels/set2_spt6_v_rpb1.pdf",
        "panels/h3_vs_rpb1_ma.pdf",
        "panels/reduced_h3_h3_metagene.pdf",
        "panels/h3_single_locus_datavis.pdf",
        "figures/figure_h3.pdf",
        # "panels/h3_metagene.pdf",
        "panels/h3_modification_datavis.pdf",
        "figures/figure_h3_mods.pdf"

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
        pdf = "figures/figure_spn1_depletion.pdf"
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

rule spn1_v_rpb1_nondepleted:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        spn1 = config["spn1_v_rpb1_nondepleted"]["spn1"],
        rpb1 = config["spn1_v_rpb1_nondepleted"]["rpb1"],
    output:
        pdf = "panels/spn1_v_rpb1_nondepleted.pdf",
        grob = "panels/spn1_v_rpb1_nondepleted.Rdata",
    params:
        fig_height = eval(str(config["spn1_v_rpb1_nondepleted"]["fig_height"])),
        fig_width = eval(str(config["spn1_v_rpb1_nondepleted"]["fig_width"])),
        panel_letter = config["spn1_v_rpb1_nondepleted"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_v_rpb1_nondepleted.R"

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
        spn1_v_rpb1_nondepleted = "panels/spn1_v_rpb1_nondepleted.Rdata",
        spn1_rpb1norm_v_rpb1_nondepleted = "panels/spn1_rpb1norm_v_rpb1_nondepleted.Rdata",
    output:
        pdf = "figures/figure_spn1_depletion_supplemental.pdf"
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
        pdf = "figures/figure_rnaseq_rpb1.pdf"
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
        "scripts/chipseq_abundance_barplots_rpb1.R"

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
        antisense_single_locus_datavis = "panels/antisense_single_locus_datavis.Rdata",
        chipseq_abundance_barplots_rpb1 = "panels/chipseq_abundance_barplots_rpb1.Rdata",
    output:
        pdf = "figures/figure_rnaseq_rpb1_supplemental.pdf"
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
    output:
        pdf = "figures/figure_splicing.pdf"
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
        pdf = "figures/figure_promoter_swap.pdf"
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
        plot_title = "Set2 ChIP-seq",
        plot_subtitle = "3087 non-overlapping coding genes",
        legend_position = [0.6, 0.3]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/single_factor_standardized_metagene.R"

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
        plot_title = "Spt6 ChIP-seq",
        plot_subtitle = "3087 non-overlapping coding genes",
        legend_position = [0.6, 0.3]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/single_factor_standardized_metagene.R"

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
        set2_abundance_chipseq_barplot = "panels/set2_abundance_chipseq_barplot.Rdata",
        spt6_abundance_chipseq_barplot = "panels/spt6_abundance_chipseq_barplot.Rdata",
    output:
        pdf = "figures/figure_set2_spt6.pdf"
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
        pdf = "figures/figure_h3.pdf"
    params:
        fig_width = eval(str(config["h3_figure"]["fig_width"])),
        fig_height = eval(str(config["h3_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_h3.R"

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

rule h3_modification_datavis:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["h3_modification_datavis"]["data"],
        annotation = config["h3_modification_datavis"]["annotation"],
    output:
        pdf = "panels/h3_modification_datavis.pdf",
        grob = "panels/h3_modification_datavis.Rdata",
    params:
        fig_height = eval(str(config["h3_modification_datavis"]["fig_height"])),
        fig_width = eval(str(config["h3_modification_datavis"]["fig_width"])),
        panel_letter = config["h3_modification_datavis"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/h3_modification_datavis.R"

rule assemble_figure_h3_mods:
    input:
        fonts = ".fonts_registered.txt",
        # h3_metagene = "panels/h3_metagene.Rdata",
        h3_modification_datavis = "panels/h3_modification_datavis.Rdata",
    output:
        pdf = "figures/figure_h3_mods.pdf"
    params:
        fig_width = eval(str(config["h3_mods_figure"]["fig_width"])),
        fig_height = eval(str(config["h3_mods_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_h3_mods.R"

