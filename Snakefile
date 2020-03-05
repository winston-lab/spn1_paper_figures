#!/usr/bin/env python

configfile: "config.yaml"

rule target:
    input:
        "panels/spn1_depletion_western.pdf",
        "panels/spn1_depletion_metagene.pdf",
        "panels/spn1_depletion_chipseq_barplot.pdf",
        "figures/figure_spn1_depletion.pdf",
        "panels/spn1_depletion_viability.pdf",
        "panels/rnaseq_maplot.pdf",
        "panels/rpb1_metagenes.pdf",
        "panels/rnaseq_vs_rpb1.pdf",
        "figures/figure_rnaseq_rpb1.pdf",
        "panels/histone_metagenes.pdf",

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
        # quant_data = config["spn1_depletion_western"]["quant_data"],
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
        panel_letter = config["spn1_depletion_chipseq_barplot"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_depletion_chipseq_barplot.R"

rule assemble_figure_spn1_depletion:
    input:
        fonts = ".fonts_registered.txt",
        spn1_depletion_western = "panels/spn1_depletion_western.Rdata",
        spn1_depletion_chipseq_barplot = "panels/spn1_depletion_chipseq_barplot.Rdata",
        spn1_depletion_metagene = "panels/spn1_depletion_metagene.Rdata",
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

rule rpb1_metagenes:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["rpb1_metagenes"]["data"],
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

rule histone_metagenes:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["histone_metagenes"]["data"],
    output:
        pdf = "panels/histone_metagenes.pdf",
        grob = "panels/histone_metagenes.Rdata",
    params:
        fig_height = eval(str(config["histone_metagenes"]["fig_height"])),
        fig_width = eval(str(config["histone_metagenes"]["fig_width"])),
        panel_letter = config["histone_metagenes"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/histone_metagenes.R"

