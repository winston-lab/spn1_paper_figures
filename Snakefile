#!/usr/bin/env python

configfile: "config.yaml"

rule target:
    input:
        "figures/spn1_depletion_metagene.pdf",
        "figures/figure_spn1_depletion.pdf",
        "figures/spn1_depletion_chipseq_barplot.pdf",

rule register_fonts:
    input:
        fonts_path = config["fonts_path"],
    output:
        output_path = ".fonts_registered.txt"
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/register_fonts.R"

rule spn1_depletion_metagene:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["spn1_depletion_metagene"]["data"],
    output:
        pdf = "figures/spn1_depletion_metagene.pdf",
        grob = "figures/spn1_depletion_metagene.Rdata",
    params:
        fig_height = eval(str(config["spn1_depletion_metagene"]["fig_height"])),
        fig_width = eval(str(config["spn1_depletion_metagene"]["fig_width"])),
        panel_letter = config["spn1_depletion_metagene"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_depletion_metagene.R"

rule assemble_figure_spn1_depletion:
    input:
        fonts = ".fonts_registered.txt",
        spn1_depletion_metagene = "figures/spn1_depletion_metagene.Rdata",
    output:
        pdf = "figures/figure_spn1_depletion.pdf"
    params:
        fig_width = eval(str(config["spn1_depletion_figure"]["fig_width"])),
        fig_height = eval(str(config["spn1_depletion_figure"]["fig_height"])),
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/assemble_figure_spn1_depletion.R"

rule spn1_depletion_chipseq_barplot:
    input:
        fonts = ".fonts_registered.txt",
        theme = config["theme_path"],
        data = config["spn1_depletion_chipseq_barplot"]["data"],
    output:
        pdf = "figures/spn1_depletion_chipseq_barplot.pdf",
        grob = "figures/spn1_depletion_chipseq_barplot.Rdata",
    params:
        fig_height = eval(str(config["spn1_depletion_chipseq_barplot"]["fig_height"])),
        fig_width = eval(str(config["spn1_depletion_chipseq_barplot"]["fig_width"])),
        panel_letter = config["spn1_depletion_chipseq_barplot"]["panel_letter"]
    conda:
        "envs/plot_figures.yaml"
    script:
        "scripts/spn1_depletion_chipseq_barplot.R"




