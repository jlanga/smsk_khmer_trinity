def get_path_forward(wildcards):
    """"""
    return config["samples_pe"][wildcards.sample]["forward"]


def get_path_reverse(wildcards):
    """"""
    return config["samples_pe"][wildcards.sample]["reverse"]



rule raw_make_links_pe_sample:
    """
    Make a link next to the original file, with a prettier name than default.
    """
    input:
        forward = get_path_forward,
        reverse = get_path_reverse
    output:
        forward = protected(raw + "{sample}_1.fq.gz"),
        reverse = protected(raw + "{sample}_2.fq.gz")
    log:
        raw + "make_links_pe_{sample}.log"
    benchmark:
        raw + "make_links_pe_{sample}.bmk"
    shell:
        """
        (ln --symbolic \
            $(readlink --canonicalize {input.forward}) \
            {output.forward}
        ln --symbolic \
            $(readlink --canonicalize {input.reverse}) \
            {output.reverse} ) \
        2>> {log}
        """


def get_path_single(wildcards):
    """"""
    return config["samples_se"][wildcards.sample]["single"]


rule raw_make_links_se_sample:
    """
    Make a link next to the original file, with a prettier name than default.
    """
    input:
        single = get_path_single,
    output:
        single = protected(raw + "{sample}_se.fq.gz")
    log:
        raw + "make_links_se_{sample}.log"
    benchmark:
        raw + "make_links_se_{sample}.bmk"
    shell:
        """
        ln --symbolic \
            $(readlink --canonicalize {input.single}) \
            {output.single} \
        2> {log}
        """


rule raw_fastqc_pe_sample:
    """FASTQC over PE files"""
    input:
        forward = raw + "{sample}_1.fq.gz",
        reverse = raw + "{sample}_2.fq.gz"
    output:
        html1 = protected(raw + "{sample}_1_fastqc.html"),
        html2 = protected(raw + "{sample}_2_fastqc.html"),
        zip1 = protected(raw + "{sample}_1_fastqc.zip"),
        zip2 = protected(raw + "{sample}_2_fastqc.zip")

    threads:
        2
    params:
        outdir = raw
    log:
        raw + "fastqc_pe_{sample}.log"
    benchmark:
        raw + "fastqc_pe_{sample}.bmk"
    conda:
        "raw.yml"
    shell:
        """
        fastqc \
            --threads {threads} \
            --nogroup \
            --outdir {params.outdir} \
            {input.forward} \
            {input.reverse} \
        2> {log} 1>&2
        """


rule raw_fastqc_se_sample:
    """FASTQC over SE files"""
    input:
        fq = raw + "{sample}_se.fq.gz",
    output:
        html = protected(raw + "{sample}_se_fastqc.html"),
        zip = protected(raw + "{sample}_se_fastqc.zip")
    params:
        outdir = raw
    log:
        raw + "fastqc_se_{sample}.log"
    benchmark:
        raw + "fastqc_se_{sample}.bmk"
    conda:
        "raw.yml"
    shell:
        """
        fastqc \
            --nogroup \
            --outdir {params.outdir} \
            {input.fq} \
        2> {log} 1>&2
        """


rule raw_multiqc:
    '''MultiQC report over the FASTQC ones'''
    input:
        pe_files = expand(
            raw + "{sample}_{pair}_fastqc.{extension}",
            sample=SAMPLES_PE,
            pair="1 2".split(),
            extension="html zip".split()
        ),
        se_files = expand(
            raw + "{sample}_se_fastqc.{extension}",
            sample=SAMPLES_SE,
            extension="html zip".split()
        )
    output:
        html = protected(raw + "multiqc_report.html")
    params:
        folder = raw
    log:
        raw + "multiqc.log"
    benchmark:
        raw + "multiqc.bmk"
    conda:
        "raw.yml"
    shell:
        """
        multiqc \
            --title Raw \
            --filename {output.html} \
            {params.folder} \
        2> {log}
        """


rule raw_results:
    """Checkpoint to generate all the links for raw data"""
    input:
        expand(
            raw + "{sample}_{end}.fq.gz",

            sample=SAMPLES_PE,
            end="1 2".split()
        ),
        expand(
            raw + "{sample}_se.fq.gz",
            sample=SAMPLES_SE
        )

        
rule raw_doc:
    """Checkpoint to generate all reports for raw data"""
    input:
        html = raw + "multiqc_report.html"


rule raw:
    """Make both results + reports for raw (raw = raw_doc)"""
    input:
        pe_files = expand(
            raw + "{sample}_{end}.fq.gz",
            sample=SAMPLES_PE,
            end="1 2".split()
        ),
        se_files = expand(
            raw + "{sample}_se.fq.gz",
            sample=SAMPLES_SE
        ),
        html = raw + "multiqc_report.html"
