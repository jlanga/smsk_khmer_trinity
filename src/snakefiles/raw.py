def get_path_forward(wildcards):
    """"""
    return samples[wildcards.sample]["forward"]


def get_path_reverse(wildcards):
    """"""
    return samples[wildcards.sample]["reverse"]


rule raw_make_links_pe:
    """
    Make a link next to the original file, with a prettier name than default.
    """
    input:
        forward = get_path_forward,
        reverse = get_path_reverse
    output:
        forward = protected(RAW + "{sample}_1.fq.gz"),
        reverse = protected(RAW + "{sample}_2.fq.gz")
    log:
        RAW + "make_links_pe_{sample}.log"
    benchmark:
        RAW + "make_links_pe_{sample}.bmk"
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
    return samples[wildcards.sample]["single"]


rule raw_make_links_se:
    """
    Make a link next to the original file, with a prettier name than default.
    """
    input:
        single = get_path_single,
    output:
        single = protected(RAW + "{sample}_se.fq.gz")
    log:
        RAW + "make_links_se_{sample}.log"
    benchmark:
        RAW + "make_links_se_{sample}.bmk"
    shell:
        """
        ln --symbolic \
            $(readlink --canonicalize {input.single}) \
            {output.single} \
        2> {log}
        """


rule raw_results:
    """Checkpoint to generate all the links for raw data"""
    input:
        expand(
            RAW + "{sample}_{end}.fq.gz",
            sample=SAMPLES_PE,
            end="1 2".split()
        ),
        expand(
            RAW + "{sample}_se.fq.gz",
            sample=SAMPLES_SE
        )


rule raw_doc:
    """Checkpoint to generate all reports for raw data"""
    input:
        html = RAW + "multiqc_report.html"


rule raw:
    """Make both results + reports for raw (raw = raw_doc)"""
    input:
        rules.raw_results.input,
        rules.raw_doc.input
