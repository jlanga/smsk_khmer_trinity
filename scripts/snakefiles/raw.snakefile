rule raw_make_links_pe_sample:
    """Make a link next to the original file, with a prettier name than default.
    """
    input:
        forward= lambda wildcards: config["samples_pe"][wildcards.sample]["forward"],
        reverse= lambda wildcards: config["samples_pe"][wildcards.sample]["reverse"]
    output:
        forward= RAW_DIR + "{sample}_1.fq.gz",
        reverse= RAW_DIR + "{sample}_2.fq.gz"
    log:
        RAW_DOC + "make_links_pe_{sample}.log"
    benchmark:
        RAW_DOC + "make_links_pe_{sample}.json"
    shell:
        "ln -s $(readlink -f {input.forward}) {output.forward} 2> {log};"
        "ln -s $(readlink -f {input.reverse}) {output.reverse} 2>> {log}"



rule raw_make_links_se_sample:
    """Make a link next to the original file, with a prettier name than default.
    """
    input:
        single= lambda wildcards: config["samples_se"][wildcards.sample]["single"],
    output:
        single= RAW_DIR + "{sample}_se.fq.gz"
    log:
        RAW_DOC + "make_links_se_{sample}.log"
    benchmark:
        RAW_DOC + "make_links_se_{sample}.json"
    shell:
        """
        ln -s $(readlink -f {input.single}) {output.single} 2>  {log}
        """


rule raw_fastqc_pe_sample:
    """FASTQC over PE files"""
    input:
        forward = RAW_DIR + "{sample}_1.fq.gz",
        reverse = RAW_DIR + "{sample}_2.fq.gz"
    output:
        html1 = RAW_DOC + "{sample}_1_fastqc.html",
        html2 = RAW_DOC + "{sample}_2_fastqc.html",
        zip1=  RAW_DOC + "{sample}_1_fastqc.zip",
        zip2=  RAW_DOC + "{sample}_2_fastqc.zip"
    threads:
        2
    params:
        outdir = RAW_DOC
    log:
        RAW_DOC + "fastqc_pe_{sample}.log"
    benchmark:
        RAW_DOC + "fastqc_pe_{sample}.json"
    shell:
        "fastqc "
            "--nogroup "
            "--outdir {params.outdir} "
            "{input.forward} {input.reverse} "
            "2> {log} 1>&2"



rule raw_fastqc_se_sample:
    """FASTQC over SE files"""
    input:
        fq = RAW_DIR + "{sample}_se.fq.gz",
    output:
        html= RAW_DOC + "{sample}_se_fastqc.html",
        zip=  RAW_DOC + "{sample}_se_fastqc.zip"
    params:
        outdir = RAW_DOC
    log:
        RAW_DOC + "fastqc_se_{sample}.log"
    benchmark:
        RAW_DOC + "fastqc_se_{sample}.json"
    shell:
        "fastqc "
            "--nogroup "
            "--outdir {params.outdir} "
            "{input.fq} "
            "2> {log} 1>&2"



rule raw_results:
    """Checkpoint to generate all the links for raw data"""
    input:
        expand(
            RAW_DIR + "{sample}_{end}.fq.gz",
            sample = SAMPLES_PE,
            end = "1 2".split()
        ),
        expand(
            RAW_DIR + "{sample}_se.fq.gz",
            sample = SAMPLES_SE
        )



rule raw_doc:
    """Checkpoint to generate all reports for raw data"""
    input:
        expand(
            RAW_DOC + "{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "1 2".split(),
            extension = ["html", "zip"]
        ),
        expand(
            RAW_DOC + "{sample}_se_fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = ["html", "zip"]
        )



rule raw:
    """Make both results + reports for raw (raw = raw_doc)"""
    input:
        expand(
            RAW_DIR + "{sample}_{end}.fq.gz",
            sample = SAMPLES_PE,
            end = "1 2".split()
        ),
        expand(
            RAW_DIR + "{sample}_se.fq.gz",
            sample = SAMPLES_SE
        ),
        expand(
            RAW_DOC + "{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "1 2".split(),
            extension = ["html", "zip"]
        ),
        expand(
            RAW_DOC + "{sample}_se_fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = ["html", "zip"]
        )
