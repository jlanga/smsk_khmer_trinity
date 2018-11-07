rule mapping_index:
    input:
        fasta = ASSEMBLY + "Trinity.fasta"
    output:
        indexes = expand(
            ASSEMBLY + "Trinity.{ending}.bt2",
            ending="1 2 3 4 rev.1 rev.2".split(" ")
        )
    params:
        index_prefix = ASSEMBLY + "Trinity"
    log:
        ASSEMBLY + "bowtie2_index.log"
    benchmark:
        ASSEMBLY + "bowtie2_index.bmk"
    conda:
        "mapping.yml"
    shell:
        "bowtie2-build  {input.fasta} {params.index_prefix} "
        "2> {log} 1>&2"


rule mapping_bowtie2:
    """
    Map reads with bowtie
    """
    input:
        fasta = ASSEMBLY + "Trinity.fasta",
        pe_pe = QC + "{sample}_pepe.fq.gz",
        pe_se = QC + "{sample}_pese.fq.gz",
        indexes = rules.mapping_index.output.indexes
    output:
        cram = MAPPING + "{sample}.cram"
    log:
        MAPPING + "{sample}.log"
    benchmark:
        MAPPING + "{sample}.bmk"
    params:
        sample = "{sample}",
        index_prefix = ASSEMBLY + "Trinity"
    conda:
        "mapping.yml"
    threads:
        4
    shell:
        "(bowtie2"
        "   --threads {threads}"
        "   --no-unal"
        "   --local"
        "   --rg {params.sample}"
        "   -x {params.index_prefix}"
        "   --interleaved {input.pe_pe}"
        "   -U {input.pe_se}"
        "| samtools sort"
        "   -l 9"
        "   -@ {threads}"
        "   --reference {input.fasta}"
        "   -o {output.cram}"
        "   --output-fmt CRAM"
        "   -)"
        "2> {log} 1>&2"


rule mapping:
    input:
        expand(
            MAPPING + "{sample}.cram",
            sample=SAMPLES_PE
        )
