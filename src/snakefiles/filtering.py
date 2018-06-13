rule filtering_index_trinity:
    """
    Make a kallisto index for the assembly
    """
    input:
        assembly = assembly + "Trinity.fasta"
    output:
        index    = filtering + "Trinity.idx"
    log:
        filtering + "index_trinity.log"
    benchmark:
        filtering + "index_trinity.json"
    shell:
        "kallisto index "
            "--index={output.index} "
            "{input.assembly} "
        "2> {log}"





rule filtering_pseudoalign_all:
    """
    Pseudomap and quantify all samples together.
    """
    input:
        forward = expand(
            qc + "{sample}_1.fq",
            sample = SAMPLES_PE
        ),
        reverse = expand(
            qc + "{sample}_2.fq",
            sample = SAMPLES_PE
        ),
        index   = filtering   + "Trinity.idx"
    output:
        abundance_tsv = filtering + "raw/abundance.tsv",
        abundance_h5  = filtering + "raw/abundance.h5",
    params:
        outdir  = filtering + "raw/",
        params  = config["kallisto_params"]
    threads:
        24
    log:
        filtering + "raw/filtering_all.log"
    benchmark:
        filtering + "raw/filtering_all.json"
    shell:
        "kallisto quant "
            "--index={input.index} "
            "--output={params.outdir} "
            "--threads={threads} "
            "{params.params} "
            "<(cat {input.forward}) "
            "<(cat {input.reverse}) "
        "2> {log}"



rule filtering_filtering_transcriptome_id:
    """
    Filter the transcripts as at least 5% per loci/gene
    """
    input:
        abundance_tsv = filtering + "raw/abundance.tsv",
        assembly_map = assembly + "Trinity_gene_to_trans.tsv"
    output:
        tsv = filtering + "filtered_transcriptome_id.tsv"
    log:
        filtering + "filtering_transcriptome_id.log"
    benchmark:
        filtering + "filtering_transcriptome_id.json"
    shell:
        "Rscript src/get_expressed_ids.R "
            "--quantification {input.abundance_tsv} "
            "--gene_to_transcript_map {input.assembly_map} "
            "--expressed_ids {output.tsv} "
            "--tpms_percent 0.05 "
        "2> {log}"



rule filtering_filtering_transcriptome_fasta:
    """
    Get the filtered fasta
    """
    input:
        assembly = assembly + "Trinity.fasta",
        fai = assembly + "Trinity.fasta.fai",
        ids = filtering + "filtered_transcriptome_id.tsv"
    output:
        fasta = filtering + "filtered_transcriptome.fasta"
    log:
        filtering + "filtering_transcriptome_fasta.log"
    benchmark:
        filtering + "filtering_transcriptome_fasta.json"
    shell:
        "cat {input.ids} "
        "| xargs samtools faidx {input.assembly} "
        "> {output.fasta} "
        "2> {log}"



rule filtering_index_filtered_transcriptome:
    """
    Index the filtered fasta (samtools)
    """
    input:
        assembly = filtering + "filtered_transcriptome.fasta"
    output:
        index    = filtering + "filtered_transcriptome.idx"
    log:
        filtering + "index_filtered_transcriptome.log"
    benchmark:
        filtering + "index_filtered_transcriptome.json"
    shell:
        "kallisto index "
            "--index={output.index} "
            "{input.assembly} "
        "2> {log}"



rule filtering_faidx_filtered_transcriptome:
    input:
        fasta = filtering + "filtered_transcriptome.fasta"
    output:
        fai = filtering + "filtered_transcriptome.fasta.fai"
    log:
        filtering + "faidx_filtered_transcriptome.log"
    benchmark:
        filtering + "faidx_filtered_transcriptome.log"
    shell:
        "samtools faidx {input.fasta} 2> {log}"