rule assembly_split_pe_files:
    """
    Split pe_pe files into _1 and _2. Note the sleep: If it is not there, the assembly_merge_right_and_left fails.
    """
    input:
        fastq_pe = norm + "{sample}.final.pe_pe.fq.gz"
    output:
        left  = assembly + "{sample}_1.fq.gz",
        right = assembly + "{sample}_2.fq.gz"
    threads:
        ALL_THREADS
    priority:
        20
    params:
        left  = "{sample}.final.pe_pe.fq.gz.1",
        right = "{sample}.final.pe_pe.fq.gz.2"
    log:
        assembly + "split_pe_files_{sample}.log"
    benchmark:
        assembly + "split_pe_files_{sample}.json"
    shell:
        "split-paired-reads.py "
            "--output-first >(pigz --best > {output.left}) "
            "--output-second >(pigz --best > {output.right}) "
            "<(gzip --decompress --stdout {input.fastq_pe}) "
        "> {log} 2>&1 && "
        "sleep 5"


rule assembly_merge_right_and_left:
    """
    Generate the left.fq and right.fq
    left.fq = /1 reads from PE + all SE reads
    right.fq = /2 reads form PE
    Forced the use of gzip because pigz doesn't have --force option.
        --force is required because /dev/null isn't a file, doesn't even
        have an end, and isn't compressed.
    Test if SAMPLES_PE is empty because gzip/pigz may because otherwise
        it may be waiting something from stdout.
    """
    input:
        forward=  expand(
            assembly + "{sample}_1.fq.gz",
            sample = SAMPLES_PE
        ) if SAMPLES_PE  else ["/dev/null"],
        reverse=  expand(
            assembly + "{sample}_2.fq.gz",
            sample = SAMPLES_PE
        ) if SAMPLES_PE  else ["/dev/null"],
        unpaired= expand( # pe_se
            norm + "{sample}.final.pe_se.fq.gz",
            sample = SAMPLES_PE
        ) + expand( # se
            norm + "{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        )
    output:
        left=  assembly + "left.fq",
        right= assembly + "right.fq"
    threads:
        1
    priority:
        50
    log:
        assembly + "merge_right_and_left.log"
    benchmark:
        assembly + "merge_right_and_left.json"
    shell:
        "gzip --decompress --stdout {input.forward} {input.unpaired} > {output.left} 2> {log} && "
        "gzip --decompress --stdout {input.reverse} > {output.right} 2>> {log} "



rule assembly_run_trinity:
    """
    Assembly reads with Trinity.
    Notes on hardcoded settings:
        - Runs on paired end mode
        - Expect fastq files as inputs (left and right)
        - Does the full cleanup so it only remains a fasta file.
    """
    input:
        left  = assembly + "left.fq",
        right = assembly + "right.fq"
    output:
        fasta = protected(assembly + "Trinity.fasta")
    threads:
        ALL_THREADS
    priority:
        50
    params:
        memory= config["trinity_params"]["memory"],
        outdir= assembly + "trinity_out_dir"
    log:
        assembly + "run_trinity.log"
    benchmark:
        assembly + "run_trinity.json"
    shell:
        "Trinity "
            "--seqType fq "
            "--max_memory {params.memory} "
            "--left {input.left} "
            "--right {input.right} "
            "--CPU {threads} "
            "--full_cleanup "
            "--output {params.outdir} "
        "> {log} ; "
        "mv {params.outdir}.Trinity.fasta {output.fasta} ; "
        "rm {input.left}.readcount {input.right}.readcount"



rule assembly_gene_to_trans_map:
    """
    Create the gene_id TAB transcript_id file
    """
    input:
        fasta = assembly + "Trinity.fasta"
    output:
        tsv = assembly + "Trinity_gene_to_trans.tsv"
    log:
        assembly + "gene_to_trans_map.log"
    benchmark:
        assembly + "gene_to_trans_map.json"
    shell:
        "get_Trinity_gene_to_trans_map.pl "
        "< {input.fasta} "
        "> {output.tsv} "
        "2> {log}"



rule assembly_index_trinity:
    """
    Create a samtool index for the assembly
    """
    input:
        fasta = assembly + "Trinity.fasta"
    output:
        fai = assembly + "Trinity.fasta.fai"
    log:
        assembly + "index.log"
    benchmark:
        assembly + "index.json"
    shell:
        "samtools faidx {input.fasta} "
        "2> {log} 1>&2"