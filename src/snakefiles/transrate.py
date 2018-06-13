rule transrate_analysis:
    """
    Compare the raw and filtered assembly with transrate
    """
    input:
        pre = assembly + "Trinity.fasta",
        post = filtering + "filtered_transcriptome.fasta",
        forwards = expand(
            qc + "{sample}_1.fq",
            sample = SAMPLES_PE
        ),
        reverses = expand(
            qc + "{sample}_2.fq",
            sample = SAMPLES_PE
        )
    output:
        csv = transrate + "assemblies.csv"
    params:
        out_dir = transrate,
        forwards = ",".join(
            expand(
                qc + "{sample}_1.fq",
                sample = SAMPLES_PE
            ),
        ),
        reverses = ",".join(
            expand(
                qc + "{sample}_2.fq",
                sample = SAMPLES_PE
            ),
        )
    threads:
        24
    log:
        transrate + "analysis.log"
    benchmark:
        transrate + "analysis.json"
    shell:
        "docker run "
            "--rm "
            "--volume `pwd`:`pwd` "
            "--workdir `pwd` "
            "--user `id -u $USER`:`id -g $USER` "
            "ycogne/transrate "
                "/usr/bin/transrate-1.0.3-linux-x86_64/transrate "
                    "--assembly {input.pre},{input.post} "
                    "--left {params.forwards} "
                    "--right {params.reverses} "
                    "--output {params.out_dir} "
                    "--threads {threads} "
        "2> {log} 1>&2"





