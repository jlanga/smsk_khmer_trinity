rule tissue_experimental_design:
    output:
        table = tissue + "experimental_design.tsv"
    threads:
        1
    log:
        tissue + "experimental_design.log"
    benchmark:
        tissue + "experimental_design.json"
    run:
        design_dict = {
            sample_name: [
                sample_name,
                "results/tissue/" + sample_name]
                for sample_name in config["samples_pe"]
        }
        with open(output.table, "w") as table:
            table.write("sample\ttissue\tpath\n")
            for sample_name in design_dict:
                table.write(
                    "\t".join(
                        [sample_name, design_dict[sample_name][0], design_dict[sample_name][1]]
                    ) + "\n"
                )



rule tissue_pseudoalign_sample:
    input:
        forward = qc + "{sample}_1.fq.gz",
        reverse = qc + "{sample}_2.fq.gz",
        index   = filtering + "filtered_transcriptome.idx"
    output:
        abundance_tsv = tissue + "{sample}/abundance.tsv",
        abundance_h5  = tissue + "{sample}/abundance.h5",
    params:
        outdir  = tissue + "{sample}/",
        params  = config["kallisto_params"]
    threads:
        24
    log:
        tissue + "{sample}/quant.log"
    benchmark:
        tissue + "{sample}/quant.json"
    shell:
        "kallisto quant "
            "--index={input.index} "
            "--output={params.outdir} "
            "--threads={threads} "
            "{params.params} "
            "<(pigz --decompress --stdout {input.forward}) "
            "<(pigz --decompress --stdout {input.reverse}) "
        "2> {log}"




rule tissue_normalised_tpms:
    input:
        tsv = expand(
            tissue + "{sample}/abundance.tsv",
            sample= config["samples_pe"]
        ),
        h5s = expand(
            tissue + "{sample}/abundance.h5",
            sample= config["samples_pe"]
        ),
        design = tissue + "experimental_design.tsv"
    output:
        rdata = tissue + "so.Rdata",
        tpms = tissue + "normalised_tpms.tsv"
    threads:
        1
    log:
        tissue + "normalised_tpms.log"
    benchmark:
        tissue + "normalised_tpms.json"
    shell:
        "Rscript src/get_normalised_tpms.R "
            "--experimental_design {input.design} "
            "--sleuth_object {output.rdata} "
            "--normalised_tpms {output.tpms} "
        "2> {log}"



rule tissue_extract_ids_sample:
    input:
        tpms = tissue + "normalised_tpms.tsv"
    output:
        tsv = tissue + "ids_{sample}.tsv"
    params:
        sample = "{sample}"
    log:
        tissue + "extract_id_{sample}.log"
    benchmark:
        tissue + "extract_id_{sample}.json"
    shell:
        "(tail -n +2 {input.tpms} "
        "| awk '{{if($4 > 0 && $2 == \"{params.sample}\") print $1 }}' "
        "> {output.tsv}) "
        "2> {log}"
