rule clean:
    shell:
        "if [ -e test.txt ]; then rm test.txt; fi"

rule all:
    output:
        "test.txt"
    shell:
        "echo Hello World! > {output}"
