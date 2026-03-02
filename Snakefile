import pandas as pd

configfile: "configs/config.yaml"

df = pd.read_csv(config["genomes_tsv"], sep="\t")
ACCESSIONS = df["Accession"].tolist()

rule all:
    input:
        "results/qc/assembly_qc.tsv",
        expand("results/proteins/{acc}/{acc}.faa", acc=ACCESSIONS),
        config["hmm_panel"] + ".h3p",
        expand("results/hmmer/{acc}/{acc}.tbl", acc=ACCESSIONS),
        "results/matrices/osmo_presence.tsv",
        "results/tables/osmo_freq_by_group.tsv",
        "results/tables/osmo_chi2_results.tsv",
        config["tableS3_out"],
        config["fig1_out"],
        config["fig2_out"]

rule assembly_qc:
    input:
        genomes=expand(config["genome_dir"] + "/{acc}/{acc}.fna", acc=ACCESSIONS)
    output:
        "results/qc/assembly_qc.tsv"
    shell:
        "mkdir -p results/qc && python workflow/scripts/assembly_qc.py {input.genomes} > {output}"

rule prodigal:
    input:
        fna=config["genome_dir"] + "/{acc}/{acc}.fna"
    output:
        faa="results/proteins/{acc}/{acc}.faa"
    shell:
        "mkdir -p results/proteins/{wildcards.acc} && prodigal -i {input.fna} -a {output.faa} -p single -q"

rule build_hmm_panel:
    input:
        osmo_kos=config["osmo_kos"]
    output:
        panel=config["hmm_panel"]
    shell:
        "python workflow/scripts/build_hmm_panel.py --osmo_kos {input.osmo_kos} --profiles {config[kofam_profiles]} --out {output.panel}"

rule hmmpress_panel:
    input:
        panel=config["hmm_panel"]
    output:
        flag=config["hmm_panel"] + ".h3p"
    shell:
        "hmmpress -f {input.panel}"

rule hmmsearch:
    input:
        panel=config["hmm_panel"],
        faa="results/proteins/{acc}/{acc}.faa",
        pressed=config["hmm_panel"] + ".h3p"
    output:
        tbl="results/hmmer/{acc}/{acc}.tbl",
        txt="results/hmmer/{acc}/{acc}.txt"
    shell:
        "mkdir -p results/hmmer/{wildcards.acc} && hmmsearch --tblout {output.tbl} --noali -E {config[hmmer_evalue]} {input.panel} {input.faa} > {output.txt}"

rule parse_matrix:
    input:
        tbls=expand("results/hmmer/{acc}/{acc}.tbl", acc=ACCESSIONS),
        osmo_kos=config["osmo_kos"],
        ko_list=config["ko_list"]
    output:
        "results/matrices/osmo_presence.tsv",
        "results/tables/osmo_hits_long.tsv"
    shell:
        "python workflow/scripts/parse_hmmer_to_matrix.py"

rule freq_by_group:
    input:
        mat="results/matrices/osmo_presence.tsv",
        genomes=config["genomes_tsv"]
    output:
        "results/tables/osmo_freq_by_group.tsv",
        "results/tables/osmo_counts_by_group.tsv"
    shell:
        "python workflow/scripts/freq_by_group.py"

rule chi2_fdr:
    input:
        mat="results/matrices/osmo_presence.tsv",
        genomes=config["genomes_tsv"]
    output:
        "results/tables/osmo_chi2_results.tsv"
    shell:
        "python workflow/scripts/chi2_fdr.py"

rule tableS3:
    input:
        chi="results/tables/osmo_chi2_results.tsv",
        freq="results/tables/osmo_freq_by_group.tsv"
    output:
        config["tableS3_out"]
    shell:
        "python workflow/scripts/tableS3.py"

rule figure1:
    input:
        mat="results/matrices/osmo_presence.tsv",
        genomes=config["genomes_tsv"]
    output:
        config["fig1_out"]
    shell:
        "python workflow/scripts/figure1.py"

rule figure2:
    input:
        chi="results/tables/osmo_chi2_results.tsv",
        freq="results/tables/osmo_freq_by_group.tsv"
    output:
        config["fig2_out"]
    shell:
        "python workflow/scripts/figure2.py"
