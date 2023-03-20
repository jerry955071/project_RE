from tinydb import TinyDB, Query
from pathlib import Path

# 
configfile: "configs/config.json"
wildcard_constraints:
    species="Ptr|Egr",
    id="\w*-\w*-\w*",
    filter_preset="\w*"

# configure `MYDB`
MYDB_PATH=Path(config["path_mydb"])
MYDB = TinyDB(config["mydb"])
# SPECIES = ["Ptr", "Egr"]
SPECIES = ["Ptr"]
Q = Query()

# parameter preset to filter reditable using awk
FILTER_PRESET = "set1"

# Requesting all output files
# ==============================================================================
rule all:
    input:
        [
            "outputs/%s/%s/%s.txt.gz.tbi" % (
            "reditools2" if record["stranded"] else "reditools",
            record["species"],
            record["id"]
            ) for record in MYDB.search(
                (Q.species.one_of(SPECIES)) & 
                (Q.type == "RNA")
            )
        ]
    run:
        print("[RULE: ALL] Requesting the following outputs:")
        for fname in input:
            print(f"[RULE: ALL] {fname}")


# Identify differential editing sites
# ==============================================================================
rule get_target_reditable:
    conda: "envs/re.yml"
    input:
        reditable=lambda wildcards: [
            "outputs/%s/%s/%s.txt.gz" % (
                "reditools2" if record["stranded"] else "reditools",
                record["species"],
                record["id"]
            ) for record in MYDB.search(
                (Q.id == wildcards.id)
            )
        ],
        target_positions="outputs/target_positions/{species}/{filter_preset}.txt"
    output:
        "outputs/target_reditable/{species}/{filter_preset}/{id}.txt"
    log:
        "logs/get_target_reditable/{species}/{filter_preset}/{id}.log"
    benchmark:
        "benchmarks/get_target_reditable/{species}/{filter_preset}/{id}.tsv"
    shell:
        """
        tabix -R {input.target_positions} {input.reditable} > {output} 2> {log}
        """
    

rule get_target_positions:
    input:
        lambda wildcards: [
            "outputs/%s/%s/%s.txt" % (
                "reditools2" if record["stranded"] else "reditools",
                record["species"],
                record["id"]
            ) for record in MYDB.search(
                (Q.species == wildcards.species) &
                (Q.type == "RNA")
            )
        ]
    output:
        outdir=directory("outputs/target_positions/{species}/{filter_preset}/"),
        fout="outputs/target_positions/{species}/{filter_preset}.txt"
    params:
        lambda wildcards: config["filter_presets"][wildcards.filter_preset]
    log:
        "logs/get_target_positions/{species}/{filter_preset}.log"
    benchmark:
        "benchmarks/get_target_positions/{species}/{filter_preset}.tsv"
    shell:
        """
        echo '' > {log}
        mkdir {output.outdir} > {log} 2>&1
        for f in {input}
        do
            out_tmp={output.outdir}/$(basename $f).tmp
            awk \
                -v FS='\t' \
                -v OFS='\t' \
                '{params}' \
                $f > $out_tmp &
            echo PID: $! $f >> {log}
        done
        wait

        cat {output.outdir}/*.tmp \
            | sort -k1,1 -k2,2n \
            | uniq > {output.fout}
        """

# Tabix indicing
# ==============================================================================
# tabix index REDItables
rule tabix_reditable:
    conda: "envs/re.yml"
    threads: 1
    input:
        "{fname}.txt"
    output:
        gz="{fname}.txt.gz",
        tbi="{fname}.txt.gz.tbi"
    log:
        "logs/tabix_reditable/{fname}.log"
    benchmark:
        "benchmarks/tabix_reditable/{fname}.tsv"
    shell:
        """
        bgzip < {input} > {output.gz} \
            && tabix -s 1 -b 2 -e 2 -c 'R' {output.gz} >> {log}
        """


# tabix index GFF
rule tabix_gff:
    conda: 
        "envs/re.yml"
    threads: 
        1
    input: 
        "{fname}.gff3"
    output:
        gz="{fname}.gff3.gz",
        tbi="{fname}.gff3.gz.tbi"
    log: 
        "logs/tabix_gff/{fname}.log"
    benchmark:
        "benchmarks/tabix_gff/{fname}.tsv"
    shell:
        """
        # NOTE: grep -v  ^"#" | sort | bgzip | tabix
        grep ^"#" {input} > {log}
        grep -v ^"#" {input} \
            | sort -k1,1 -k4,4n \
            | bgzip 1> {output.gz} 2>> {log} \
            && tabix -p gff {output.gz} >> {log}
        """


# Identify RNA editing events using REDItools 1 & 2.
# Following the tutorial on https://github.com/BioinfoUNIBA/REDItools2.
# ==============================================================================
# Using REDItools v1.3 for unstranded data:
rule reditools1:
    threads: 8
    input:
        dna="outputs/sorted_bam/{species}/DNA/{species}.sorted.bam",
        rna="outputs/sorted_bam/{species}/RNA/{id}.sorted.bam",
        asm="genomic_data/{species}/{species}.fa",
        ant="genomic_data/{species}/{species}.gff3.gz"
    output:
        outdir=directory("outputs/reditools/{species}/{id}/"),
        reditable="outputs/reditools/{species}/{id}.txt"
    params:
        mapq_dna=config["mapq_cutoff"]["bwa"],
        mapq_rna=config["mapq_cutoff"]["star"]
    log:
        "logs/reditools/{species}/{id}.log"
    benchmark:
        "benchmarks/reditools/{species}/{id}.tsv"
    shell:
        """
        # NOTE: REDItoolDnaRna.py
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            --name reditools_{wildcards.id} \
            ccc/reditool:v1.3 \
                REDItoolDnaRna.py \
                -i {input.rna} \
                -j {input.dna} \
                -f {input.asm} \
                -t {threads} \
                -o {output.outdir} \
                -m {params.mapq_dna},{params.mapq_rna} \
                -G {input.ant} \
                2>&1 \
                > {log} && \
        cp {output.outdir}/DnaRna_*/outTable_* {output.reditable}  
        """


# ==============================================================================
# Using REDITools 2 for stranded data:
# Step 1: REDITools2 for RNA-seq
rule reditools2_RNA:
    conda: "envs/re.yml"
    threads: 1
    input:
        rna="outputs/sorted_bam/{species}/RNA/{id}.sorted.bam",
        asm="genomic_data/{species}/{species}.fa"
    output:
        "outputs/reditools2/{species}/RNA/{id}.txt"
    params:
        strand=1
    log:
        "logs/reditools2/{species}/RNA/{id}.log"
    benchmark:
        "benchmarks/reditools2/{species}/RNA/{id}.tsv"
    shell:
        """
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            --name reditools_rna_{wildcards.id} \
            ccc/reditools2:latest \
                /reditools2.0/src/cineca/reditools.py \
                -f {input.rna} \
                -o {output} \
                -r {input.asm} \
                -s {params.strand} \
                -H \
                > {log}
        """



# Step 2: REDITools2 for DNA-seq
# NOTE: passing directory as input
rule parallel_reditools2_DNA_with_merged_RNA:
    conda: "envs/re.yml"
    threads: 32
    input:
        fa="genomic_data/{species}/{species}.fa",
        fai="genomic_data/{species}/{species}.fa.fai",
        rna="outputs/reditools2/{species}/RNA/{species}.merged.bed",
        dna="outputs/sorted_bam/{species}/DNA/{species}.sorted.bam",
        cov_dna="outputs/reditools2/{species}/DNA/coverage/{species}.sorted.cov",
        cov_dna_dir="outputs/reditools2/{species}/DNA/coverage/"
    output:
        tmp_dir=directory("outputs/reditools2/{species}/DNA/tmp/"),
    log:
        "logs/reditools2/{species}/DNA/parallel_reditools_dna.log",
    benchmark:
        "benchmarks/reditools2/{species}/DNA/parallel_reditools_dna.log"
    shell:
        """
        # NOTE: parallel_reditools.py --dna
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            ccc/reditools2:latest \
            mpirun -np {threads} \
                /reditools2.0/src/cineca/parallel_reditools.py \
                    --dna \
                    -r {input.fa} \
                    -f {input.dna} \
                    -B {input.rna} \
                    -G {input.cov_dna} \
                    -D {input.cov_dna_dir}/ \
                    -Z {input.fai} \
                    -t {output.tmp_dir}/ \
                    -H > {log}
        """

rule merge_parallel_reditools2_DNA_with_merged_RNA:
    conda: "envs/re.yml"
    threads: 2
    input:
        tmp_dir="outputs/reditools2/{species}/DNA/tmp/"
    output:
        output="outputs/reditools2/{species}/DNA/merged_RNA_DNA.txt.gz"
    log:
        "logs/reditools2/{species}/DNA/merge.log"
    benchmark:
        "benchmarks/reditools2/{species}/DNA/merge.log"
    shell:
        """
        # NOTE: /reditools2.0/merge.sh
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            ccc/reditools2:latest \
                bash /reditools2.0/merge.sh \
                {input.tmp_dir} \
                {output.output} \
                {threads} > {log}
        """



# Step 3: annotate RNA.txt with merged_RNA_DNA.txt.gz
rule reditools2_annot_with_merged_RNA_DNA:
    threads: 1
    input:
        rna="outputs/reditools2/{species}/RNA/{id}.txt",
        dna="outputs/reditools2/{species}/DNA/merged_RNA_DNA.txt.gz",
        fai="genomic_data/{species}/{species}.fa.fai"
    output:
        "outputs/reditools2/{species}/{id}.txt"
    log:
        "logs/reditools2/{species}/annot/{id}.log"
    benchmark:
        "benchmarks/reditools2/{species}/annot/{id}.tsv"
    shell:
        """
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            --name reditools_annot_{wildcards.id} \
            ccc/reditools2:latest \
                python /reditools2.0/src/cineca/annotate_with_DNA.py \
                    -R {input.fai} \
                    -r {input.rna} \
                    -d {input.dna} \
                    -H \
                    2> {log} 1> {output}
        """


# # Step 2: REDITools2 for DNA-seq
# rule reditools2_dna_parallel_singleRNA:
#     conda: "envs/re.yml"
#     threads: 16
#     input:
#         target_pos="outputs/reditools2/{species}/RNA/{id}.bed",
#         dna="outputs/sorted_bam/{species}/DNA/merged.sorted.bam",
#         cov_dna="outputs/reditools2/{species}/DNA/coverage/merged.sorted.cov"
#     output:
#         "outputs/reditools2/{species}/DNA/{id}.txt"
#     params:
#         ref="genomic_data/{species}/{species}.fa",
#         size_file="genomic_data/{species}/{species}.fa.fai",
#         tmp_dir=temp("outputs/reditools2/{species}/DNA/{id}/tmp/"),
#         cov_dir="outputs/reditools2/{species}/DNA/coverage/"
#     log:
#         "logs/reditools2/{species}/DNA/{id}.log"
#     benchmark:
#         "benchmarks/reditools2/{species}/DNA/{id}.tsv"
#     shell:
#         """
#         # NOTE: reditools2
#         docker run \
#             --cpus {threads} \
#             --rm \
#             -u $(id -u) \
#             -v $(pwd):/data:rw \
#             -w /data \
#             -t \
#             --name reditools_parallel_dna_{wildcards.id} \
#             ccc/reditools2:latest \
#             mpirun -np {threads} \
#                 /reditools2.0/src/cineca/parallel_reditools.py \
#                 --dna \
#                 -r {params.ref} \
#                 -f {input.dna} \
#                 -B {input.target_pos} \
#                 -G {input.cov_dna} \
#                 -D {params.cov_dir} \
#                 -Z {params.size_file} \
#                 -t {params.tmp_dir} \
#                 -o {output} > {log}

#         # NOTE: merge parallel outputs
#         docker run \
#             --cpus {threads} \
#             --rm \
#             -u $(id -u) \
#             -v $(pwd):/data:rw \
#             -w /data \
#             -t \
#             --name reditools_merge_{wildcards.id} \
#             ccc/reditools2:latest \
#                 bash /reditools2.0/merge.sh \
#                 {params.tmp_dir} \
#                 {output}.gz \
#                 {threads}
            
#         # NOTE: gunzip
#         gunzip {output}.gz
#         """
#
#
# # Step 3: annotate RNA.txt with DNA.txt
# rule reditools2_annot:
#     threads: 1
#     input:
#         rna="outputs/reditools2/{species}/RNA/{id}.txt",
#         dna="outputs/reditools2/{species}/DNA/{id}.txt",
#         fai="genomic_data/{species}/{species}.fa.fai"
#     output:
#         "outputs/reditools2/{species}/{id}.txt"
#     log:
#         "logs/reditools2/{species}/Final/{id}.log"
#     benchmark:
#         "benchmarks/reditools2/{species}/Final/{id}.tsv"
#     shell:
#         """
#         docker run \
#             --cpus {threads} \
#             --rm \
#             -u $(id -u) \
#             -v $(pwd):/data:rw \
#             -w /data \
#             --name reditools_annot_{wildcards.id} \
#             ccc/reditools2:latest \
#                 python /reditools2.0/src/cineca/annotate_with_DNA.py \
#                     -R {input.fai} \
#                     -r {input.rna} \
#                     -d {input.dna} \
#                     -H \
#                     2> {log} 1> {output}
#         """


# Accessory rules to reditools 1&2
# ==============================================================================
# merge RNA.bed files
rule get_merge_rna_bed:
    threads: 4
    input:
        lambda wildcards:
            [
                "outputs/reditools2/%s/RNA/%s.bed" % (
                r["species"],
                r["id"]
                ) for r in MYDB.search(
                    (Q.species == wildcards.species) &
                    (Q.type == "RNA") &
                    (Q.stranded == True)
                )
            ]
    output:
        tmp=temp("outputs/reditools2/{species}/RNA/{species}.bed.tmp"),
        merged="outputs/reditools2/{species}/RNA/{species}.merged.bed"
    log:
        "logs/merge_rna/{species}.log"
    benchmark:
        "benchmarks/merge_rna/{species}.tsv"
    shell:
        """
        # NOTE: cat and sort all input files (low-mem)
        touch {output.tmp}
        for file in {input}
        do 
            # NOTE: cat and sort 2 bed files
            cat {output.tmp} "$file" | sort -k1,1 -k2,2n > {output.tmp}
        done

        # NOTE: bedtools merge
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 \
            bedtools merge -i {output.tmp} 1> {output.merged} 2> {log}
        """
    


# reditable to bed
rule reditable2bed:
    threads: 1
    input:
        "outputs/reditools2/{species}/RNA/{id}.txt"
    output:
        "outputs/reditools2/{species}/RNA/{id}.bed"
    log:
        "logs/reditools2/reditable2bed/{species}/{id}.log"
    benchmark:
        "benchmarks/reditools2/reditable2bed/{species}/{id}.tsv"
    shell:
        """
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            --name reditable2bed_{wildcards.id} \
            ccc/reditools2:latest \
            python /reditools2.0/src/cineca/reditools_table_to_bed.py \
                -i {input} \
                -o {output} > {log}
        """


# merge dna.sorted.bam files
rule get_merge_dna_sorted_bam:
    threads: 32
    conda: "envs/re.yml"
    input:
        lambda wildcards: [
            f"outputs/sorted_bam/{wildcards.species}/DNA/{i['id']}.sorted.bam" 
                for i in MYDB.search(
                    (Q.species == wildcards.species) &
                    (Q.type == "DNA")
                )
            ]
    output:
        "outputs/sorted_bam/{species}/DNA/{species}.sorted.bam"
    log:
        "logs/merge_dna/{species}.log"
    benchmark:
        "benchmarks/merge_dna/{species}.tsv"
    shell:
        """
        # NOTE: merge DNA.bam files
        samtools merge -@ {threads} {output} {input} 
        
        # NOTE: index merge.bam file
        samtools index -@ {threads} -b {output}
        """


# extract coverage from DNA files for `parallel_reditools.py`
rule extract_cov_merge_dna_sorted_bam:
    threads: 25
    input: 
        merged="outputs/sorted_bam/{species}/DNA/{species}.sorted.bam",
        fai="genomic_data/{species}/{species}.fa.fai"
    output: 
        cov="outputs/reditools2/{species}/DNA/coverage/{species}.sorted.cov",
        cov_dir=directory("outputs/reditools2/{species}/DNA/coverage/")
    log:
        "logs/reditools2/{species}/DNA/extract_coverage.log"
    benchmark:
        "benchmarks/reditools2/{species}/DNA/extract_coverage.tsv"
    shell:
        """
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            ccc/reditools2:latest \
            bash /reditools2.0/extract_coverage.sh \
                {input.merged} {output.cov_dir}/ {input.fai} 1> {log} 2> {log}
        """


# Align DNA-seq using BWA-MEM
# ==============================================================================
rule qc_bwa:
    conda: "envs/re.yml"
    threads: 2
    input:
        "outputs/sorted_bam/{species}/DNA/{id}.bam"
    output:
        "outputs/sorted_bam/{species}/DNA/{id}.sorted.bam"
    params:
        mapq=config["mapq_cutoff"]["bwa"]
    log:
        "logs/qc_bwa/{species}/{id}.log"
    benchmark:
        "benchmarks/qc_bwa/{species}/{id}.tsv"
    shell:
        """
        samtools view \
            -b \
            -h \
            -q {params.mapq} \
            -F 0x100 \
            -F 0x800 \
            -@ {threads} \
            {input} 2> {log} \
        | samtools sort \
            -o {output} \
            -O BAM \
            -@ {threads} 1>> {log} 2>> {log} 
        """

# NOTE: passing directory as input
rule bwa_mem_align:
    conda: "envs/re.yml"
    threads: 16
    input:
        fqdir="outputs/fastp/{id}/",
        ref="genomic_data/{species}/{species}.fa"
    output:
        "outputs/sorted_bam/{species}/DNA/{id}.bam"
    log:
        "logs/bwa_mem/{species}/{id}.log"
    benchmark:
        "benchmarks/bwa_mem/{species}/{id}.tsv"
    shell:
        """
        # NOTE: bwa-mem | samtools-view | samtools-sort
        bwa mem \
            -t {threads} \
            {input.ref} \
            $( ls {input.fqdir}/*.fq ) 2>> {log} \
        | samtools view \
            -b \
            -h \
            -o {output} \
            -@ {threads} 2>> {log}
        """


# rule bwa_mem_index:
#     conda: "envs/re.yml"
#     threads: 16
#     input:
#         fnames=lambda wildcards: 
#             [
#                 MYDB_PATH / i for i in 
#                 MYDB.search(
#                     (Q.species == wildcards.species) &
#                     (Q.type == "DNA") &
#                     (Q.id == wildcards.id)
#                 )[0]["fnames"]
#             ],
#         ref="genomic_data/{species}/{species}.fa"
#     output:
#         "outputs/sorted_bam/{species}/DNA/{id}.sorted.bam"
#     params:
#         mapq=config["mapq_cutoff"]["bwa"]
#     log:
#         "logs/bwa-mem/{species}/{id}.log"
#     benchmark:
#         "benchmarks/bwa_mem/{species}/{id}.tsv"
#     shell:
#         """
#         # NOTE: check index file exist
#         if [[ ! -e {params.ref}.bwt ]]; then
#             bwa index {params.ref} 1>> {log} 2>> {log}
#         fi

#         # NOTE: bwa-mem | samtools-view | samtools-sort
#         bwa mem \
#             -t {threads} \
#             {input.ref} \
#             {input.fnames} \
#             2>> {log} \
#         | samtools view \
#             --bam \
#             --with-header \
#             --min-MQ {params.mapq} \
#             --exclude-flags 0x100 \
#             --exclude-flags 0x800 \
#             -@ {threads} \
#         | samtools sort \
#             -o {output} \
#             -O BAM \
#             -@ {threads} \
#             1>> {log} 2>> {log}
#         """


# Align RNA-seq data to genome using STAR
# ==============================================================================
rule qc_star:
    conda: "envs/re.yml"
    threads: 2
    input:
        "outputs/sorted_bam/{species}/RNA/{id}/Aligned.out.sam"
    output:
        "outputs/sorted_bam/{species}/RNA/{id}.sorted.bam"
    params:
        mapq=config["mapq_cutoff"]["star"]
    log:
        "logs/qc_star/{species}/{id}.log"
    benchmark:
        "benchmarks/qc_star/{species}/{id}.tsv"
    shell:
        """
        samtools view \
            -b \
            -h \
            -q {params.mapq} \
            -F 0x100 \
            -F 0x800 \
            -@ {threads} \
            {input} 2> {log} \
        | samtools sort \
            -o {output} \
            -O BAM \
            -@ {threads} 1>> {log} 2>> {log} \
        && samtools index {output} 1>> {log} 2>> {log}
        """


rule star_align:
    conda: "envs/re.yml"
    threads: 16
    input:
        fqdir="outputs/fastp/{id}",
        asm="genomic_data/{species}/{species}.fa",
        fai="genomic_data/{species}/{species}.fa.fai",
        ant="genomic_data/{species}/{species}.gff3",
        idx=lambda wildcards:
            "genomic_data/%s/star_%s" % (
                wildcards.species,
                int(MYDB.search(Q.id == wildcards.id)[0]["nbases"]) - 1
            )
    output:
        outdir=directory("outputs/sorted_bam/{species}/RNA/{id}/"),
        outsam="outputs/sorted_bam/{species}/RNA/{id}/Aligned.out.sam"
    log:
        "logs/star/align/{species}/{id}.log"
    benchmark:
        "benchmarks/star/align/{species}/{id}.tsv"
    shell:
        """
        # NOTE: STAR-align
        docker run \
            --rm \
            --cpus {threads} \
            -v $(pwd):/data:rw \
            -v {MYDB_PATH}:{MYDB_PATH}:rw \
            -w /data \
            -u $(id -u):$(id -g) \
            ccc/star \
            STAR \
                --runThreadN {threads} \
                --genomeDir {input.idx} \
                --readFilesIn $( ls {input.fqdir}/*.fq ) \
                --outFileNamePrefix {output.outdir}/ > {log}
        """


rule star_index:
    threads: 1
    input:
        asm="genomic_data/{species}/{species}.fa",
        ant="genomic_data/{species}/{species}.gff3"
    output:
        directory("genomic_data/{species}/star_{overhang}")
    log:
        "logs/star/index/{species}_{overhang}.log"
    benchmark:
        "benchmarks/star/{species}_{overhang}.tsv"
    shell:
        """
        # NOTE: STAR-index
        docker run \
            --rm \
            --cpus {threads} \
            -v $(pwd):/working:rw \
            -w /working \
            -u $(id -u):$(id -g) \
            ccc/star \
            STAR \
                --runThreadN {threads} \
                --runMode genomeGenerate \
                --genomeDir {output} \
                --genomeFastaFiles {input.asm} \
                --sjdbGTFfile {input.ant} \
                --sjdbOverhang {wildcards.overhang} 1>> {log} 2>> {log}
        """


# fastp
# ====================================================================
def get_fastp_params(wildcards):
    records =  MYDB.search(Q.id == wildcards.id)[0]
    id=records["id"]
    f1=records["fnames"][0]
    f2=records["fnames"][1] if len(records["fnames"]) == 2 else None
    # 
    if records["layout"] == "paired":
        return f"""\
            --in1 {MYDB_PATH}/{f1} \
            --in2 {MYDB_PATH}/{f2} \
            --out1 outputs/fastp/{id}/{id}_1.fq \
            --out2 outputs/fastp/{id}/{id}_2.fq \
            --detect_adapter_for_pe \
            --cut_front \
            --cut_tail \
            --correction \
            --json outputs/fastp/{id}/report.json \
            --html outputs/fastp/{id}/report.html\
        """
    if records["layout"] == "single":
        return f"""\
            --in1 {MYDB_PATH}/{f1} \
            --out1 outputs/fastp/{id}/{id}.fq \
            --cut_front \
            --cut_tail \
            --json outputs/fastp/{id}/report.json \
            --html outputs/fastp/{id}/report.html \
        """ 


rule fastp:
    threads: 1
    conda:
        "envs/re.yml"
    input:
        lambda wildcards: [
            MYDB_PATH / i for i in MYDB.search(
                (Q.id == wildcards.id)
            )[0]["fnames"]
        ]
    output:
        directory("outputs/fastp/{id}/")
    params:
        get_fastp_params
    log:
        "logs/fastp/{id}.log"
    benchmark:
        "benchmarks/fastp/{id}.tsv"
    shell:
        """
        mkdir -p {output}
        fastp {params} 1> {log} 2> {log}
        """

# ====================================================================
# template
rule template:
    conda:
        "envs/re.yml"
    input:
    output:
    log:
    shell:
        """
        echo 'Hello world!'
        """


# # Step 2: REDITools2 for DNA-seq
# rule reditools2_dna_parallel_allRNA:
#     conda: "envs/re.yml"
#     threads: 32
#     input:
#         target_pos="outputs/reditools2/{species}/RNA/target_pos.bed",
#         dna="outputs/sorted_bam/{species}/DNA/merged.sorted.bam",
#         cov_dna="outputs/reditools2/{species}/DNA/coverage/merged.sorted.cov"
#     output:
#         "outputs/reditools2/{species}/DNA/dna.txt"
#     params:
#         ref="genomic_data/{species}/{species}.fa",
#         size_file="genomic_data/{species}/{species}.fa.fai",
#         tmp_dir="outputs/reditools2/{species}/DNA/tmp/",
#         cov_dir="outputs/reditools2/{species}/DNA/coverage/"
#     log:
#         "logs/reditools2/{species}/DNA/dna.log"
#     shell:
#         """
#         # reditools2
#         docker run \
#             --cpus {threads} \
#             --rm \
#             -u $(id -u) \
#             -v $(pwd):/data:rw \
#             -w /data \
#             -t \
#             ccc/reditools2:latest \
#             mpirun -np {threads} \
#                 /reditools2.0/src/cineca/parallel_reditools.py \
#                     --dna \
#                     -r {params.ref} \
#                     -f {input.dna} \
#                     -B {input.target_pos} \
#                     -G {input.cov_dna} \
#                     -D {params.cov_dir} \
#                     -Z {params.size_file} \
#                     -t {params.tmp_dir} \
#                     -o {output} > {log}

#         # merge parallel outputs
#         docker run \
#             --cpus {threads} \
#             --rm \
#             -u $(id -u) \
#             -v $(pwd):/data:rw \
#             -w /data \
#             -t \
#             ccc/reditools2:latest \
#                 bash /reditools2.0/merge.sh \
#                 {params.tmp_dir} \
#                 {output}.gz \
#                 {threads}
            
#         # gunzip
#         gunzip {output}.gz
#         """


# # Step 2: REDITools2 for DNA-seq
# rule reditools2_dna:
#     conda: "envs/re.yml"
#     threads: 24
#     input:
#         target_pos="outputs/reditools2/{species}/RNA/{id}.bed",
#         dna="outputs/sorted_bam/{species}/DNA/merged.sorted.bam",
#         cov_dna="outputs/reditools2/{species}/DNA/coverage/merged.sorted.cov"
#     output:
#         "outputs/reditools2/{species}/DNA/{id}.txt"
#     params:
#         ref="genomic_data/{species}/{species}.fa",
#         size_file="genomic_data/{species}/{species}.fa.fai",
#         tmp_dir="outputs/reditools2/{species}/DNA/{id}/tmp/",
#         cov_dir="outputs/reditools2/{species}/DNA/coverage/"
#     log:
#         "logs/reditools2/{species}/DNA/{id}.log"
#     shell:
#         """
#         # remove tmp_dir if exist
#         if [[ -e {params.tmp_dir} ]]; then
#             rm -r {params.tmp_dir}
#         fi

#         # reditools2
#         docker run \
#             --cpus {threads} \
#             --rm \
#             -u $(id -u) \
#             -v $(pwd):/data:rw \
#             -w /data \
#             -t \
#             ccc/reditools2:latest \
#             mpirun -np {threads} \
#                 /reditools2.0/src/cineca/parallel_reditools.py \
#                     --dna \
#                     -r {params.ref} \
#                     -f {input.dna} \
#                     -B {input.target_pos} \
#                     -G {input.cov_dna} \
#                     -D {params.cov_dir} \
#                     -Z {params.size_file} \
#                     -t {params.tmp_dir} \
#                     -o {output} > {log}

#         # merge parallel outputs
#         docker run \
#             --cpus {threads} \
#             --rm \
#             -u $(id -u) \
#             -v $(pwd):/data:rw \
#             -w /data \
#             -t \
#             ccc/reditools2:latest \
#                 bash /reditools2.0/merge.sh \
#                 {params.tmp_dir} \
#                 {output}.gz \
#                 {threads}
            
#         # gunzip
#         gunzip {output}.gz
#         """


# # Step 3: annotate RNA.txt with DNA.txt
# rule reditools2_annot:
#     threads: 1
#     input:
#         rna="outputs/reditools2/{species}/RNA/{id}.txt",
#         dna="outputs/reditools2/{species}/DNA/{id}.txt"
#     output:
#         "outputs/reditools2/{species}/Final/{id}.txt"
#     params:
#         ref="genomic_data/{species}/{species}.fa.fai"
#     log:
#         "logs/reditools2/{species}/Final/{id}.log"
#     shell:
#         """
#         docker run \
#             --cpus {threads} \
#             --rm \
#             -u $(id -u) \
#             -v $(pwd):/data:rw \
#             -w /data \
#             --name reditools_annot_{wildcards.id} \
#             ccc/reditools2:latest \
#                 python /reditools2.0/src/cineca/annotate_with_DNA.py \
#                     -R {params.ref} \
#                     -r {input.rna} \
#                     -d {input.dna} \
#                     2> {log} 1> {output}
#         """
