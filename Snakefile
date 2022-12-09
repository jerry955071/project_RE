from tinydb import TinyDB, Query
from pathlib import Path


# Configure `mydb`
MYDB_PATH = Path("/home/b05b01002/HDD3/project_RE/mydb/")
mydb = TinyDB(MYDB_PATH / "db.json")
Q = Query()


# Query all records to process
_RECORDS = []
_SPECIES = ["Ptr", "Egr"]
for spe in _SPECIES:
    _RECORDS += mydb.search((Q.species == spe) & (Q.type == "RNA"))


# Checking db.search results
rule check_query:
    run:
        print("[RULE: check_query] Found the following files:")
        for r in _RECORDS:
            print("[RULE: check_query] " + r["id"])    
       

# Requesting all output files
rule all:
    input:
        [
            "processed/reditools2/%s/Final/%s/%s.txt" % (
            r["species"],
            r["tissue"],
            r["id"]
            ) for r in _RECORDS
        ]
    run:
        print("[RULE: ALL] Requesting the following outputs:")
        for fname in input:
            print(f"[RULE: ALL] {fname}")


# Identify RNA editing events using REDItools2.
# Following the tutorial on https://github.com/BioinfoUNIBA/REDItools2.
# ==============================================================================
# Step 1: REDITools2 for RNA-seq
rule reditools2_rna:
    """
    Following the tutorial on https://github.com/BioinfoUNIBA/REDItools2
    """
    conda: "envs/re.yml"
    threads: 1
    input:
        sam="processed/sorted_bam/{species}/RNA/{tissue}/{id}/Aligned.out.sam",
        ref="genomic_data/{species}/{species}.fa"
    output:
        "processed/reditools2/{species}/RNA/{tissue}/{id}.txt"
    params:
        bam="processed/sorted_bam/{species}/RNA/{tissue}/{id}.sorted.bam",
        strand=lambda wildcards: 1 if mydb.search(Q.id == wildcards.id)[0]["stranded"] else 0
    log:
        "logs/reditools2/{species}/RNA/{tissue}/{id}.log"
    shell:
        """
        # sam -> sorted.bam
        samtools sort \
            -@ {threads} \
            -O BAM \
            -o {params.bam} \
            {input.sam} \
            1>> {log} 2>> {log}

        # index bam file
        if [[ ! -e {params.bam}.bai ]];then
            samtools index -@ {threads} -b {params.bam}
        fi
        
        # docker reditools.py
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            ccc/reditools2:latest \
            /reditools2.0/src/cineca/reditools.py \
                -s {params.strand} \
                -f {params.bam} \
                -r {input.ref} \
                -o {output} > {log}
        """

# Step 2: REDITools2 for DNA-seq
rule reditools2_dna_parallel_singleRNA:
    conda: "envs/re.yml"
    threads: 32
    input:
        target_pos="processed/reditools2/{species}/RNA/{tissue}/{id}.bed",
        dna="processed/sorted_bam/{species}/DNA/merged.sorted.bam",
        cov_dna="processed/reditools2/{species}/DNA/coverage/merged.sorted.cov"
    output:
        "processed/reditools2/{species}/DNA/{tissue}/{id}.txt"
    params:
        ref="genomic_data/{species}/{species}.fa",
        size_file="genomic_data/{species}/{species}.fa.fai",
        tmp_dir="processed/reditools2/{species}/DNA/{tissue}/{id}/tmp/",
        cov_dir="processed/reditools2/{species}/DNA/coverage/"
    log:
        "logs/reditools2/{species}/DNA/dna.log"
    shell:
        """
        # reditools2
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
                    -r {params.ref} \
                    -f {input.dna} \
                    -B {input.target_pos} \
                    -G {input.cov_dna} \
                    -D {params.cov_dir} \
                    -Z {params.size_file} \
                    -t {params.tmp_dir} \
                    -o {output} > {log}

        # merge parallel outputs
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            ccc/reditools2:latest \
                bash /reditools2.0/merge.sh \
                {params.tmp_dir} \
                {output}.gz \
                {threads}
            
        # gunzip
        gunzip {output}.gz
        """


# # Step 2: REDITools2 for DNA-seq
# rule reditools2_dna_parallel_allRNA:
#     conda: "envs/re.yml"
#     threads: 32
#     input:
#         target_pos="processed/reditools2/{species}/RNA/target_pos.bed",
#         dna="processed/sorted_bam/{species}/DNA/merged.sorted.bam",
#         cov_dna="processed/reditools2/{species}/DNA/coverage/merged.sorted.cov"
#     output:
#         "processed/reditools2/{species}/DNA/dna.txt"
#     params:
#         ref="genomic_data/{species}/{species}.fa",
#         size_file="genomic_data/{species}/{species}.fa.fai",
#         tmp_dir="processed/reditools2/{species}/DNA/tmp/",
#         cov_dir="processed/reditools2/{species}/DNA/coverage/"
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
#         target_pos="processed/reditools2/{species}/RNA/{tissue}/{id}.bed",
#         dna="processed/sorted_bam/{species}/DNA/merged.sorted.bam",
#         cov_dna="processed/reditools2/{species}/DNA/coverage/merged.sorted.cov"
#     output:
#         "processed/reditools2/{species}/DNA/{tissue}/{id}.txt"
#     params:
#         ref="genomic_data/{species}/{species}.fa",
#         size_file="genomic_data/{species}/{species}.fa.fai",
#         tmp_dir="processed/reditools2/{species}/DNA/{tissue}/{id}/tmp/",
#         cov_dir="processed/reditools2/{species}/DNA/coverage/"
#     log:
#         "logs/reditools2/{species}/DNA/{tissue}/{id}.log"
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


# Step 3: annotate RNA.txt with DNA.txt
rule reditools2_annot:
    threads: 1
    input:
        rna="processed/reditools2/{species}/RNA/{tissue}/{id}.txt",
        dna="processed/reditools2/{species}/DNA/{tissue}/{id}.txt"
    output:
        "processed/reditools2/{species}/Final/{tissue}/{id}.txt"
    params:
        ref="genomic_data/{species}/{species}.fa.fai"
    log:
        "logs/reditools2/{species}/Final/{tissue}/{id}.log"
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
                    -R {params.ref} \
                    -r {input.rna} \
                    -d {input.dna} \
                    2> {log} 1> {output}
        """


# # Step 3: annotate RNA.txt with DNA.txt
# rule reditools2_annot:
#     threads: 1
#     input:
#         rna="processed/reditools2/{species}/RNA/{tissue}/{id}.txt",
#         dna="processed/reditools2/{species}/DNA/{tissue}/{id}.txt"
#     output:
#         "processed/reditools2/{species}/Final/{tissue}/{id}.txt"
#     params:
#         ref="genomic_data/{species}/{species}.fa.fai"
#     log:
#         "logs/reditools2/{species}/Final/{tissue}/{id}.log"
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


# RNA.txt -> RNA.bed
rule reditools_table_to_bed:
    threads: 1
    input:
        "processed/reditools2/{species}/RNA/{tissue}/{id}.txt"
    output:
        "processed/reditools2/{species}/RNA/{tissue}/{id}.bed"
    log:
        "logs/reditools2/{species}/to_bed/{tissue}/{id}.log"
    shell:
        """
        # awk
        cat {input} | \
            awk -v OFS='\\t' 'NR != 1 {{print}}' > {input}.tmp

        # docker run reditools
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            ccc/reditools2:latest \
            python /reditools2.0/src/cineca/reditools_table_to_bed.py \
                -i {input}.tmp \
                -o {output} > {log}

        # rm .tmp file
        rm {input}.tmp
        """


# merge RNA.bed files
rule merge_rna:
    threads: 4
    input:
        lambda wildcards:
            [
                "processed/reditools2/%s/RNA/%s/%s.bed" % (
                r["species"],
                r["tissue"],
                r["id"]
                ) for r in mydb.search(
                    (Q.species == wildcards.species) &
                    (Q.type == "RNA")
                )
            ]
    output:
        "processed/reditools2/{species}/RNA/target_pos.bed"
    log:
        "logs/bedtools/{species}.log"
    shell:
        """
        # cat and sort all input files (low-mem)
        touch {output}.tmp
        for file in {input}
        do 
            # cat and sort 2 bed files
            cat {output}.tmp "$file" | sort -k1,1 -k2,2n > {output}.tmp
        done

        # bedtools merge
        docker run \
            --cpus {threads} \
            --rm \
            -u $(id -u) \
            -v $(pwd):/data:rw \
            -w /data \
            -t \
            biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 \
            bedtools merge -i {output}.tmp 1> {output} 2> {log}
        
        # cleanup
        rm {output}.tmp
        """
    


# merge dna.sorted.bam files
rule merge_dna:
    threads: 32
    conda: "envs/re.yml"
    input:
        lambda wildcards: [
            f"processed/sorted_bam/{wildcards.species}/DNA/{i['id']}.sorted.bam" 
                for i in mydb.search(
                    (Q.species == wildcards.species) &
                    (Q.type == "DNA")
                )
            ]
    output:
        "processed/sorted_bam/{species}/DNA/merged.sorted.bam"
    shell:
        """
        # merge DNA.bam files
        samtools merge \
            -@ {threads} \
            {output} \
            {input} 
        
        # index merge.bam file
        samtools index -@ {threads} -b {output}
        """


# extract coverage from DNA files for `parallel_reditools.py`
rule extract_cov_dna:
    threads: 25
    input: "processed/sorted_bam/{species}/DNA/merged.sorted.bam"
    output: "processed/reditools2/{species}/DNA/coverage/merged.sorted.cov"
    params:
        size_file="genomic_data/{species}/{species}.fa.fai",
        cov_dir="processed/reditools2/{species}/DNA/coverage/"
    log:
        "logs/reditools2/{species}/DNA/extract_coverage.log"
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
                {input} {params.cov_dir} {params.size_file} \
                1> {log} 2> {log}
        """


# Align DNA-seq using BWA-MEM
# ==============================================================================
# BWA-MEM
rule bwa_mem:
    conda: "envs/re.yml"
    threads: 16
    input:
        qry=lambda wildcards: 
            [
                MYDB_PATH / i for i in 
                mydb.search(
                    (Q.species == wildcards.species) &
                    (Q.type == "DNA") &
                    (Q.id == wildcards.id)
                )[0]["fnames"]
            ]
    output:
        "processed/sorted_bam/{species}/DNA/{id}.sorted.bam"
    params:
        ref="genomic_data/{species}/{species}.fa"
    log:
        "logs/bwa-mem/{species}/{id}.log"
    shell:
        """
        # check index file exist
        if [[ ! -e {params.ref}.bwt ]]; then
            bwa index {params.ref} 1>> {log} 2>> {log}
        fi
        bwa mem \
            -t {threads} \
            {params.ref} \
            {input.qry} \
            2>> {log} | \
        samtools sort \
            -@ {threads} \
            -O BAM \
            -o {output} \
            1>> {log} 2>> {log}
        """


# Align RNA-seq data using STAR
# ==============================================================================
rule star:
    conda: "envs/re.yml"
    threads: 16
    input:
        qry=lambda wildcards: 
            [
                MYDB_PATH / i for i in 
                mydb.search(
                    (Q.species == wildcards.species) &
                    (Q.tissue == wildcards.tissue) &
                    (Q.type == "RNA") &
                    (Q.id == wildcards.id)
                )[0]["fnames"]
            ]
    output:
        "processed/sorted_bam/{species}/RNA/{tissue}/{id}/Aligned.out.sam"
    params:
        ref_fa="genomic_data/{species}/{species}.fa",
        ref_annt="genomic_data/{species}/{species}.gff3",
        out_prefix="processed/sorted_bam/{species}/RNA/{tissue}/{id}/",
        ovhang=lambda wildcards: 
            int(mydb.search(Q.id == wildcards.id)[0]["nbases"]) - 1
    log:
        "logs/star/{species}/{tissue}/{id}.log"
    shell:
        """
        # STAR indicing
        if [[ ! -e genomic_data/{wildcards.species}/star_{params.ovhang} ]]; then 
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
                    --genomeDir genomic_data/{wildcards.species}/star_{params.ovhang} \
                    --genomeFastaFiles {params.ref_fa} \
                    --sjdbGTFfile {params.ref_annt} \
                    --sjdbOverhang {params.ovhang} 1>> {log} 2>> {log}
        fi
        
        # STAR mapping
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
                --genomeDir genomic_data/{wildcards.species}/star_{params.ovhang} \
                --readFilesIn {input.qry} \
                --outFileNamePrefix {params.out_prefix} 1>> {log} 2>> {log}
        
        # # sam -> sorted.bam
        # samtools sort \
        #     -@ {threads} \
        #     -O BAM \
        #     -o {output} \
        #     {params.out_prefix}/Aligned.out.sam \
        #     1>> {log} 2>> {log}
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
