
rule protein_all:
    input:
        config['lr']['cpat']['protein'],
        config['lr']['cpat']['orf'],
        config['lr']['cpat']['protein_fa'],
        config['lr']['sqanti_protein']['data_exons'],
        config['lr']['sqanti_protein']['ref_cds_renamed'],
        config['lr']['sqanti_protein']['best_orf'],
        config['lr']['sqanti_protein']['classified'],
        config['lr']['protein']['summary'],
        config['lr']['blast']['out']


# rule download_all:
#     output:
#         human_novel="data-raw/human.gtf",
#         human_ref="data-raw/human_annotation.gtf",
#         human_genome="data-raw/human_genome.fa",
#         human_ref_translations="data-raw/gencode.v29.pc_translations.fa",
#         human_hexamer=config['ref']['cpat']['hexamer'],
#         human_logit_model=config['ref']['cpat']['model'],
#         mouse_novel="data-raw/mouse.gtf",
#         mouse_ref="data-raw/mouse_annotation.gtf",
#         mouse_genome="data-raw/mouse_genome.fa",
#         mouse_ref_translations="data-raw/gencode.vM21.pc_translations.fa",
#         mouse_hexamer="data-raw/Mouse_Hexamer.tsv",
#         mouse_logit_model="data-raw/Mouse_logitModel.RData",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     shell:
#         """
#         wget https://zenodo.org/records/10407864/files/cerberus.gtf;
#         mv cerberus.gtf {output.human_novel};
#         wget https://zenodo.org/records/10408250/files/cerberus.gtf;
#         mv cerberus.gtf {output.mouse_novel};
#         wget -O - https://www.encodeproject.org/files/ENCFF991WIA/@@download/ENCFF991WIA.gtf.gz | gunzip -c > {output.human_ref};
#         wget -O - https://www.encodeproject.org/files/gencode.vM21.primary_assembly.annotation_UCSC_names/@@download/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz | gunzip -c > {output.mouse_ref};
#         wget -O - https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz | gunzip -c > {output.human_genome};
#         wget -O - https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz | gunzip -c >  {output.mouse_genome};
#         wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_translations.fa.gz | gunzip -c >  {output.human_ref_translations};
#         wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.pc_translations.fa.gz | gunzip -c >  {output.mouse_ref_translations};
#         wget -O - https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Human_Hexamer.tsv >  {output.human_hexamer};
#         wget -O - https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Mouse_Hexamer.tsv >  {output.mouse_hexamer};
#         wget -O - https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Human_logitModel.RData  >  {output.human_logit_model};
#         wget -O - https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Mouse_logitModel.RData >  {output.mouse_logit_model};
#         """


rule prep_make_query:
    input:
        query=config['lr']['gtf_no_spike'],
        genome=config['ref']['fa'],
    output:
        config['lr']['orfanage']['query'],
        config['lr']['gtf_no_spike'],
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        cat {input.query} | gffread -g {input.genome} -T -o {output} -
        """


# rule prep_filter_spikeins_annotation:
#     input:
#         config['ref']['gtf'],
#     output:
#         "data-raw/{sample}_annotation_filtered.gtf",
#     container:
#         f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
#     shell:
#         """
#         awk '$1 ~ "chr"' {input} > {output}
#         """


rule prep_make_annotation:
    input:
        annotation=config['ref']['gtf'],
        genome=config['ref']['fa'],
    output:
        config['lr']['orfanage']['ready_annot']
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        cat {input.annotation} | gffread -g {input.genome} --adj-stop -T -F -J -o {output}
        """


rule orf_prediction_run_orfanage:
    input:
        query=config['lr']['gtf_no_spike'],
        genome=config['ref']['fa'],
        annotation=config['lr']['orfanage']['ready_annot'],
    output:
        config['lr']['orfanage']['cds']
    resources:
        threads = 8,
        nodes = 2
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/orfanage
        orfanage \
                    --cleanq \
                    --mode LONGEST_MATCH \
                    --reference {input.genome} \
                    --query {input.query} \
                    --output {output} \
                    --threads {threads} \
                    {input.annotation} \
                    1>../../data/protein/transcripts_orfanage.output \
                    2>../../data/protein/transcripts_orfanage.error"""


rule orf_prediction_filter_orfanage:
    input:
        orfanage_cds=config['lr']['orfanage']['cds'],
    output:
        to_be_predicted=config['lr']['orfanage']['cds_to_be_predicted'],
        filtered=config['lr']['orfanage']['cds_filtered'],
    params:
        min_orf_len=config['params']['orfanage']['min_orf_len_nt'],
        scripts_dir = p
    resources:
        threads = 1,
        nodes = 1
    shell:
        """python {params.scripts_dir}/filter_orfanage.py \
            --orfanage_gtf_file_path {input.orfanage_cds} \
            --output_path_to_be_predicted {output.to_be_predicted} \
            --output_path_filtered {output.filtered} \
            --minimum_orf_length {params.min_orf_len}"""


rule orf_prediction_orf_sequence_orfanage_with_stop_codon:
    input:
        orfanage_cds=config['lr']['orfanage']['cds_filtered'],
        genome=config['ref']['fa'],
    output:
        config['lr']['orfanage']['orfs']
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        gffread -x {output} -g {input.genome} {input.orfanage_cds}
        """

rule orf_prediction_correct_stop_codon_orfanage:
    input:
        config['lr']['orfanage']['cds_filtered'],
    output:
        config['lr']['orfanage']['stop_codon_orf'],
    params:
        scripts_dir = p
    resources:
        threads = 1,
        nodes = 1
    shell:
        """python {params.scripts_dir}/correct_stop_codon_orfanage.py \
            --orfanage_gtf_file_path {input} \
            --output_path {output}"""


rule orf_prediction_extract_sequence_for_cpat:
    input:
        genome=config['ref']['fa'],
        query=config['lr']['orfanage']['cds_to_be_predicted'],
    output:
        config['lr']['orfanage']['missing_cds'],
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        gffread -w {output} -g {input.genome} {input.query}
        """


rule orf_prediction_run_cpat_human:
    input:
        hexamer=config['ref']['cpat']['hexamer'],
        logit_model=config['ref']['cpat']['model'],
        query=config['lr']['orfanage']['missing_cds'],
    output:
        config['lr']['cpat']['r'],
        config['lr']['cpat']['orf_seqs'],
        config['lr']['cpat']['prob_tsv'],
        config['lr']['cpat']['prob_best'],
        config['lr']['cpat']['log'],
        config['lr']['cpat']['err'],
        config['lr']['cpat']['no_orf']
    conda:
        'base'
    resources:
        threads = 1,
        nodes = 1
    params:
        min_orf=config['params']['orfanage']['min_orf_len_nt'],
        top_orf=config['params']['cpat']['top_orf'],
    shell:
        """
        module load gcc
        module load R/4.3.0
        cpat \
                -x {input.hexamer} \
                -d {input.logit_model} \
                -g {input.query} \
                --min-orf={params.min_orf} \
                --top-orf={params.top_orf} \
                -o ../../data/protein/transcripts \
                1> ../../data/protein/transcripts_cpat.output \
                2>../../data/protein/transcripts_cpat.error"""

rule postprocess_check_orf_completeness:
    input:
        cpat_seqs=config['lr']['cpat']['orf_seqs'],
        orfanage_seqs=config['lr']['orfanage']['orfs'],
        cpat_info=config['lr']['cpat']['orf_remaining'],
        orfanage_info=config['lr']['orfanage']['stop_codon_orf'],
    output:
        config['lr']['cpat']['orf_completeness'],
    resources:
        threads = 1,
        nodes = 1
    params:
        scripts_dir = p
    shell:
        """python {params.scripts_dir}/check_orf_completeness.py \
            --cpat_seqs {input.cpat_seqs} \
            --orfanage_seqs {input.orfanage_seqs} \
            --cpat_info {input.cpat_info} \
            --orfanage_info {input.orfanage_info} \
            --output_path {output}"""

rule orf_prediction_filter_cpat_human:
    input:
        input_file_path=config['lr']['cpat']['prob_tsv'],
        orf_input_seq_path=config['lr']['cpat']['orf_seqs'],
    output:
        config['lr']['cpat']['orf_remaining'],
    params:
        first_cutoff=config['params']['cpat']['cutoff_1'],
        second_cutoff=config['params']['cpat']['cutoff_2'],
        scripts_dir = p
    resources:
        threads = 1,
        nodes = 1
    shell:
        """python {params.scripts_dir}/filter_cpat.py \
            --input_file_path {input.input_file_path} \
            --orf_input_seq_path {input.orf_input_seq_path} \
            --output_path {output} \
            --first_cutoff {params.first_cutoff} \
            --second_cutoff {params.second_cutoff}"""


# rule postprocess_check_orf_completeness:
#     input:
#         cpat_seqs=config['lr']['cpat']['orf_seqs'],
#         orfanage_seqs=config['lr']['orfanage']['orfs'],
#         cpat_info=config['lr']['cpat']['orf_remaining'],
#         orfanage_info=config['lr']['orfanage']['stop_codon_orf'],
#     output:
#         config['lr']['cpat']['orf_completeness'],
#     resources:
#         threads = 1,
#         nodes = 1
#     shell:
#         """python {params.scripts_dir}/check_orf_completeness.py \
#             --cpat_seqs {input.cpat_seqs} \
#             --orfanage_seqs {input.orfanage_seqs} \
#             --cpat_info {input.cpat_info} \
#             --orfanage_info {input.orfanage_info} \
#             --output_path {output}"""


rule postprocess_create_cpat_cds_coordinates:
    input:
        sample_gtf=config['lr']['orfanage']['cds_to_be_predicted'],
        called_orfs=config['lr']['cpat']['orf_remaining'],
    output:
        config['lr']['cpat']['cds_coords']
    params:
        opref = get_odir_and_pref_from_fname(config['lr']['cpat']['cds_coords']),
        scripts_dir = p
    resources:
        threads = 1,
        nodes = 1
    shell:
        """python {params.scripts_dir}/create_cpat_CDS_coordinates.py \
                                --name {params.opref} \
                                --sample_gtf {input.sample_gtf} \
                                --called_orfs {input.called_orfs}"""


rule postprocess_combine_cds_gtf:
    input:
        cpat_cds=config['lr']['cpat']['cds_coords'],
        orfanage_cds=config['lr']['orfanage']['stop_codon_orf'],
    output:
        config['lr']['cpat']['unsourced'],
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        cat {input.cpat_cds} {input.orfanage_cds} | gffread -T - | sort -k1,1V -k4,4n -k5,5rn -k3,3r | gffread -T - > {output}
        """


rule postprocess_amend_cds_source:
    input:
        cds=config['lr']['cpat']['unsourced'],
        source_gtf=config['lr']['gtf_no_spike'],
        cpat_cds=config['lr']['cpat']['cds_coords'],
        orfanage_cds=config['lr']['orfanage']['stop_codon_orf'],
    output:
        config['lr']['cpat']['protein']
    resources:
        threads = 1,
        nodes = 1
    params:
        scripts_dir = p
    shell:
        """python {params.scripts_dir}/recover_source.py \
                                --combined_cds_path {input.cds} \
                                --source_gtf_path {input.source_gtf} \
                                --cpat_cds_path {input.cpat_cds} \
                                --orfanage_cds_path {input.orfanage_cds} \
                                --output_path {output}"""


rule postprocess_extract_orf_fasta:
    input:
        protein_gtf=config['lr']['cpat']['protein'],
        genome=config['ref']['fa'],
    output:
        config['lr']['cpat']['orf'],
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        gffread -x {output} -g {input.genome} {input.protein_gtf}
        """


rule postprocess_extract_protein_fasta:
    input:
        protein_gtf=config['lr']['cpat']['protein'],
        genome=config['ref']['fa'],
    output:
        config['lr']['cpat']['protein_fa'],
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        conda activate /gpfs/projects/bsc83/utils/conda_envs/gffread
        gffread -y {output} -g {input.genome} {input.protein_gtf}
        """


rule postprocess_prepare_sqanti_protein_gtf:
    input:
        protein_gtf=config['lr']['cpat']['protein'],
        annotation=config['ref']['gtf'],
    output:
        config['lr']['sqanti_protein']['ref_exons'],
        config['lr']['sqanti_protein']['ref_cds_renamed'],
        config['lr']['sqanti_protein']['data_exons'],
        config['lr']['sqanti_protein']['data_cds_renamed'],
    params:
        opref = get_odir_and_pref_from_fname(config['lr']['cpat']['protein']),
        scripts_dir = p
    resources:
        threads = 1,
        nodes = 3
    shell:
        """python {params.scripts_dir}/rename_cds_to_exon.py \
            --sample_gtf {input.protein_gtf} \
            --sample_name {params.opref} \
            --reference_gtf {input.annotation} \
            --num_cores {resources.threads}"""


rule postprocess_prepare_sqanti_protein_tsv:
    input:
        transcript_only_exons=config['lr']['sqanti_protein']['data_exons'],
        cds_renamed=config['lr']['sqanti_protein']['data_cds_renamed'],
    output:
        config['lr']['sqanti_protein']['best_orf']
    resources:
        threads = 1,
        nodes = 1
    params:
        scripts_dir=p
    shell:
        """python {params.scripts_dir}/create_orf_table_for_sqanti_protein.py \
                    --transcript_exons_path {input.transcript_only_exons} \
                    --cds_only_path {input.cds_renamed} \
                    --output_prefix ../../data/protein/transcripts"""


rule postprocess_run_sqanti_protein:
    input:
        best_orfs=config['lr']['sqanti_protein']['best_orf'],
        renamed_exons=config['lr']['sqanti_protein']['data_exons'],
        cds_only=config['lr']['sqanti_protein']['data_cds_renamed'],
        gencode_renamed_exons=config['lr']['sqanti_protein']['ref_exons'],
        gencode_cds_only=config['lr']['sqanti_protein']['ref_cds_renamed'],
    output:
        config['lr']['sqanti_protein']['classified'],
    params:
        opref = get_odir_and_pref_from_fname(config['lr']['sqanti_protein']['classified']),
        scripts_dir=p
    resources:
        threads = 1,
        nodes = 1
    conda:
        'ucsctools'
    shell:
        """
        python {params.scripts_dir}/sqanti3_protein.py \
                    {input.renamed_exons} \
                    {input.cds_only} \
                    {input.best_orfs} \
                    {input.gencode_renamed_exons} \
                    {input.gencode_cds_only} \
                    -d ./ \
                    -p {params.opref}"""


rule postprocess_summarize_all:
    input:
        best_orf=config['lr']['sqanti_protein']['best_orf'],
        protein_classification=config['lr']['sqanti_protein']['classified'],
        orf_completeness=config['lr']['cpat']['orf_completeness'],
        original_gtf=config['lr']['gtf_no_spike'],
        gtf_predicted=config['lr']['cpat']['protein'],
        protein_fasta=config['lr']['cpat']['protein_fa'],
        blastp=config['lr']['blast']['out'],
    output:
        config['lr']['protein']['summary'],
    params:
        scripts_dir=p
    resources:
        threads = 1,
        nodes = 1
    shell:
        """python {params.scripts_dir}/create_protein_overview_table.py \
            --best_orf_path {input.best_orf} \
            --sqanti_protein_path {input.protein_classification} \
            --orf_completeness_path {input.orf_completeness} \
            --output_name {output} \
            --gtf_original_path {input.original_gtf} \
            --gtf_predicted_path {input.gtf_predicted} \
            --protein_fasta_path {input.protein_fasta} \
            --blastp_path {input.blastp}
            """


rule postprocess_prepare_protein_fasta_for_blast:
    input:
        protein_fasta=config['ref']['pc'],
    output:
        config['ref']['pc_renamed']
    resources:
        threads = 1,
        nodes = 1
    shell:
        "sed -r 's/\|[^\|]*//2g' {input.protein_fasta} > {output}"


rule postprocess_create_blast_db:
    input:
        protein_fasta=config['ref']['pc_renamed'],
    output:
        config['lr']['blast']['blast_db'],
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load blast
        makeblastdb \
                -in {input.protein_fasta} \
                -dbtype prot \
                -parse_seqids \
                -out ../../data/protein/gencode.pc_translations_renamed
        """


rule postprocess_run_blast:
    input:
        protein_fasta=config['lr']['cpat']['protein_fa'],
        protein_reference=config['ref']['pc_renamed'],
        dbs=config['lr']['blast']['blast_db'],
    output:
        config['lr']['blast']['out'],
    resources:
        threads = 8,
        nodes = 2
    params:
        blast_evalue=config['params']['blast']['eval'],
    shell:
        """
        module load blast
        blastp \
            -evalue {params.blast_evalue} \
            -num_threads {resources.threads} \
            -outfmt 6 \
            -db ../../data/protein/gencode.pc_translations_renamed \
            -query {input.protein_fasta} > \
            {output}
        """
