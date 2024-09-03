rule meryl_ref_kmer_count:
    resources:
        mem_gb = 64,
        threads = 8
    shell:
        """
        meryl count \
            k={params.kmer_len} \
            output {output.db} \
            {input.fa}
        """

rule meryl_ref_greater_than:
    resources:
        mem_gb = 64,
        threads = 8
    shell:
        """
        meryl print greater-than \
            distinct={params.distinct} \
            {input.db} > {output.kmers}
        """

rule winnowmap_pb:
    resources:
        threads = 8,
        mem_gb = 64
    shell:
        """
        winnowmap \
            -W {input.kmers} \
            -ax map-pb \
            {input.fa} \
            {input.fq} > {output.sam}
        """
