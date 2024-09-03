rule dl_encid:
    resources:
        nodes = 1,
        threads = 1
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.{params.file_ext} -O {output.out}"

rule dl_encid_gz:
    resources:
        nodes = 1,
        threads = 1
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.{params.file_ext}.gz -O {output.out}"

rule wget:
    resources:
        nodes = 1,
        threads = 1
    shell:
        """
        wget "{params.link}" -O {output.out}
        """
