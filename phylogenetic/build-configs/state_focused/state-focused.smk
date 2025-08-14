
# Define global parameters
STRAIN_ID = config["strain_id_field"]

"""
This file provides custom rules to substitute the default sampling scheme.
These rules will supercede those included in the 'rules' directory.
"""
rule filter_genomic_state:
    """
    Filtering {params.division} sequences to
      - {params.sequences_per_group} sequence(s) per {params.group_by!s}
      - from {params.min_date} onwards
      - excluding strains in {input.exclude}
      - minimum genome length of {params.min_length}

    """
    input:
        sequences = "data/sequences.fasta",
        metadata = "data/metadata.tsv",
        exclude = resolve_config_path(config["files"]["exclude"])
    output:
        strains = "results/genome/{params.division}_filtered.txt"
    params:
        strain_id=STRAIN_ID,
        division=config["filter"]["division"],
        group_by=config["filter"]["group_by_state"],
        sequences_per_group=config["filter"]["sequences_per_group_state"],
        min_date=config["filter"]["min_date"],
        min_length=config["filter"]["min_length"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --query "division == '{params.division}'" \
            --exclude {input.exclude} \
            --output-strains results/genome/{params.division}_filtered.txt \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length}
        """

rule filter_genomic_continental:
    """
    Filtering USA and Canada sequences to
      - {params.sequences_per_group} sequence(s) per {params.group_by!s}
      - from {params.min_date} onwards
      - excluding strains in {input.exclude}
      - minimum genome length of {params.min_length}
    """
    input:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv",
        exclude = resolve_config_path(config["files"]["exclude"]),
    output:
        strains="results/genome/continental_filtered.txt"
    params:
        strain_id=STRAIN_ID,
        group_by=config["filter"]["group_by_continental"],
        sequences_per_group=config["filter"]["sequences_per_group_continental"],
        min_date=config["filter"]["min_date"],
        min_length=config["filter"]["min_length"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --query "division != 'Washington' & (country == 'USA' | country == 'Canada')" \
            --exclude {input.exclude} \
            --output-strains {output.strains} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length}
        """

rule combine_sequences_genomic:
"""

"""
    input:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv",
        include = resolve_config_path(config["files"]["include"])({"gene": "genome"}),
    output:
        sequences="results/genome/filtered.fasta"
    params:
        strain_id=STRAIN_ID,
        division=config["filter"]["division"],
        min_date=config["filter"]["min_date"],
        min_length=config["filter"]["min_length"]
    shell:
        """
        augur filter \
        --sequences {input.sequences} \
        --metadata-id-columns {params.strain_id} \
        --metadata {input.metadata} \
        --exclude-all \
        --include results/genome/{params.division}_filtered.txt \
            results/genome/continental_filtered.txt \
            {input.include} \
        --output-sequences {output.sequences} \
        """

"""
This file provides custom rules to substitute the default sampling scheme.
These rules will supercede those included in the 'rules' directory
These rules apply only to the genomic data set (see state_focused_sampling_N450.smk for the N450 data set.)
"""
rule filter_N450_state:
    """
    Filtering {params.division} sequences to
      - {params.sequences_per_group} sequence(s) per {params.group_by!s}
      - from {params.min_date} onwards
      - excluding strains in {input.exclude}
      - minimum genome length of {params.min_length}

    """
    input:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv",
        exclude=resolve_config_path(config["files"]["exclude"]),
    output:
        strains="results/genome/{params.division}_filtered.txt"
    params:
        strain_id=STRAIN_ID,
        division=config["filter"]["division"],
        group_by=config["filter_N450"]["group_by_state"],
        sequences_per_group=config["filter_N450"]["sequences_per_group_state"],
        min_date=config["filter_N450"]["min_date"],
        min_length=config["filter_N450"]["min_length"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --query "division == '{params.division}'" \
            --exclude {input.exclude} \
            --output-strains {output.strains} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length}
        """

rule filter_N450_continental:
    """
    Filtering USA and Canada sequences to
      - {params.sequences_per_group} sequence(s) per {params.group_by!s}
      - from {params.min_date} onwards
      - excluding strains in {input.exclude}
      - minimum genome length of {params.min_length}
    """
    input:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv",
        exclude=resolve_config_path(config["files"]["exclude"])
    output:
        strains="results/genome/continental_filtered.txt",
    params:
        strain_id=STRAIN_ID,
        group_by=config["filter_N450"]["group_by_continental"],
        sequences_per_group=config["filter_N450"]["sequences_per_group_continental"],
        min_date=config["filter_N450"]["min_date"],
        min_length=config["filter_N450"]["min_length"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --query "division != 'Washington' & (country == 'USA' | country == 'Canada')" \
            --exclude {input.exclude} \
            --output-strains {output.strains} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length}
        """

rule combine_sequences_N450:
"""

"""
    input:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv",
        include=resolve_config_path(config["files"]["include"])({"gene": "N450"}),
    output:
        sequences="results/N450/aligned.fasta"
    params:
        strain_id=STRAIN_ID,
        division=config["filter"]["division"],
        min_date=config["filter_N450"]["min_date"],
        min_length=config["filter_N450"]["min_length"]
    shell:
        """
        augur filter \
        --sequences {input.sequences} \
        --metadata-id-columns {params.strain_id} \
        --metadata {input.metadata} \
        --exclude-all \
        --include results/N450/{params.division}_filtered.txt \
            results/N450/continental_filtered.txt \
        --output-sequences {output.sequences} \
        """

ruleorder: filter_genomic_state > filter_genomic_continental > combine_sequences_genomic
ruleorder: filter_N450_state > filter_N450_continental > combine_sequences_N450