"""
BVBRC Service MCP Prompts

This module contains all the MCP prompt functions for the BVBRC Service MCP Server.
All prompts are registered with the FastMCP server instance.
"""

from fastmcp import FastMCP
from typing import Any, Dict, Optional
from datetime import datetime


def register_prompts(mcp: FastMCP) -> None:
    """
    Register all MCP prompts with the FastMCP server instance.
    
    Args:
        mcp: FastMCP server instance
    """
    
    @mcp.tool(name="date_app_parameters", description="Get the parameters for the date app.")
    def prompt_date() -> str:
        date_prompt = """
        The parameters for the date app are:
        - output_path: The path for the output files (required)
        - output_file: The name of the output file (required)
        """
        return date_prompt

    @mcp.tool(name="genome_assembly_parameters", description="Get the parameters for the genome assembly app.")
    def prompt_genome_assembly() -> str:
        genome_assembly_prompt = """
        The parameters for the genome assembly app are:
          - paired_end_libs: Paired-end read libraries (optional)
          - single_end_libs: Single-end read libraries (optional) 
          - srr_ids: SRA Run IDs (optional)
          - max_bases: Maximum bases before downsampling (default: 10B)
          - recipe: Assembly algorithm (default: 'auto')
            - options: auto, unicycler, flye, meta-flye, canu, spades, meta-spades, plasmid-spades, single-cell, megahit
          - racon_iter: Racon polishing iterations (default: 2)
          - pilon_iter: Pilon polishing iterations (default: 2)
          - trim: Trim reads before assembly (default: False)
          - target_depth: Target depth for normalization (default: 200)
          - normalize: Normalize reads with BBNorm (default: False)
          - filtlong: Filter long reads (default: False)
          - genome_size: Estimated genome size (default: 5M)
          - min_contig_len: Minimum contig length (default: 300)
          - min_contig_cov: Minimum contig coverage (default: 5.0)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        Either paired_end_libs or single_end_libs or srr_ids must be provided.
        """
        return genome_assembly_prompt

    @mcp.tool(name="genome_annotation_parameters", description="Get the parameters for the genome annotation app.")
    def prompt_genome_annotation() -> str:
        genome_annotation_prompt = """
          - contigs: Input contig file (required)
          - scientific_name: Scientific name of organism (required)
          - taxonomy_id: NCBI Taxonomy ID (optional)
          - public: Make genome public (default: False)
          - output_path: Output directory (required)
          - output_file: Output basename (required)
          - reference_genome_id: Reference genome ID (optional)
          - lowvan_min_contig_length: Minimum contig length for lowvan (default: 300)
          - lowvan_max_contig_length: Maximum contig length for lowvan (default: 35000)
          - reference_virus_name: Reference virus name (optional)
          - fix_errors: Automatically try to fix annotation errors: can involve removing candidate genes (optional)
          - fix_frameshifts: Fix frameshifts (optional)
          - analyze_quality: Enable quality analysis (optional)
          - assembly_output: Workspace path for assembly output (optional)
        """
        return genome_annotation_prompt

    @mcp.tool(name="comprehensive_genome_analysis_parameters", description="Get the parameters for the comprehensive genome analysis app.")
    def prompt_comprehensive_genome_analysis() -> str:
        comprehensive_genome_analysis_prompt = """
        The parameters for the comprehensive genome analysis app are:
          - input_type: Input type - reads/contigs/genbank (required)
          - output_path: Output directory (required)
          - output_file: Output basename (required)
          - paired_end_libs: Paired-end read libraries (optional)
          - single_end_libs: Single-end read libraries (optional)
          - srr_ids: SRA Run IDs (optional)
          - reference_assembly: Reference contig file (optional)
          - recipe: Assembly algorithm (default: 'auto')
            - options: auto, unicycler, flye, meta-flye, canu, spades, meta-spades, plasmid-spades, single-cell, megahit
          - racon_iter: Racon polishing iterations (default: 2)
          - pilon_iter: Pilon polishing iterations (default: 2)
          - trim: Trim reads before assembly (default: False)
          - normalize: Normalize reads with BBNorm (default: False)
          - filtlong: Filter long reads (default: False)
          - target_depth: Target depth for normalization (default: 200)
          - genome_size: Estimated genome size (default: 5M)
          - min_contig_len: Minimum contig length (default: 300)
          - min_contig_cov: Minimum contig coverage (default: 5.0)
          - gto: Preannotated genome object (optional)
          - genbank_file: GenBank file (optional)
          - contigs: Contig file for annotation (optional)
          - scientific_name: Scientific name of organism (required)
          - taxonomy_id: NCBI Taxonomy ID (optional)
          - code: Genetic code for translation (required)
          - domain: Organism domain - Bacteria/Archaea/Viruses/auto (required)
          - public: Make genome public (default: False)
          - queue_nowait: Don't wait for indexing (default: False)
          - skip_indexing: Don't index genome (default: False)
          - reference_genome_id: Reference genome ID (optional)
          - analyze_quality: Enable quality analysis (optional)

        At least one input source (paired_end_libs, single_end_libs, srr_ids, contigs, genbank_file, or gto) must be provided.
        """
        return comprehensive_genome_analysis_prompt

    @mcp.tool(name="blast_app_parameters", description="Get the parameters for the blast app.")
    def prompt_blast_app() -> str:
        blast_app_prompt = """
        The parameters for the blast app are:
          - input_type: Input sequence type - dna/aa (required)
          - input_source: Input source type (required)
          - input_fasta_data: Input sequences in FASTA format (optional)
          - input_id_list: Input sequence identifiers (optional)
          - input_fasta_file: Input FASTA file (optional)
          - input_feature_group: Input feature group (optional)
          - input_genome_group: Input genome group (optional)
          - db_type: Database type - faa/ffn/frn/fna (required)
          - db_source: Database source type (required)
          - db_fasta_data: Database sequences in FASTA format (optional)
          - db_fasta_file: Database FASTA file (optional)
          - db_id_list: Database sequence identifiers (optional)
          - db_feature_group: Database feature group (optional)
          - db_genome_group: Database genome group (optional)
          - db_genome_list: Database genome list (optional)
          - db_taxon_list: Database taxon list (optional)
          - db_precomputed_database: Precomputed database name (optional)
          - blast_program: BLAST program to use (optional)
          - blast_evalue_cutoff: E-value cutoff (default: 1e-5)
          - blast_max_hits: Maximum hits to return (default: 300)
          - blast_min_coverage: Minimum coverage percentage (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one input source (input_fasta_data, input_id_list, input_fasta_file, input_feature_group, or input_genome_group) must be provided.
        At least one database source (db_fasta_data, db_fasta_file, db_id_list, db_feature_group, db_genome_group, db_genome_list, db_taxon_list, or db_precomputed_database) must be provided.
        """
        return blast_app_prompt

    @mcp.tool(name="primer_design_app_parameters", description="Get the parameters for the primer design app.")
    def prompt_primer_design_app() -> str:
        primer_design_app_prompt = """
        The parameters for the primer design app are:
          - output_file: Output basename (required)
          - output_path: Output directory (required)
          - input_type: Input format - sequence_text/workplace_fasta/database_id (required)
          - sequence_input: DNA sequence data (required)
          - SEQUENCE_ID: Sequence identifier (optional)
          - SEQUENCE_TARGET: Amplification target region (optional)
          - SEQUENCE_INCLUDED_REGION: Region where primers can be picked (optional)
          - SEQUENCE_EXCLUDED_REGION: Region where primers cannot overlap (optional)
          - SEQUENCE_OVERLAP_JUNCTION_LIST: Junction regions primers must flank (optional)
          - PRIMER_PICK_INTERNAL_OLIGO: Pick internal oligo/probe (optional)
          - PRIMER_PRODUCT_SIZE_RANGE: Min/max product size (optional)
          - PRIMER_NUM_RETURN: Max number of primer pairs to report (optional)
          - PRIMER_MIN_SIZE: Minimum primer length (optional)
          - PRIMER_OPT_SIZE: Optimal primer length (optional)
          - PRIMER_MAX_SIZE: Maximum primer length (optional)
          - PRIMER_MAX_TM: Maximum primer melting temperature (optional)
          - PRIMER_MIN_TM: Minimum primer melting temperature (optional)
          - PRIMER_OPT_TM: Optimal primer melting temperature (optional)
          - PRIMER_PAIR_MAX_DIFF_TM: Max Tm difference of paired primers (optional)
          - PRIMER_MAX_GC: Maximum primer GC percentage (optional)
          - PRIMER_MIN_GC: Minimum primer GC percentage (optional)
          - PRIMER_OPT_GC: Optimal primer GC percentage (optional)
          - PRIMER_SALT_MONOVALENT: Monovalent cation concentration (mM) (optional)
          - PRIMER_SALT_DIVALENT: Divalent cation concentration (mM) (optional)
          - PRIMER_DNA_CONC: Annealing oligo concentration (nM) (optional)
          - PRIMER_DNTP_CONC: dNTP concentration (mM) (optional)
        """
        return primer_design_app_prompt

    @mcp.tool(name="variation_app_parameters", description="Get the parameters for the variation app.")
    def prompt_variation_app() -> str:
        variation_app_prompt = """
        The parameters for the variation app are:
          - reference_genome_id: Reference genome ID (required)
          - paired_end_libs: Paired-end read libraries (optional)
          - single_end_libs: Single-end read libraries (optional)
          - srr_ids: SRA Run IDs (optional)
          - mapper: Short read mapper - BWA-mem/BWA-mem-strict/Bowtie2/MOSAIK/LAST/minimap2 (default: BWA-mem)
          - caller: SNP caller - FreeBayes/BCFtools (default: FreeBayes)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one input source (paired_end_libs, single_end_libs, or srr_ids) must be provided.
        """
        return variation_app_prompt

    @mcp.tool(name="tnseq_app_parameters", description="Get the parameters for the tnseq app.")
    def prompt_tnseq_app() -> str:
        tnseq_app_prompt = """
        The parameters for the tnseq app are:
          - experimental_conditions: Experimental conditions for analysis (optional)
          - contrasts: Statistical contrasts for comparison (optional)
          - read_files: Read file groups (optional)
          - reference_genome_id: Reference genome ID (required)
          - recipe: TnSeq analysis recipe - gumbel/griffin/tn5gaps/rankproduct/hmm/binomial/resampling (default: gumbel)
          - protocol: TnSeq protocol - sassetti/tn5/mme1 (default: sassetti)
          - primer: Primer DNA string for read trimming (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)
        """
        return tnseq_app_prompt

    @mcp.tool(name="bacterial_genome_tree_app_parameters", description="Get the parameters for the bacterial genome tree app.")
    def prompt_bacterial_genome_tree_app() -> str:
        bacterial_genome_tree_app_prompt = """
        The parameters for the bacterial genome tree app are:
          - output_path: Output directory (required)
          - output_file: Output basename (required)
          - genome_ids: Main genomes for analysis (optional)
          - genome_groups: Main genome groups (optional)
          - optional_genome_ids: Optional genomes (not penalized for missing genes) (optional)
          - genome_metadata_fields: Genome metadata fields to retrieve (optional)
          - number_of_genes: Desired number of genes for analysis (default: 20)
          - bootstraps: Number of bootstrap replicates (default: 100)
          - max_genomes_missing: Max main genomes missing from any PGFam (default: 0)
          - max_allowed_dups: Max duplications allowed for main genomes (default: 0)

        At least one genome input (genome_ids or genome_groups) must be provided.
        """
        return bacterial_genome_tree_app_prompt

    @mcp.tool(name="gene_tree_app_parameters", description="Get the parameters for the gene tree app.")
    def prompt_gene_tree_app() -> str:
        gene_tree_app_prompt = """
        The parameters for the gene tree app are:
          - sequences: Sequence data inputs (required)
          - alignment_program: Alignment program - muscle/mafft (optional)
          - trim_threshold: Alignment end-trimming threshold (optional)
          - gap_threshold: Delete gappy sequences threshold (optional)
          - alphabet: Sequence alphabet - DNA/Protein (required)
          - substitution_model: Substitution model (optional)
          - bootstrap: Number of bootstrap replicates (optional)
          - recipe: Tree construction recipe - RAxML/PhyML/FastTree (default: RAxML)
          - tree_type: Type of tree - viral_genome/gene (optional)
          - feature_metadata_fields: Gene metadata fields (optional)
          - genome_metadata_fields: Genome metadata fields (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)
        """
        return gene_tree_app_prompt

    @mcp.tool(name="core_genome_mlst_app_parameters", description="Get the parameters for the core genome MLST app.")
    def prompt_core_genome_mlst_app() -> str:
        core_genome_mlst_app_prompt = """
        The parameters for the core genome MLST app are:
          - input_genome_type: Input type - genome_group/genome_fasta (required, default: genome_group)
          - analysis_type: Analysis type - chewbbaca (required, default: chewbbaca)
          - input_genome_group: Genome group name (optional)
          - input_genome_fasta: Genome FASTA data (optional)
          - schema_location: Path to schema directory (optional)
          - input_schema_selection: Schema species selection (required)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one genome input (input_genome_group or input_genome_fasta) must be provided.
        """
        return core_genome_mlst_app_prompt

    @mcp.tool(name="whole_genome_snp_app_parameters", description="Get the parameters for the whole genome SNP app.")
    def prompt_whole_genome_snp_app() -> str:
        whole_genome_snp_app_prompt = """
        The parameters for the whole genome SNP app are:
          - input_genome_type: Input type - genome_group/genome_fasta (required)
          - majority_threshold: Majority fraction for SNP inclusion (default: 0.5)
          - min_mid_linkage: Minimum value for mid linkage (default: 10)
          - max_mid_linkage: Maximum value for mid linkage (default: 40)
          - analysis_type: Analysis type - Whole Genome SNP Analysis (required)
          - input_genome_group: Genome group name (optional)
          - input_genome_fasta: Genome FASTA data (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one genome input (input_genome_group or input_genome_fasta) must be provided.
        """
        return whole_genome_snp_app_prompt

    @mcp.tool(name="taxonomic_classification_app_parameters", description="Get the parameters for the taxonomic classification app.")
    def prompt_taxonomic_classification_app() -> str:
        taxonomic_classification_app_prompt = """
        The parameters for the taxonomic classification app are:
          - host_genome: Host genome for filtering (required, default: no_host)
          - analysis_type: Analysis type - pathogen/microbiome/16S (required, default: 16S)
          - paired_end_libs: Paired-end read libraries with sample IDs (optional)
          - single_end_libs: Single-end read libraries with sample IDs (optional)
          - srr_libs: SRA datasets with sample IDs (optional)
          - database: Reference database - bvbrc/greengenes/silva/standard (required, default: SILVA)
          - save_classified_sequences: Save classified sequences (default: False)
          - save_unclassified_sequences: Save unclassified sequences (default: False)
          - confidence_interval: Classification confidence threshold (default: 0.1)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one input source (paired_end_libs, single_end_libs, or srr_libs) must be provided.
        """
        return taxonomic_classification_app_prompt

    @mcp.tool(name="metagenomic_binning_app_parameters", description="Get the parameters for the metagenomic binning app.")
    def prompt_metagenomic_binning_app() -> str:
        metagenomic_binning_app_prompt = """
        The parameters for the metagenomic binning app are:
          - paired_end_libs: Paired-end read libraries (optional)
          - single_end_libs: Single-end read libraries (optional)
          - srr_ids: SRA Run IDs (optional)
          - contigs: Input contig file (optional)
          - genome_group: Output genome group name (optional)
          - skip_indexing: Don't index bins in Solr (default: False)
          - recipe: Custom annotation recipe (optional)
          - viral_recipe: Custom viral annotation recipe (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)
          - force_local_assembly: Force local assembly (required, default: False)
          - force_inline_annotation: Force inline annotation (default: True)
          - perform_bacterial_binning: Perform bacterial binning (default: True)
          - perform_viral_binning: Perform viral binning (default: False)
          - perform_viral_annotation: Perform viral annotation (default: False)
          - perform_bacterial_annotation: Perform bacterial annotation (default: True)
          - assembler: Assembler to use (optional)
          - danglen: Dangling element length (default: 50)
          - min_contig_len: Minimum contig length (default: 400)
          - min_contig_cov: Minimum contig coverage (default: 4.0)

        At least one input source (paired_end_libs, single_end_libs, srr_ids, contigs, or genome_group) must be provided.
        """
        return metagenomic_binning_app_prompt

    @mcp.tool(name="metagenomic_read_mapping_app_parameters", description="Get the parameters for the metagenomic read mapping app.")
    def prompt_metagenomic_read_mapping_app() -> str:
        metagenomic_read_mapping_app_prompt = """
        The parameters for the metagenomic read mapping app are:
          - gene_set_type: Input type - predefined_list/fasta_file/feature_group (required)
          - gene_set_name: Predefined gene set - VFDB/CARD/feature_group/fasta_file (optional)
          - gene_set_fasta: Gene set FASTA data (optional)
          - gene_set_feature_group: Gene set feature group name (optional)
          - paired_end_libs: Paired-end read libraries (optional)
          - single_end_libs: Single-end read libraries (optional)
          - srr_ids: SRA Run IDs (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one gene set source (gene_set_name, gene_set_fasta, or gene_set_feature_group) must be provided.
        At least one input source (paired_end_libs, single_end_libs, or srr_ids) must be provided.
        """
        return metagenomic_read_mapping_app_prompt

    @mcp.tool(name="rnaseq_app_parameters", description="Get the parameters for the RNA-seq app.")
    def prompt_rnaseq_app() -> str:
        rnaseq_app_prompt = """
        The parameters for the RNA-seq app are:
          - experimental_conditions: Experimental conditions for analysis (optional)
          - contrasts: Statistical contrasts for comparison (optional)
          - strand_specific: Strand-specific reads (default: True)
          - paired_end_libs: Paired-end read libraries with sample IDs (optional)
          - single_end_libs: Single-end read libraries with sample IDs (optional)
          - srr_libs: SRA datasets with sample IDs (optional)
          - reference_genome_id: Reference genome ID (required)
          - genome_type: Genome type - bacteria/host (required)
          - recipe: Analysis recipe - HTSeq-DESeq/cufflinks/Host (required, default: HTSeq-DESeq)
          - host_ftp: Host FTP prefix for file access (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)
          - trimming: Run TrimGalore on reads (default: False)
          - unit_test: Unit test JSON file path (optional)
          - skip_sampling: Skip sampling step in alignment (optional)

        At least one input source (paired_end_libs, single_end_libs, or srr_libs) must be provided.
        """
        return rnaseq_app_prompt

    @mcp.tool(name="expression_import_app_parameters", description="Get the parameters for the expression import app.")
    def prompt_expression_import_app() -> str:
        expression_import_app_prompt = """
        The parameters for the expression import app are:
          - xfile: Experiment data file with comparison values (required)
          - mfile: Metadata template file (optional)
          - ustring: User information JSON string (required)
          - output_path: Output directory (optional)
          - output_file: Output basename (optional)
        """
        return expression_import_app_prompt

    @mcp.tool(name="sars_wastewater_analysis_app_parameters", description="Get the parameters for the SARS wastewater analysis app.")
    def prompt_sars_wastewater_analysis_app() -> str:
        sars_wastewater_analysis_app_prompt = """
        The parameters for the SARS wastewater analysis app are:
          - paired_end_libs: Paired-end read libraries with sample metadata (optional)
          - single_end_libs: Single-end read libraries with sample metadata (optional)
          - srr_libs: SRA datasets with sample metadata (optional)
          - recipe: Assembly recipe - onecodex (default: auto)
          - primers: Primer set - ARTIC/midnight/qiagen/swift/varskip/varskip-long (required, default: ARTIC)
          - minimum_base_quality_score: Minimum base quality score (default: 20)
          - minimum_genome_coverage: Minimum genome coverage (default: 60)
          - agg_minimum_lineage_abundance: Minimum lineage abundance for plots (default: 0.01)
          - minimum_coverage_depth: Minimum coverage depth (default: 0)
          - confirmedonly: Exclude unconfirmed lineages (default: False)
          - minimum_lineage_abundance: Minimum lineage abundance (default: 0.001)
          - coverage_estimate: Coverage estimate value (default: 10)
          - timeseries_plot_interval: Timeseries plot interval (default: 0)
          - primer_version: Primer version number (optional)
          - barcode_csv: Custom barcode path (optional)
          - sample_metadata_csv: Sample metadata CSV (optional)
          - keep_intermediates: Keep intermediate outputs (default: True)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one input source (paired_end_libs, single_end_libs, or srr_libs) must be provided.
        """
        return sars_wastewater_analysis_app_prompt

    @mcp.tool(name="sequence_submission_app_parameters", description="Get the parameters for the sequence submission app.")
    def prompt_sequence_submission_app() -> str:
        sequence_submission_app_prompt = """
        The parameters for the sequence submission app are:
          - input_source: Input source type - id_list/fasta_data/fasta_file/genome_group (required)
          - input_fasta_data: Input sequences in FASTA format (optional)
          - input_fasta_file: Input FASTA file (optional)
          - input_genome_group: Input genome group (optional)
          - metadata: Metadata CSV file (required)
          - affiliation: Submitter affiliation (optional)
          - first_name: Submitter first name (required)
          - last_name: Submitter last name (required)
          - email: Submitter email (required)
          - consortium: Submitter consortium (optional)
          - country: Submitter country (optional)
          - phoneNumber: Submitter phone number (optional)
          - street: Submitter street address (optional)
          - postal_code: Submitter postal code (optional)
          - city: Submitter city (optional)
          - state: Submitter state (optional)
          - numberOfSequences: Number of sequences in submission (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one sequence input (input_fasta_data, input_fasta_file, or input_genome_group) must be provided.
        """
        return sequence_submission_app_prompt

    @mcp.tool(name="influenza_ha_subtype_conversion_app_parameters", description="Get the parameters for the influenza HA subtype conversion app.")
    def prompt_influenza_ha_subtype_conversion_app() -> str:
        influenza_ha_subtype_conversion_app_prompt = """
        The parameters for the influenza HA subtype conversion app are:
          - input_source: Input source type - feature_list/fasta_data/fasta_file/feature_group (required)
          - input_fasta_data: Input sequences in FASTA format (optional)
          - input_fasta_file: Input FASTA file (optional)
          - input_feature_group: Input feature group (optional)
          - input_feature_list: Input feature ID list (optional)
          - types: Selected types in the submission (required)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one sequence input (input_fasta_data, input_fasta_file, input_feature_group, or input_feature_list) must be provided.
        """
        return influenza_ha_subtype_conversion_app_prompt

    @mcp.tool(name="subspecies_classification_app_parameters", description="Get the parameters for the subspecies classification app.")
    def prompt_subspecies_classification_app() -> str:
        subspecies_classification_app_prompt = """
        The parameters for the subspecies classification app are:
          - input_source: Input source type - id_list/fasta_data/fasta_file/genome_group (required)
          - input_fasta_data: Input sequences in FASTA format (optional)
          - input_fasta_file: Input FASTA file (optional)
          - input_genome_group: Input genome group (optional)
          - ref_msa_fasta: Reference multiple sequence alignment (optional)
          - virus_type: Virus type for classification (required)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one sequence input (input_fasta_data, input_fasta_file, or input_genome_group) must be provided.
        """
        return subspecies_classification_app_prompt

    @mcp.tool(name="viral_assembly_app_parameters", description="Get the parameters for the viral assembly app.")
    def prompt_viral_assembly_app() -> str:
        viral_assembly_app_prompt = """
        The parameters for the viral assembly app are:
          - paired_end_lib: Paired-end read library (optional)
          - single_end_lib: Single-end read library (optional)
          - srr_id: SRA Run ID (optional)
          - recipe: Assembly strategy - auto/irma (default: auto)
          - module: Virus module - FLU/CoV/RSV/EBOLA/FLU_AD (optional)
          - viral_size: Estimated viral genome size (default: 5M)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one input source (paired_end_lib, single_end_lib, or srr_id) must be provided.
        """
        return viral_assembly_app_prompt

    @mcp.tool(name="genome_alignment_app_parameters", description="Get the parameters for the genome alignment app.")
    def prompt_genome_alignment_app() -> str:
        genome_alignment_app_prompt = """
        The parameters for the genome alignment app are:
          - genome_ids: Genome IDs for alignment (required)
          - recipe: Alignment algorithm - progressiveMauve (default: progressiveMauve)
          - seedWeight: Seed weight for alignment (optional)
          - maxGappedAlignerLength: Maximum gapped aligner length (optional)
          - maxBreakpointDistanceScale: Maximum breakpoint distance scale (optional)
          - conservationDistanceScale: Conservation distance scale (optional)
          - weight: Alignment weight (optional)
          - minScaledPenalty: Minimum scaled penalty (optional)
          - hmmPGoHomologous: HMM probability for homologous regions (optional)
          - hmmPGoUnrelated: HMM probability for unrelated regions (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)
        """
        return genome_alignment_app_prompt

    @mcp.tool(name="sars_genome_analysis_app_parameters", description="Get the parameters for the SARS genome analysis app.")
    def prompt_sars_genome_analysis_app() -> str:
        sars_genome_analysis_app_prompt = """
        The parameters for the SARS genome analysis app are:
          - paired_end_libs: Paired-end read libraries (optional)
          - single_end_libs: Single-end read libraries (optional)
          - srr_ids: SRA Run IDs (optional)
          - recipe: Assembly recipe - auto (default: auto)
          - primers: Primer set - ARTIC/midnight/qiagen/swift/varskip/varskip-long (required, default: ARTIC)
          - primer_version: Primer version number (optional)
          - min_depth: Minimum depth for consensus (default: 100)
          - max_depth: Maximum depth for consensus (default: 8000)
          - keep_intermediates: Keep intermediate outputs (default: 0)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one input source (paired_end_libs, single_end_libs, or srr_ids) must be provided.
        """
        return sars_genome_analysis_app_prompt

    @mcp.tool(name="msa_snp_analysis_app_parameters", description="Get the parameters for the MSA SNP analysis app.")
    def prompt_msa_snp_analysis_app() -> str:
        msa_snp_analysis_app_prompt = """
        The parameters for the MSA SNP analysis app are:
          - input_status: Input alignment status - unaligned/aligned (default: unaligned)
          - input_type: Input type - input_group (default: input_group)
          - fasta_files: FASTA file inputs (optional)
          - select_genomegroup: Genome group selection (optional)
          - feature_groups: Feature group selection (optional)
          - feature_list: Feature list selection (optional)
          - genome_list: Genome list selection (optional)
          - aligner: Alignment program - Muscle/MAFFT/ClustalW (default: Muscle)
          - alphabet: Sequence alphabet - dna/aa (default: dna)
          - fasta_keyboard_input: Direct FASTA input (optional)
          - ref_type: Reference type - none (default: none)
          - strategy: Analysis strategy - auto (default: auto)
          - ref_string: Reference sequence string (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one sequence input (fasta_files, select_genomegroup, feature_groups, feature_list, genome_list, or fasta_keyboard_input) must be provided.
        """
        return msa_snp_analysis_app_prompt

    @mcp.tool(name="metacats_app_parameters", description="Get the parameters for the metacats app.")
    def prompt_metacats_app() -> str:
        metacats_app_prompt = """
        The parameters for the metacats app are:
          - output_path: Output directory (required)
          - output_file: Output basename (required)
          - p_value: P-value threshold for significance (default: 0.05)
          - year_ranges: Year range filters (optional)
          - metadata_group: Metadata group for analysis (optional)
          - input_type: Input data type (optional)
          - alphabet: Sequence alphabet - na (default: na)
          - groups: Analysis groups (optional)
          - alignment_file: Input alignment file (optional)
          - group_file: Group definition file (optional)
          - alignment_type: Alignment type (optional)
          - auto_groups: Automatic group detection (optional)
        """
        return metacats_app_prompt

    @mcp.tool(name="proteome_comparison_app_parameters", description="Get the parameters for the proteome comparison app.")
    def prompt_proteome_comparison_app() -> str:
        proteome_comparison_app_prompt = """
        The parameters for the proteome comparison app are:
          - genome_ids: Reference genome IDs (optional)
          - user_genomes: User-provided genome IDs (optional)
          - user_feature_groups: User-provided feature groups (optional)
          - reference_genome_index: Index of reference genome (default: 1)
          - min_seq_cov: Minimum sequence coverage (default: 0.30)
          - max_e_val: Maximum E-value threshold (default: 1e-5)
          - min_ident: Minimum identity threshold (default: 0.1)
          - min_positives: Minimum positives threshold (default: 0.2)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one genome input (genome_ids, user_genomes, or user_feature_groups) must be provided.
        """
        return proteome_comparison_app_prompt

    @mcp.tool(name="comparative_systems_app_parameters", description="Get the parameters for the comparative systems app.")
    def prompt_comparative_systems_app() -> str:
        comparative_systems_app_prompt = """
        The parameters for the comparative systems app are:
          - output_path: Output directory (required)
          - output_file: Output basename (required)
          - genome_ids: Genome IDs for comparison (optional)
          - genome_groups: Genome groups for comparison (optional)

        At least one genome input (genome_ids or genome_groups) must be provided.
        """
        return comparative_systems_app_prompt

    @mcp.tool(name="docking_app_parameters", description="Get the parameters for the docking app.")
    def prompt_docking_app() -> str:
        docking_app_prompt = """
        The parameters for the docking app are:
          - protein_input_type: Protein input type (required)
          - input_pdb: PDB structure IDs (optional)
          - user_pdb_file: User-provided PDB files (optional)
          - ligand_library_type: Ligand library type (required)
          - ligand_named_library: Named ligand library (optional)
          - ligand_smiles_list: SMILES string list (optional)
          - ligand_ws_file: Workspace ligand file (optional)
          - top_n: Number of top results to return (optional)
          - batch_size: Batch processing size (default: 10)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one protein input (input_pdb or user_pdb_file) must be provided.
        At least one ligand input (ligand_named_library, ligand_smiles_list, or ligand_ws_file) must be provided.
        """
        return docking_app_prompt

    @mcp.tool(name="similar_genome_finder_app_parameters", description="Get the parameters for the similar genome finder app.")
    def prompt_similar_genome_finder_app() -> str:
        similar_genome_finder_app_prompt = """
        The parameters for the similar genome finder app are:
          - selectedGenomeId: Reference genome ID (optional)
          - fasta_file: FASTA file input (optional)
          - max_pvalue: Maximum p-value threshold (optional)
          - max_distance: Maximum distance threshold (optional)
          - max_hits: Maximum number of hits (optional)
          - include_reference: Include reference genomes (optional)
          - include_representative: Include representative genomes (optional)
          - include_bacterial: Include bacterial genomes (optional)
          - include_viral: Include viral genomes (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)

        At least one query input (selectedGenomeId or fasta_file) must be provided.
        """
        return similar_genome_finder_app_prompt

    @mcp.tool(name="fastqutils_app_parameters", description="Get the parameters for the fastqutils app.")
    def prompt_fastqutils_app() -> str:
        fastqutils_app_prompt = """
        The parameters for the fastqutils app are:
          - reference_genome_id: Reference genome ID for QC (optional)
          - paired_end_libs: Paired-end read libraries (optional)
          - single_end_libs: Single-end read libraries (optional)
          - srr_libs: SRA datasets (optional)
          - output_path: Output directory (required)
          - output_file: Output basename (required)
          - recipe: QC recipe options (optional)
        """
        return fastqutils_app_prompt