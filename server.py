from flask import Flask
from flaskmcp import create_app, tool, register_resource, register_prompt
from json_rpc import JsonRpcCaller
from service_functions import (
    enumerate_apps, start_date_app, start_genome_annotation_app, query_tasks,
    start_genome_assembly_app, start_comprehensive_genome_analysis_app, start_blast_app,
    start_primer_design_app, start_variation_app, start_tnseq_app, start_bacterial_genome_tree_app,
    start_gene_tree_app, start_core_genome_mlst_app, start_whole_genome_snp_app,
    start_taxonomic_classification_app, start_metagenomic_binning_app, start_metagenomic_read_mapping_app,
    start_rnaseq_app, start_expression_import_app, start_sars_wastewater_analysis_app,
    start_sequence_submission_app, start_influenza_ha_subtype_conversion_app,
    start_subspecies_classification_app, start_viral_assembly_app, start_fastqutils_app,
    start_genome_alignment_app, start_sars_genome_analysis_app, start_msa_snp_analysis_app,
    start_metacats_app, start_proteome_comparison_app, start_comparative_systems_app,
    start_docking_app, start_similar_genome_finder_app
)
import json
import sys
from typing import Any, List, Dict

with open('config.json', 'r') as f:
    config = json.load(f)

service_api_url = config['service_api_url']
similar_genome_finder_api_url = config.get('similar_genome_finder_api_url', service_api_url)
port = config.get("port", 5001)

api = JsonRpcCaller(service_api_url)
similar_genome_finder_api = JsonRpcCaller(similar_genome_finder_api_url)

app = create_app({'DEBUG': True})

def extract_userid_from_token(token: str = None) -> str:
    """
    Extract user ID from JWT token.
    Returns a default user ID if token is None or invalid.
    """
    if not token:
        return None
    
    try:
        user_id = token.split('|')[0].replace('un=','')
        return user_id
            
    except Exception as e:
        print(f"Error extracting user ID from token: {e}")
        return None

# Optional: Add a health check route
@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    from flask import jsonify
    return jsonify({"status": "healthy", "service": "bvbrc-app-service-mcp"})

@tool(name="service_enumerate_apps", description="Enumerate all apps. Parameters: token (str, optional) - Authentication token for API access")
def service_enumerate_apps(token: str = None) -> List[str]:
    user_id = extract_userid_from_token(token)
    return enumerate_apps(api, token, user_id=user_id)

@tool(name="service_start_date_app", description="Start the date app. Parameters: token (str, optional) - Authentication token for API access, output_path (str, required) - Path for output files, output_file (str, required) - Name of output file")
def service_start_date_app(token: str = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_date_app(api, token=token, user_id=user_id, output_path=output_path, output_file=output_file)

""" Genomics Analysis Services """

@tool(name="service_start_genome_assembly_app", 
      description="""Assemble whole genome sequencing (WGS) reads into contigs using various algorithms.
      
      Supports multiple sequencing platforms (Illumina, PacBio, Nanopore) and handles both 
      paired-end and single-end reads. Automatically selects optimal assembly recipe or allows 
      manual selection from Unicycler, Flye, Canu, SPAdes, and others.
      
      Features:
      - Read preprocessing (trimming, normalization, filtering)
      - Post-assembly polishing (Racon for long reads, Pilon for short reads)
      - Quality filtering (minimum contig length and coverage)
      - Support for SRA datasets
      
      Parameters:
      - paired_end_libs: Paired-end read libraries (optional)
      - single_end_libs: Single-end read libraries (optional) 
      - srr_ids: SRA Run IDs (optional)
      - max_bases: Maximum bases before downsampling (default: 10B)
      - recipe: Assembly algorithm (default: 'auto')
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
      """)
def service_start_genome_assembly_app(token: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_ids: List[str] = None, max_bases: int = 10000000000, recipe: str = "auto", racon_iter: int = 2, pilon_iter: int = 2, trim: bool = False, target_depth: int = 200, normalize: bool = False, filtlong: bool = False, genome_size: int = 5000000, min_contig_len: int = 300, min_contig_cov: float = 5.0, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_genome_assembly_app(api, token=token, user_id=user_id, paired_end_libs=paired_end_libs, single_end_libs=single_end_libs, srr_ids=srr_ids, max_bases=max_bases, recipe=recipe, racon_iter=racon_iter, pilon_iter=pilon_iter, trim=trim, target_depth=target_depth, normalize=normalize, filtlong=filtlong, genome_size=genome_size, min_contig_len=min_contig_len, min_contig_cov=min_contig_cov, output_path=output_path, output_file=output_file)

@tool(name="service_start_genome_annotation_app", 
      description="""Call genes and functionally annotate input contig set using RASTtk pipeline.
      
      This tool performs comprehensive genome annotation including gene calling, functional
      annotation, and feature identification. It uses the RASTtk pipeline with customizable
      workflows for different organism types and annotation requirements.
      
      Key Features:
      - Gene calling with multiple algorithms (Glimmer3, Prodigal, GeneMark)
      - Functional annotation using k-mer analysis and similarity searches
      - Feature identification (rRNA, tRNA, CRISPR, prophages, etc.)
      - Support for Bacteria, Archaea, and Viruses
      - Customizable annotation pipeline
      - Quality analysis and error correction
      - Frameshift detection and correction
      
      Annotation Pipeline Components:
      - rRNA and tRNA identification
      - Repeat region detection
      - CRISPR and prophage identification
      - CDS prediction and annotation
      - Protein functional annotation
      - Overlapping feature resolution
      - Close neighbor analysis
      
      Parameters:
      - contigs: Input contig file (required)
      - scientific_name: Scientific name of organism (required)
      - taxonomy_id: NCBI Taxonomy ID (optional)
      - code: Genetic code for translation (required, default: 0)
      - domain: Organism domain - Bacteria/Archaea/Viruses/auto (required)
      - public: Make genome public (default: False)
      - queue_nowait: Don't wait for indexing (default: False)
      - skip_indexing: Don't index genome (default: False)
      - skip_workspace_output: Don't write to workspace (default: False)
      - output_path: Output directory (optional)
      - output_file: Output basename (optional)
      - reference_genome_id: Reference genome ID (optional)
      - lowvan_min_contig_length: Minimum contig length for lowvan (default: 300)
      - lowvan_max_contig_length: Maximum contig length for lowvan (default: 35000)
      - reference_virus_name: Reference virus name (optional)
      - fix_errors: Automatically fix annotation errors (optional)
      - fix_frameshifts: Fix frameshifts (optional)
      - verbose_level: Verbosity level (optional)
      - workflow: Custom workflow document (optional)
      - recipe: Non-default annotation recipe (optional)
      - disable_replication: Run from scratch even if identical to previous job (optional)
      - analyze_quality: Enable quality analysis (optional)
      - assembly_output: Workspace path for assembly output (optional)
      - custom_pipeline: Customize RASTtk pipeline components (optional)""")
def service_start_genome_annotation_app(token: str = None, genome_id: str = None, contigs: str = None, scientific_name: str = None, tax_id: str = None, my_label: str = None, taxonomy_id: int = None, code: int = 0, domain: str = "auto", public: bool = False, queue_nowait: bool = False, skip_indexing: bool = False, skip_workspace_output: bool = False, output_path: str = None, output_file: str = None, reference_genome_id: str = None, lowvan_min_contig_length: int = 300, lowvan_max_contig_length: int = 35000, reference_virus_name: str = None, fix_errors: bool = None, fix_frameshifts: bool = None, verbose_level: int = None, workflow: str = None, recipe: str = None, disable_replication: bool = None, analyze_quality: bool = None, assembly_output: str = None, custom_pipeline: Dict = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_genome_annotation_app(api, token=token, user_id=user_id, genome_id=genome_id, contigs=contigs, scientific_name=scientific_name, tax_id=tax_id, my_label=my_label, taxonomy_id=taxonomy_id, code=code, domain=domain, public=public, queue_nowait=queue_nowait, skip_indexing=skip_indexing, skip_workspace_output=skip_workspace_output, output_path=output_path, output_file=output_file, reference_genome_id=reference_genome_id, lowvan_min_contig_length=lowvan_min_contig_length, lowvan_max_contig_length=lowvan_max_contig_length, reference_virus_name=reference_virus_name, fix_errors=fix_errors, fix_frameshifts=fix_frameshifts, verbose_level=verbose_level, workflow=workflow, recipe=recipe, disable_replication=disable_replication, analyze_quality=analyze_quality, assembly_output=assembly_output, custom_pipeline=custom_pipeline)

@tool(name="service_start_comprehensive_genome_analysis_app", 
      description="""Perform comprehensive genome analysis from reads, contigs, or GenBank files.
      
      This tool provides end-to-end genome analysis including assembly (if needed), annotation,
      and quality assessment. It can process multiple input types and automatically determines
      the best analysis pipeline based on the input data.
      
      Key Features:
      - Multi-input support: reads, contigs, or GenBank files
      - Automatic assembly with multiple algorithms (Unicycler, Canu, SPAdes, Flye)
      - Comprehensive genome annotation with taxonomic classification
      - Quality analysis and assessment
      - Support for multiple sequencing platforms (Illumina, PacBio, Nanopore)
      - Read preprocessing and polishing options
      
      Input Types:
      - Reads: WGS reads for de novo assembly and annotation
      - Contigs: Pre-assembled contigs for annotation only
      - GenBank: Pre-annotated genomes for quality analysis
      
      Parameters:
      - input_type: Input type - reads/contigs/genbank (required)
      - output_path: Output directory (required)
      - output_file: Output basename (required)
      - paired_end_libs: Paired-end read libraries (optional)
      - single_end_libs: Single-end read libraries (optional)
      - srr_ids: SRA Run IDs (optional)
      - reference_assembly: Reference contig file (optional)
      - recipe: Assembly algorithm (default: 'auto')
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
      """)
def service_start_comprehensive_genome_analysis_app(token: str = None, input_type: str = None, output_path: str = None, output_file: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_ids: List[str] = None, reference_assembly: str = None, recipe: str = "auto", racon_iter: int = 2, pilon_iter: int = 2, trim: bool = False, normalize: bool = False, filtlong: bool = False, target_depth: int = 200, genome_size: int = 5000000, min_contig_len: int = 300, min_contig_cov: float = 5.0, gto: str = None, genbank_file: str = None, contigs: str = None, scientific_name: str = None, taxonomy_id: int = None, code: int = 0, domain: str = "auto", public: bool = False, queue_nowait: bool = False, skip_indexing: bool = False, reference_genome_id: str = None, analyze_quality: bool = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_comprehensive_genome_analysis_app(api, token=token, user_id=user_id, input_type=input_type, output_path=output_path, output_file=output_file, paired_end_libs=paired_end_libs, single_end_libs=single_end_libs, srr_ids=srr_ids, reference_assembly=reference_assembly, recipe=recipe, racon_iter=racon_iter, pilon_iter=pilon_iter, trim=trim, normalize=normalize, filtlong=filtlong, target_depth=target_depth, genome_size=genome_size, min_contig_len=min_contig_len, min_contig_cov=min_contig_cov, gto=gto, genbank_file=genbank_file, contigs=contigs, scientific_name=scientific_name, taxonomy_id=taxonomy_id, code=code, domain=domain, public=public, queue_nowait=queue_nowait, skip_indexing=skip_indexing, reference_genome_id=reference_genome_id, analyze_quality=analyze_quality)

@tool(name="service_start_blast_app", 
      description="""Perform homology searches on sequence data using BLAST algorithms.
      
      This tool performs comprehensive homology searches using various BLAST programs
      (blastp, blastn, blastx, tblastn, tblastx) against multiple database types and
      sources. Supports flexible input and database configurations for different
      sequence analysis scenarios.
      
      Key Features:
      - Multiple BLAST programs: blastp, blastn, blastx, tblastn, tblastx
      - Flexible input sources: FASTA data, files, feature groups, genome groups
      - Multiple database types: protein (faa), DNA (ffn), RNA (frn), contigs (fna)
      - Various database sources: custom sequences, genome lists, taxon lists
      - Configurable search parameters: E-value, coverage, max hits
      - Support for precomputed databases
      
      Input Options:
      - FASTA data: Direct sequence input
      - FASTA files: Workspace file references
      - Feature groups: Annotated feature collections
      - Genome groups: Complete genome collections
      - ID lists: Sequence identifier lists
      
      Database Options:
      - Custom FASTA databases
      - Genome collections
      - Taxonomic groups
      - Feature collections
      - Precomputed databases
      
      Parameters:
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
      - output_file: Output basename (required)""")
def service_start_blast_app(token: str = None, input_type: str = None, input_source: str = None, input_fasta_data: str = None, input_id_list: List[str] = None, input_fasta_file: str = None, input_feature_group: str = None, input_genome_group: str = None, db_type: str = None, db_source: str = None, db_fasta_data: str = None, db_fasta_file: str = None, db_id_list: List[str] = None, db_feature_group: str = None, db_genome_group: str = None, db_genome_list: List[str] = None, db_taxon_list: List[str] = None, db_precomputed_database: str = None, blast_program: str = None, blast_evalue_cutoff: float = 1e-5, blast_max_hits: int = 300, blast_min_coverage: int = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_blast_app(api, token=token, user_id=user_id, input_type=input_type, input_source=input_source, input_fasta_data=input_fasta_data, input_id_list=input_id_list, input_fasta_file=input_fasta_file, input_feature_group=input_feature_group, input_genome_group=input_genome_group, db_type=db_type, db_source=db_source, db_fasta_data=db_fasta_data, db_fasta_file=db_fasta_file, db_id_list=db_id_list, db_feature_group=db_feature_group, db_genome_group=db_genome_group, db_genome_list=db_genome_list, db_taxon_list=db_taxon_list, db_precomputed_database=db_precomputed_database, blast_program=blast_program, blast_evalue_cutoff=blast_evalue_cutoff, blast_max_hits=blast_max_hits, blast_min_coverage=blast_min_coverage, output_path=output_path, output_file=output_file)

@tool(name="service_start_primer_design_app", 
      description="""Use Primer3 to design primers for given DNA sequences.
      
      This tool designs PCR primers using the Primer3 algorithm with comprehensive
      parameter control for primer selection, melting temperature optimization,
      and product size constraints. Supports multiple input formats and advanced
      primer design criteria.
      
      Key Features:
      - Multiple input formats: text sequence, workspace FASTA, or database ID
      - Advanced primer design with Primer3 algorithm
      - Melting temperature optimization and GC content control
      - Product size range specification
      - Target region and exclusion zone definition
      - Internal oligo (probe) design capability
      - Salt concentration and buffer condition optimization
      
      Design Parameters:
      - Primer length control (min/opt/max)
      - Melting temperature constraints (min/opt/max)
      - GC content requirements (min/opt/max)
      - Product size range specification
      - Target and exclusion region definition
      - Buffer condition optimization
      
      Parameters:
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
      - PRIMER_DNTP_CONC: dNTP concentration (mM) (optional)""")
def service_start_primer_design_app(token: str = None, output_file: str = None, output_path: str = None, input_type: str = None, sequence_input: str = None, SEQUENCE_ID: str = None, SEQUENCE_TARGET: List[List[int]] = None, SEQUENCE_INCLUDED_REGION: List[int] = None, SEQUENCE_EXCLUDED_REGION: List[int] = None, SEQUENCE_OVERLAP_JUNCTION_LIST: List[List[int]] = None, PRIMER_PICK_INTERNAL_OLIGO: int = None, PRIMER_PRODUCT_SIZE_RANGE: List[List[int]] = None, PRIMER_NUM_RETURN: int = None, PRIMER_MIN_SIZE: int = None, PRIMER_OPT_SIZE: int = None, PRIMER_MAX_SIZE: int = None, PRIMER_MAX_TM: float = None, PRIMER_MIN_TM: float = None, PRIMER_OPT_TM: float = None, PRIMER_PAIR_MAX_DIFF_TM: float = None, PRIMER_MAX_GC: float = None, PRIMER_MIN_GC: float = None, PRIMER_OPT_GC: float = None, PRIMER_SALT_MONOVALENT: float = None, PRIMER_SALT_DIVALENT: float = None, PRIMER_DNA_CONC: float = None, PRIMER_DNTP_CONC: float = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_primer_design_app(api, token=token, user_id=user_id, output_file=output_file, output_path=output_path, input_type=input_type, sequence_input=sequence_input, SEQUENCE_ID=SEQUENCE_ID, SEQUENCE_TARGET=SEQUENCE_TARGET, SEQUENCE_INCLUDED_REGION=SEQUENCE_INCLUDED_REGION, SEQUENCE_EXCLUDED_REGION=SEQUENCE_EXCLUDED_REGION, SEQUENCE_OVERLAP_JUNCTION_LIST=SEQUENCE_OVERLAP_JUNCTION_LIST, PRIMER_PICK_INTERNAL_OLIGO=PRIMER_PICK_INTERNAL_OLIGO, PRIMER_PRODUCT_SIZE_RANGE=PRIMER_PRODUCT_SIZE_RANGE, PRIMER_NUM_RETURN=PRIMER_NUM_RETURN, PRIMER_MIN_SIZE=PRIMER_MIN_SIZE, PRIMER_OPT_SIZE=PRIMER_OPT_SIZE, PRIMER_MAX_SIZE=PRIMER_MAX_SIZE, PRIMER_MAX_TM=PRIMER_MAX_TM, PRIMER_MIN_TM=PRIMER_MIN_TM, PRIMER_OPT_TM=PRIMER_OPT_TM, PRIMER_PAIR_MAX_DIFF_TM=PRIMER_PAIR_MAX_DIFF_TM, PRIMER_MAX_GC=PRIMER_MAX_GC, PRIMER_MIN_GC=PRIMER_MIN_GC, PRIMER_OPT_GC=PRIMER_OPT_GC, PRIMER_SALT_MONOVALENT=PRIMER_SALT_MONOVALENT, PRIMER_SALT_DIVALENT=PRIMER_SALT_DIVALENT, PRIMER_DNA_CONC=PRIMER_DNA_CONC, PRIMER_DNTP_CONC=PRIMER_DNTP_CONC)

@tool(name="service_start_similar_genome_finder_app", 
      description="""Find similar genomes based on genomic similarity analysis.
      
      This tool identifies genomes with high similarity to a reference genome using
      various genomic distance metrics and statistical methods. It supports filtering
      by organism type, reference status, and similarity thresholds.
      
      Parameters:
      - selectedGenomeId: Reference genome ID for similarity search (optional)
      - fasta_file: FASTA/FASTQ file for similarity search (optional)
      - max_pvalue: Maximum p-value threshold for similarity (optional)
      - max_distance: Maximum distance threshold for similarity (optional) 
      - max_hits: Maximum number of similar genomes to return (optional)
      - include_reference: Include reference genomes in results (optional)
      - include_representative: Include representative genomes in results (optional)
      - include_bacterial: Include bacterial genomes in results (optional)
      - include_viral: Include viral genomes in results (optional)
      - output_path: Output directory path (required)
      - output_file: Output file basename (required)
      """)
def service_start_similar_genome_finder_app(token: str = None, selectedGenomeId: str = None, fasta_file: str = None, max_pvalue: float = None, max_distance: float = None, max_hits: int = None, include_reference: bool = None, include_representative: bool = None, include_bacterial: bool = None, include_viral: bool = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    if not selectedGenomeId and not fasta_file:
        return "Error: selectedGenomeId or fasta_file is required"
    return start_similar_genome_finder_app(similar_genome_finder_api, token=token, user_id=user_id, selectedGenomeId=selectedGenomeId, fasta_file=fasta_file, max_pvalue=max_pvalue, max_distance=max_distance, max_hits=max_hits, include_reference=include_reference, include_representative=include_representative, include_bacterial=include_bacterial, include_viral=include_viral, output_path=output_path, output_file=output_file)

@tool(name="service_start_genome_alignment_app", 
      description="""Perform multiple whole genome alignment with rearrangements using Mauve.
      
      This tool performs comprehensive whole genome alignment analysis using the Mauve
      algorithm to identify conserved regions, rearrangements, and structural variations
      across multiple genomes. It supports advanced alignment parameters and statistical
      models for robust genome comparison and evolutionary analysis.
      
      Key Features:
      - Multiple whole genome alignment with Mauve algorithm
      - Detection of genomic rearrangements and structural variations
      - Support for progressiveMauve and mauveAligner methods
      - Advanced statistical parameters for alignment optimization
      - HMM-based homologous region detection
      - Breakpoint penalty and conservation distance scaling
      - Configurable seed weight and gapped aligner parameters
      
      Alignment Methods:
      - progressiveMauve: Progressive multiple genome alignment
      - mauveAligner: Standard Mauve alignment algorithm
      
      Statistical Parameters:
      - Seed weight: Controls initial anchor calculation
      - Gapped aligner length: Maximum base pairs for gapped alignment
      - Breakpoint distance scaling: Weight scaling by breakpoint distance (0-1)
      - Conservation distance scaling: Scale conservation distances (0-1)
      - Minimum pairwise LCB score: Threshold for locally collinear blocks
      - Scaled penalty: Minimum breakpoint penalty after scaling
      
      HMM Parameters:
      - P(Go|Homologous): Transition probability from unrelated to homologous state
      - P(Go|Unrelated): Transition probability from homologous to unrelated state
      
      Parameters:
      - genome_ids: List of genome IDs to align (required)
      - recipe: Mauve method - progressiveMauve/mauveAligner (default: progressiveMauve)
      - seedWeight: Seed weight for calculating initial anchors (optional)
      - maxGappedAlignerLength: Maximum base pairs for gapped aligner (optional)
      - maxBreakpointDistanceScale: Maximum weight scaling by breakpoint distance (0-1) (optional)
      - conservationDistanceScale: Scale conservation distances (0-1) (optional)
      - weight: Minimum pairwise LCB score (optional)
      - minScaledPenalty: Minimum breakpoint penalty after scaling (optional)
      - hmmPGoHomologous: Probability of transitioning from unrelated to homologous state (optional)
      - hmmPGoUnrelated: Probability of transitioning from homologous to unrelated state (optional)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_genome_alignment_app(token: str = None, genome_ids: List[str] = None, recipe: str = "progressiveMauve", seedWeight: float = None, maxGappedAlignerLength: float = None, maxBreakpointDistanceScale: float = None, conservationDistanceScale: float = None, weight: float = None, minScaledPenalty: float = None, hmmPGoHomologous: float = None, hmmPGoUnrelated: float = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_genome_alignment_app(api, token=token, user_id=user_id, genome_ids=genome_ids, recipe=recipe, seedWeight=seedWeight, maxGappedAlignerLength=maxGappedAlignerLength, maxBreakpointDistanceScale=maxBreakpointDistanceScale, conservationDistanceScale=conservationDistanceScale, weight=weight, minScaledPenalty=minScaledPenalty, hmmPGoHomologous=hmmPGoHomologous, hmmPGoUnrelated=hmmPGoUnrelated, output_path=output_path, output_file=output_file)

@tool(name="service_start_variation_app", 
      description="""Identify and annotate small nucleotide variations (SNVs) relative to a reference genome.
      
      This tool performs comprehensive variant calling analysis using multiple mapping
      and variant calling algorithms. It identifies single nucleotide variants (SNVs)
      and other small variations by comparing sequencing reads against a reference
      genome using state-of-the-art bioinformatics tools.
      
      Key Features:
      - Multiple short read mappers: BWA-mem, Bowtie2, MOSAIK, LAST, minimap2
      - Advanced variant callers: FreeBayes, BCFtools
      - Support for paired-end and single-end reads
      - SRA dataset integration
      - Flexible mapping and calling parameters
      - Comprehensive variant annotation
      
      Mapping Options:
      - BWA-mem: Fast and accurate alignment
      - BWA-mem-strict: Strict alignment parameters
      - Bowtie2: Versatile aligner with sensitive options
      - MOSAIK: High-quality alignment
      - LAST: Long read alignment
      - minimap2: Fast long read mapping
      
      Variant Calling:
      - FreeBayes: Bayesian variant caller
      - BCFtools: Comprehensive variant calling suite
      
      Parameters:
      - reference_genome_id: Reference genome ID (required)
      - paired_end_libs: Paired-end read libraries (optional)
      - single_end_libs: Single-end read libraries (optional)
      - srr_ids: SRA Run IDs (optional)
      - mapper: Short read mapper - BWA-mem/BWA-mem-strict/Bowtie2/MOSAIK/LAST/minimap2 (default: BWA-mem)
      - caller: SNP caller - FreeBayes/BCFtools (default: FreeBayes)
      - output_path: Output directory (required)
      - output_file: Output basename (required)
      """)
def service_start_variation_app(token: str = None, reference_genome_id: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_ids: List[str] = None, mapper: str = "BWA-mem", caller: str = "FreeBayes", output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_variation_app(api, token=token, user_id=user_id, reference_genome_id=reference_genome_id, paired_end_libs=paired_end_libs, single_end_libs=single_end_libs, srr_ids=srr_ids, mapper=mapper, caller=caller, output_path=output_path, output_file=output_file)

@tool(name="service_start_tnseq_app", 
      description="""Analyze TnSeq (Transposon Sequencing) data using TRANSIT pipeline.
      
      This tool performs comprehensive analysis of transposon insertion sequencing data
      to identify essential genes and functional elements in bacterial genomes. It uses
      the TRANSIT pipeline with multiple statistical methods for robust gene essentiality
      analysis and comparative studies.
      
      Key Features:
      - Multiple analysis recipes: Gumbel, Griffin, Tn5gaps, Rank Product, HMM, Binomial, Resampling
      - Support for different TnSeq protocols: Sassetti, Tn5, MME1
      - Comparative analysis with experimental conditions and contrasts
      - Read file processing and primer trimming
      - Statistical significance testing
      - Essential gene identification
      
      Analysis Recipes:
      - Gumbel: Maximum likelihood estimation of essentiality
      - Griffin: Bayesian approach for essential gene detection
      - Tn5gaps: Gap-based analysis for Tn5 insertions
      - Rank Product: Non-parametric ranking method
      - HMM: Hidden Markov Model approach
      - Binomial: Statistical testing for essentiality
      - Resampling: Bootstrap resampling analysis
      
      Protocols:
      - Sassetti: Standard TnSeq protocol
      - Tn5: Tn5 transposon protocol
      - MME1: MME1 transposon protocol
      
      Parameters:
      - experimental_conditions: Experimental conditions for analysis (optional)
      - contrasts: Statistical contrasts for comparison (optional)
      - read_files: Read file groups (optional)
      - reference_genome_id: Reference genome ID (required)
      - recipe: TnSeq analysis recipe - gumbel/griffin/tn5gaps/rankproduct/hmm/binomial/resampling (default: gumbel)
      - protocol: TnSeq protocol - sassetti/tn5/mme1 (default: sassetti)
      - primer: Primer DNA string for read trimming (optional)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_tnseq_app(token: str = None, experimental_conditions: List[str] = None, contrasts: List[str] = None, read_files: List[Dict] = None, reference_genome_id: str = None, recipe: str = "gumbel", protocol: str = "sassetti", primer: str = "", output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_tnseq_app(api, token=token, user_id=user_id, experimental_conditions=experimental_conditions, contrasts=contrasts, read_files=read_files, reference_genome_id=reference_genome_id, recipe=recipe, protocol=protocol, primer=primer, output_path=output_path, output_file=output_file)

""" Phylogenomics Services """

@tool(name="service_start_bacterial_genome_tree_app", 
      description="""Compute phylogenetic tree from PGFam protein and DNA sequences.
      
      This tool constructs phylogenetic trees based on protein and DNA sequences
      of PGFams (Protein Family Groups) for a set of genomes. It uses codon-based
      phylogenetic analysis to build robust evolutionary relationships between
      bacterial genomes using shared protein families.
      
      Key Features:
      - PGFam-based phylogenetic analysis using protein and DNA sequences
      - Codon-aware tree construction for accurate evolutionary inference
      - Bootstrap support for tree reliability assessment
      - Flexible genome selection (main and optional genomes)
      - Gene family filtering and quality control
      - Metadata integration for genome annotation
      
      Analysis Options:
      - Main genomes: Primary genomes for tree construction
      - Optional genomes: Additional genomes (not penalized for missing genes)
      - Gene filtering: Control for missing genes and duplications
      - Bootstrap analysis: Statistical support for tree branches
      - Metadata fields: Genome information integration
      
      Quality Control:
      - Maximum missing genomes per PGFam
      - Maximum duplications allowed per PGFam
      - Desired number of genes for analysis
      - Bootstrap replicate control
      
      Parameters:
      - output_path: Output directory (required)
      - output_file: Output basename (required)
      - genome_ids: Main genomes for analysis (optional)
      - genome_groups: Main genome groups (optional)
      - optional_genome_ids: Optional genomes (not penalized for missing genes) (optional)
      - genome_metadata_fields: Genome metadata fields to retrieve (optional)
      - number_of_genes: Desired number of genes for analysis (default: 20)
      - bootstraps: Number of bootstrap replicates (default: 100)
      - max_genomes_missing: Max main genomes missing from any PGFam (default: 0)
      - max_allowed_dups: Max duplications allowed for main genomes (default: 0)""")
def service_start_bacterial_genome_tree_app(token: str = None, output_path: str = None, output_file: str = None, genome_ids: List[str] = None, genome_groups: List[str] = None, optional_genome_ids: List[str] = None, genome_metadata_fields: str = None, number_of_genes: int = 20, bootstraps: int = 100, max_genomes_missing: int = 0, max_allowed_dups: int = 0) -> str:
    user_id = extract_userid_from_token(token)
    return start_bacterial_genome_tree_app(api, token=token, user_id=user_id, output_path=output_path, output_file=output_file, genome_ids=genome_ids, genome_groups=genome_groups, optional_genome_ids=optional_genome_ids, genome_metadata_fields=genome_metadata_fields, number_of_genes=number_of_genes, bootstraps=bootstraps, max_genomes_missing=max_genomes_missing, max_allowed_dups=max_allowed_dups)

@tool(name="service_start_viral_genome_tree_app", description="")
def service_start_viral_genome_tree_app(token: str = None) -> str:
    return null

@tool(name="service_start_gene_tree_app", 
      description="""Estimate phylogeny of gene or other sequence features.
      
      This tool constructs phylogenetic trees for genes or sequence features using
      multiple alignment and tree-building algorithms. It supports both DNA and
      protein sequences with various substitution models and tree construction
      methods for robust phylogenetic inference.
      
      Key Features:
      - Multiple sequence input formats: aligned DNA/protein FASTA, feature groups, genome groups
      - Advanced alignment programs: MUSCLE, MAFFT
      - Multiple tree-building algorithms: RAxML, PhyML, FastTree
      - Comprehensive substitution models for DNA and protein sequences
      - Bootstrap support for tree reliability assessment
      - Alignment quality control and trimming
      - Metadata integration for genes and genomes
      
      Alignment Options:
      - MUSCLE: Multiple sequence alignment by log-expectation
      - MAFFT: Multiple alignment using fast Fourier transform
      - End trimming: Remove poorly aligned regions
      - Gap filtering: Remove sequences with excessive gaps
      
      Tree Construction:
      - RAxML: Randomized axelerated maximum likelihood
      - PhyML: Phylogenetic analysis by maximum likelihood
      - FastTree: Approximate maximum likelihood trees
      
      Substitution Models:
      - DNA: HKY85, JC69, K80, F81, F84, TN93, GTR
      - Protein: LG, WAG, JTT, MtREV, Dayhoff, DCMut, RtREV, CpREV, VT, AB, Blosum62, MtMam, MtArt, HIVw, HIVb
      
      Parameters:
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
      - output_file: Output basename (required)""")
def service_start_gene_tree_app(token: str = None, sequences: List[str] = None, alignment_program: str = None, trim_threshold: float = None, gap_threshold: float = None, alphabet: str = None, substitution_model: str = None, bootstrap: int = None, recipe: str = "RAxML", tree_type: str = None, feature_metadata_fields: str = None, genome_metadata_fields: str = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_gene_tree_app(api, token=token, user_id=user_id, sequences=sequences, alignment_program=alignment_program, trim_threshold=trim_threshold, gap_threshold=gap_threshold, alphabet=alphabet, substitution_model=substitution_model, bootstrap=bootstrap, recipe=recipe, tree_type=tree_type, feature_metadata_fields=feature_metadata_fields, genome_metadata_fields=genome_metadata_fields, output_path=output_path, output_file=output_file)

@tool(name="service_start_core_genome_mlst_app", 
      description="""Evaluate core genomes from a set of genome groups of the same species.
      
      This tool performs core genome Multi-Locus Sequence Typing (cgMLST) analysis
      using the chewBBACA pipeline to identify and analyze core genome loci across
      multiple genomes of the same species. It provides standardized typing schemes
      for bacterial strain classification and epidemiological studies.
      
      Key Features:
      - Core genome identification and analysis
      - chewBBACA-based SNP calling and analysis
      - Species-specific schema support
      - Multiple input formats: genome groups and FASTA files
      - Standardized MLST typing schemes
      - Comparative genomic analysis
      
      Analysis Pipeline:
      - Core genome extraction from multiple genomes
      - Locus identification and alignment
      - SNP calling and variant analysis
      - MLST profile generation
      - Phylogenetic relationship inference
      
      Input Options:
      - Genome groups: Pre-defined genome collections
      - FASTA files: Direct genome sequence input
      - Schema selection: Species-specific typing schemes
      
      Parameters:
      - input_genome_type: Input type - genome_group/genome_fasta (required, default: genome_group)
      - analysis_type: Analysis type - chewbbaca (required, default: chewbbaca)
      - input_genome_group: Genome group name (optional)
      - input_genome_fasta: Genome FASTA data (optional)
      - schema_location: Path to schema directory (optional)
      - input_schema_selection: Schema species selection (required)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_core_genome_mlst_app(token: str = None, input_genome_type: str = "genome_group", analysis_type: str = "chewbbaca", input_genome_group: str = None, input_genome_fasta: str = None, schema_location: str = None, input_schema_selection: str = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_core_genome_mlst_app(api, token=token, user_id=user_id, input_genome_type=input_genome_type, analysis_type=analysis_type, input_genome_group=input_genome_group, input_genome_fasta=input_genome_fasta, schema_location=schema_location, input_schema_selection=input_schema_selection, output_path=output_path, output_file=output_file)

@tool(name="service_start_whole_genome_snp_app", 
      description="""Identify SNP differences in a genome group with genomes of the same species.
      
      This tool performs comprehensive whole genome SNP analysis to identify single
      nucleotide polymorphisms across multiple genomes of the same species. It uses
      advanced algorithms to detect and analyze genetic variations for comparative
      genomics and phylogenetic studies.
      
      Key Features:
      - Whole genome SNP identification and analysis
      - Multiple analysis algorithms: chewBBACA, kSNP4
      - Flexible input formats: genome groups and FASTA files
      - Linkage analysis with configurable thresholds
      - Majority threshold filtering for robust SNP detection
      - Comparative genomic analysis
      - Phylogenetic relationship inference
      
      Analysis Pipeline:
      - Genome alignment and comparison
      - SNP detection and validation
      - Linkage analysis and filtering
      - Phylogenetic tree construction
      - Statistical analysis and reporting
      
      Linkage Analysis:
      - Minimum mid-linkage: Threshold for strong linkage (default: 10)
      - Maximum mid-linkage: Threshold for weak linkage (default: 40)
      - Configurable SNP distance thresholds
      - Genome difference rating system
      
      Majority Threshold Options:
      - 0.0 to 1.0: Fraction of genomes required for SNP inclusion
      - Tree construction based on consensus SNPs
      - Quality control for SNP reliability
      
      Parameters:
      - input_genome_type: Input type - genome_group/genome_fasta (required)
      - majority_threshold: Majority fraction for SNP inclusion (default: 0.5)
      - min_mid_linkage: Minimum value for mid linkage (default: 10)
      - max_mid_linkage: Maximum value for mid linkage (default: 40)
      - analysis_type: Analysis type - Whole Genome SNP Analysis (required)
      - input_genome_group: Genome group name (optional)
      - input_genome_fasta: Genome FASTA data (optional)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_whole_genome_snp_app(token: str = None, input_genome_type: str = None, majority_threshold: float = 0.5, min_mid_linkage: int = 10, max_mid_linkage: int = 40, analysis_type: str = "Whole Genome SNP Analysis", input_genome_group: str = None, input_genome_fasta: str = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_whole_genome_snp_app(api, token=token, user_id=user_id, input_genome_type=input_genome_type, majority_threshold=majority_threshold, min_mid_linkage=min_mid_linkage, max_mid_linkage=max_mid_linkage, analysis_type=analysis_type, input_genome_group=input_genome_group, input_genome_fasta=input_genome_fasta, output_path=output_path, output_file=output_file)

""" Protein Tools """

@tool(name="service_start_msa_snp_analysis_app", 
      description="""Compute multiple sequence alignment and analyze SNP/variance.
      
      This tool performs comprehensive multiple sequence alignment (MSA) analysis
      to identify single nucleotide polymorphisms (SNPs) and other variations
      across multiple sequences. It supports various input formats and alignment
      algorithms for robust sequence comparison and variant analysis.
      
      Key Features:
      - Multiple sequence alignment using MUSCLE, MAFFT, or progressiveMauve
      - Support for DNA and protein sequences
      - Multiple input formats: FASTA files, genome groups, feature groups, ID lists
      - Reference sequence support for alignment optimization
      - SNP and variation analysis
      - Flexible alignment strategies and parameters
      
      Input Options:
      - FASTA files: Direct sequence file input with format specification
      - Genome groups: Pre-defined genome collections
      - Feature groups: Annotated feature collections
      - Feature lists: Direct feature ID collections
      - Genome lists: Direct genome ID collections
      - Keyboard input: Direct FASTA sequence input
      
      Alignment Algorithms:
      - MUSCLE: Multiple sequence alignment by log-expectation
      - MAFFT: Multiple alignment using fast Fourier transform
      - progressiveMauve: Progressive multiple alignment
      
      MAFFT Strategies:
      - auto: Automatic strategy selection
      - fftns1, fftns2, fftnsi: FFT-based methods
      - einsi, linsi, ginsi: Iterative refinement methods
      
      Reference Options:
      - none: No reference sequence
      - string: Direct sequence input
      - feature_id: Feature identifier
      - genome_id: Genome identifier
      - first: Use first sequence as reference
      
      Parameters:
      - input_status: Input alignment status - unaligned/aligned (default: unaligned)
      - input_type: Input type - input_group/input_fasta/input_sequence/input_genomegroup/input_featurelist/input_genomelist (default: input_group)
      - fasta_files: FASTA file groups with file and type (optional)
      - select_genomegroup: Genome groups for analysis (optional)
      - feature_groups: Feature groups for analysis (optional)
      - feature_list: Feature ID lists (optional)
      - genome_list: Genome ID lists (optional)
      - aligner: Alignment algorithm - Muscle/Mafft/progressiveMauve (default: Muscle)
      - alphabet: Sequence alphabet - dna/protein (default: dna)
      - fasta_keyboard_input: Direct FASTA sequence input (optional)
      - ref_type: Reference type - none/string/feature_id/genome_id/first (default: none)
      - strategy: MAFFT strategy - auto/fftns1/fftns2/fftnsi/einsi/linsi/ginsi (default: auto)
      - ref_string: Reference sequence identity (optional)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_msa_snp_analysis_app(token: str = None, input_status: str = "unaligned", input_type: str = "input_group", fasta_files: List[Dict] = None, select_genomegroup: List[str] = None, feature_groups: List[str] = None, feature_list: List[str] = None, genome_list: List[str] = None, aligner: str = "Muscle", alphabet: str = "dna", fasta_keyboard_input: str = "", ref_type: str = "none", strategy: str = "auto", ref_string: str = "", output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_msa_snp_analysis_app(api, token=token, user_id=user_id, input_status=input_status, input_type=input_type, fasta_files=fasta_files, select_genomegroup=select_genomegroup, feature_groups=feature_groups, feature_list=feature_list, genome_list=genome_list, aligner=aligner, alphabet=alphabet, fasta_keyboard_input=fasta_keyboard_input, ref_type=ref_type, strategy=strategy, ref_string=ref_string, output_path=output_path, output_file=output_file)

@tool(name="service_start_metacats_app", 
      description="""Metadata-driven Comparative Analysis Tool (meta-CATS) for identifying significant differences between sequence groups.
      
      This tool performs comprehensive comparative analysis to identify positions that
      significantly differ between user-defined groups of sequences. It uses statistical
      methods to detect variations and associations across different sequence groups
      based on metadata and alignment data.
      
      Key Features:
      - Statistical analysis of sequence positions between groups
      - Support for nucleotide and amino acid sequences
      - Multiple input types: feature groups, alignment files, group files
      - Metadata-driven group definition and analysis
      - Configurable p-value thresholds for significance testing
      - Year range filtering and metadata grouping
      - Flexible group assignment and analysis
      
      Input Types:
      - auto: Automatic input type detection
      - groups: Feature group-based analysis
      - files: File-based analysis with alignment and group files
      
      Sequence Alphabets:
      - na: Nucleotide sequences (DNA/RNA)
      - aa: Amino acid sequences (proteins)
      
      Analysis Options:
      - Feature groups: Pre-defined feature collections
      - Alignment files: Direct sequence alignment input
      - Group files: TSV files defining sequence groups
      - Auto groups: Manual group assignment with metadata
      
      Statistical Parameters:
      - P-value cutoff: Significance threshold for analysis (default: 0.05)
      - Year ranges: Temporal filtering for sequences
      - Metadata grouping: Category-based analysis
      
      Group Assignment:
      - Feature/sequence IDs: Unique identifiers for sequences
      - Group strings: Descriptive group assignments
      - Genome IDs: Associated genome identifiers
      - Metadata categories: Original metadata classification
      
      Parameters:
      - output_path: Output directory (required)
      - output_file: Output basename (required)
      - p_value: P-value cutoff for significance testing (default: 0.05)
      - year_ranges: Year ranges for temporal filtering (optional)
      - metadata_group: Metadata type for reference (optional)
      - input_type: Input type - auto/groups/files (required)
      - alphabet: Sequence alphabet - na/aa (default: na)
      - groups: Feature groups for analysis (optional)
      - alignment_file: Alignment file location (optional)
      - group_file: Group file location (optional)
      - alignment_type: Alignment format - aligned_dna_fasta/aligned_protein_fasta (optional)
      - auto_groups: Manual group assignments with id, grp, g_id, metadata (optional)""")
def service_start_metacats_app(token: str = None, output_path: str = None, output_file: str = None, p_value: float = 0.05, year_ranges: str = None, metadata_group: str = None, input_type: str = None, alphabet: str = "na", groups: List[str] = None, alignment_file: str = None, group_file: str = None, alignment_type: str = None, auto_groups: List[Dict] = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_metacats_app(api, token=token, user_id=user_id, output_path=output_path, output_file=output_file, p_value=p_value, year_ranges=year_ranges, metadata_group=metadata_group, input_type=input_type, alphabet=alphabet, groups=groups, alignment_file=alignment_file, group_file=group_file, alignment_type=alignment_type, auto_groups=auto_groups)

@tool(name="service_start_proteome_comparison_app", 
      description="""Compare proteome sets from multiple genomes using BLAST-based analysis.
      
      This tool performs comprehensive proteome comparison across multiple genomes
      using BLAST algorithms to identify homologous proteins, functional similarities,
      and evolutionary relationships. It supports various input formats and provides
      detailed comparative analysis with configurable similarity thresholds.
      
      Key Features:
      - BLAST-based proteome comparison across multiple genomes
      - Support for genome IDs, user genomes, and feature groups
      - Configurable reference genome selection
      - Advanced BLAST parameters for similarity detection
      - Coverage and identity threshold controls
      - E-value and positive scoring position filtering
      - Comprehensive comparative analysis results
      
      Input Options:
      - Genome IDs: Pre-defined genome collections from database
      - User genomes: Custom protein sequence files in FASTA format
      - Feature groups: Annotated feature collections
      - Reference genome: Configurable reference genome selection
      
      BLAST Parameters:
      - Minimum sequence coverage: Threshold for query and subject coverage (default: 0.30)
      - Maximum E-value: Statistical significance threshold (default: 1e-5)
      - Minimum identity: Fraction identity threshold (default: 0.1)
      - Minimum positives: Positive-scoring positions threshold (default: 0.2)
      
      Analysis Features:
      - Homologous protein identification
      - Functional similarity assessment
      - Evolutionary relationship analysis
      - Comparative proteome statistics
      - Reference-based comparison framework
      
      Parameters:
      - genome_ids: Genome IDs for comparison (optional)
      - user_genomes: Custom protein sequence files in FASTA format (optional)
      - user_feature_groups: Feature groups for analysis (optional)
      - reference_genome_index: Index of reference genome (1-based) (default: 1)
      - min_seq_cov: Minimum coverage of query and subject (default: 0.30)
      - max_e_val: Maximum E-value threshold (default: 1e-5)
      - min_ident: Minimum fraction identity (default: 0.1)
      - min_positives: Minimum fraction positive-scoring positions (default: 0.2)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_proteome_comparison_app(token: str = None, genome_ids: List[str] = None, user_genomes: List[str] = None, user_feature_groups: List[str] = None, reference_genome_index: int = 1, min_seq_cov: float = 0.30, max_e_val: float = 1e-5, min_ident: float = 0.1, min_positives: float = 0.2, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_proteome_comparison_app(api, token=token, user_id=user_id, genome_ids=genome_ids, user_genomes=user_genomes, user_feature_groups=user_feature_groups, reference_genome_index=reference_genome_index, min_seq_cov=min_seq_cov, max_e_val=max_e_val, min_ident=min_ident, min_positives=min_positives, output_path=output_path, output_file=output_file)

@tool(name="service_start_protein_family_sorter_app", description="")
def service_start_protein_family_sorter_app(token: str = None) -> str:
    return None

@tool(name="service_start_comparative_systems_app", 
      description="""Create data structures to decompose genomes for comparative analysis.
      
      This tool creates specialized data structures to decompose and analyze genomes
      for comparative systems biology studies. It processes genome collections and
      generates structured data representations that facilitate downstream comparative
      analysis and systems-level investigations.
      
      Key Features:
      - Genome decomposition and data structure creation
      - Support for individual genome IDs and genome groups
      - Structured data representation for comparative analysis
      - Systems-level genome organization and analysis
      - Data structure optimization for large-scale comparisons
      - Flexible input handling for various genome collections
      
      Input Options:
      - Genome IDs: Individual genome identifiers for analysis
      - Genome Groups: Pre-defined genome collections
      - Mixed inputs: Combination of individual genomes and groups
      
      Analysis Capabilities:
      - Genome structure decomposition
      - Comparative data structure generation
      - Systems-level organization
      - Data structure optimization
      - Comparative framework preparation
      
      Use Cases:
      - Systems biology research
      - Comparative genomics studies
      - Large-scale genome analysis
      - Data structure preparation for downstream tools
      - Genome collection organization
      
      Parameters:
      - output_path: Output directory (required)
      - output_file: Output basename (required)
      - genome_ids: Individual genome IDs for analysis (optional)
      - genome_groups: Genome groups for analysis (optional)""")
def service_start_comparative_systems_app(token: str = None, output_path: str = None, output_file: str = None, genome_ids: List[str] = None, genome_groups: List[str] = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_comparative_systems_app(api, token=token, user_id=user_id, output_path=output_path, output_file=output_file, genome_ids=genome_ids, genome_groups=genome_groups)

@tool(name="service_start_docking_app", 
      description="""Dock small molecules into protein structures using advanced molecular docking algorithms.
      
      This tool performs comprehensive molecular docking analysis to predict the binding
      of small molecules to protein structures. It supports multiple input formats for
      proteins and ligands, with configurable docking parameters and result filtering.
      
      Key Features:
      - Small molecule docking into protein structures
      - Multiple protein input types: PDB identifiers and user-provided files
      - Flexible ligand library options: named libraries, SMILES lists, workspace files
      - Advanced DiffDock algorithm with configurable batch processing
      - Top-N result filtering and ranking
      - Debug mode for detailed analysis output
      - High-performance batch processing optimization
      
      Protein Input Options:
      - input_pdb: PDB identifiers from precomputed structures
      - user_pdb_file: User-provided PDB files from workspace
      
      Ligand Library Options:
      - ws_file: Workspace files containing SMILES strings
      - named_library: Pre-defined ligand libraries
      - smiles_list: Direct SMILES string arrays
      
      Docking Parameters:
      - Batch size: DiffDock batch processing size (default: 10)
      - Top N results: Number of top-ranked results to return
      - Performance optimization: Automatic batch size adjustment based on protein size
      
      Analysis Features:
      - Binding pose prediction
      - Binding affinity estimation
      - Molecular interaction analysis
      - Conformational sampling
      - Result ranking and filtering
      
      Parameters:
      - protein_input_type: Type of PDB input - input_pdb/user_pdb_file (required)
      - input_pdb: PDB identifiers for precomputed structures (optional)
      - user_pdb_file: User-provided PDB files (optional)
      - ligand_library_type: Ligand library type - ws_file/named_library/smiles_list (required)
      - ligand_named_library: Name of ligand library (optional)
      - ligand_smiles_list: Array of SMILES strings (optional)
      - ligand_ws_file: Workspace file with SMILES strings (optional)
      - top_n: Return top N results (optional)
      - batch_size: DiffDock batch size (default: 10)
      - output_path: Output directory (required)
      - output_file: Output basename (required)
      - enable_debug: Enable debug statements (optional)""")
def service_start_docking_app(token: str = None, protein_input_type: str = None, input_pdb: List[str] = None, user_pdb_file: List[str] = None, ligand_library_type: str = None, ligand_named_library: str = None, ligand_smiles_list: List[str] = None, ligand_ws_file: str = None, top_n: int = None, batch_size: int = 10, output_path: str = None, output_file: str = None, enable_debug: bool = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_docking_app(api, token=token, user_id=user_id, protein_input_type=protein_input_type, input_pdb=input_pdb, user_pdb_file=user_pdb_file, ligand_library_type=ligand_library_type, ligand_named_library=ligand_named_library, ligand_smiles_list=ligand_smiles_list, ligand_ws_file=ligand_ws_file, top_n=top_n, batch_size=batch_size, output_path=output_path, output_file=output_file, enable_debug=enable_debug)

""" Metagenomics Services """

@tool(name="service_start_taxonomic_classification_app", 
      description="""Compute taxonomic classification for read data using metagenomic analysis.
      
      This tool performs comprehensive taxonomic classification of sequencing reads
      from metagenomic samples, supporting multiple analysis types including pathogen
      detection, microbiome analysis, and 16S rRNA gene classification. It uses
      advanced classification algorithms with multiple reference databases.
      
      Key Features:
      - Multiple analysis types: pathogen detection, microbiome analysis, 16S rRNA
      - Host genome filtering for various model organisms
      - Multiple input formats: paired-end, single-end, and SRA datasets
      - Comprehensive reference databases: BVBRC, Greengenes, SILVA, standard
      - Confidence interval control for classification accuracy
      - Sequence preservation options for classified/unclassified reads
      
      Analysis Types:
      - Pathogen: Detection and classification of pathogenic organisms
      - Microbiome: Comprehensive microbial community analysis
      - 16S: 16S rRNA gene-based taxonomic classification
      
      Host Genome Options:
      - Human (Homo sapiens)
      - Mouse (Mus musculus)
      - Rat (Rattus norvegicus)
      - C. elegans, Drosophila, Zebrafish
      - Chicken, Macaque, Ferret, Pig
      - No host filtering
      
      Reference Databases:
      - BVBRC: Bacterial and viral reference database
      - Greengenes: 16S rRNA gene database
      - SILVA: Comprehensive ribosomal RNA database
      - Standard: General reference database
      
      Parameters:
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
      - output_file: Output basename (required)""")
def service_start_taxonomic_classification_app(token: str = None, host_genome: str = "no_host", analysis_type: str = "16S", paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_libs: List[Dict] = None, database: str = "SILVA", save_classified_sequences: bool = False, save_unclassified_sequences: bool = False, confidence_interval: float = 0.1, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_taxonomic_classification_app(api, token=token, user_id=user_id, host_genome=host_genome, analysis_type=analysis_type, paired_end_libs=paired_end_libs, single_end_libs=single_end_libs, srr_libs=srr_libs, database=database, save_classified_sequences=save_classified_sequences, save_unclassified_sequences=save_unclassified_sequences, confidence_interval=confidence_interval, output_path=output_path, output_file=output_file)

@tool(name="service_start_metagenomic_binning_app", 
      description="""Assemble, bin, and annotate metagenomic sample data.
      
      This tool performs comprehensive metagenomic analysis including assembly,
      binning, and annotation of microbial communities from sequencing data.
      It supports both bacterial and viral binning with automated annotation
      and genome reconstruction from complex environmental samples.
      
      Key Features:
      - Metagenomic assembly from sequencing reads
      - Automated binning for bacterial and viral genomes
      - Comprehensive annotation of binned genomes
      - Multiple input formats: paired-end, single-end, and SRA datasets
      - Quality filtering and contig optimization
      - Integration with PATRIC database
      - Flexible assembly and annotation options
      
      Analysis Pipeline:
      - Read preprocessing and quality control
      - Metagenomic assembly using various algorithms
      - Contig binning and genome reconstruction
      - Taxonomic classification of bins
      - Functional annotation of reconstructed genomes
      - Database integration and indexing
      
      Binning Options:
      - Bacterial binning: Reconstruction of bacterial genomes
      - Viral binning: Identification and reconstruction of viral genomes
      - Dangling element analysis: k-mer based contig binning
      - Quality filtering: Length and coverage thresholds
      
      Annotation Features:
      - Bacterial genome annotation and loading
      - Viral genome annotation and loading
      - Custom annotation recipes
      - PATRIC database integration
      - Solr indexing for analysis
      
      Parameters:
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
      - min_contig_cov: Minimum contig coverage (default: 4.0)""")
def service_start_metagenomic_binning_app(token: str = None, paired_end_libs: Dict = None, single_end_libs: Dict = None, srr_ids: str = None, contigs: str = None, genome_group: str = None, skip_indexing: bool = False, recipe: str = None, viral_recipe: str = None, output_path: str = None, output_file: str = None, force_local_assembly: bool = False, force_inline_annotation: bool = True, perform_bacterial_binning: bool = True, perform_viral_binning: bool = False, perform_viral_annotation: bool = False, perform_bacterial_annotation: bool = True, assembler: str = "", danglen: str = "50", min_contig_len: int = 400, min_contig_cov: float = 4.0) -> str:
    user_id = extract_userid_from_token(token)
    return start_metagenomic_binning_app(api, token=token, user_id=user_id, paired_end_libs=paired_end_libs, single_end_libs=single_end_libs, srr_ids=srr_ids, contigs=contigs, genome_group=genome_group, skip_indexing=skip_indexing, recipe=recipe, viral_recipe=viral_recipe, output_path=output_path, output_file=output_file, force_local_assembly=force_local_assembly, force_inline_annotation=force_inline_annotation, perform_bacterial_binning=perform_bacterial_binning, perform_viral_binning=perform_viral_binning, perform_viral_annotation=perform_viral_annotation, perform_bacterial_annotation=perform_bacterial_annotation, assembler=assembler, danglen=danglen, min_contig_len=min_contig_len, min_contig_cov=min_contig_cov)

@tool(name="service_start_metagenomic_read_mapping_app", 
      description="""Map metagenomic reads to a defined gene set for functional analysis.
      
      This tool performs comprehensive mapping of metagenomic sequencing reads
      to predefined gene sets for functional analysis and gene abundance quantification.
      It supports multiple gene set types including virulence factors, antibiotic
      resistance genes, and custom gene collections.
      
      Key Features:
      - Multiple gene set types: predefined lists, FASTA files, feature groups
      - Predefined gene databases: VFDB (virulence factors), CARD (antibiotic resistance)
      - Flexible input formats: paired-end, single-end, and SRA datasets
      - Read mapping and abundance quantification
      - Functional gene analysis
      - Custom gene set support
      
      Gene Set Options:
      - Predefined lists: VFDB, CARD databases
      - FASTA files: Custom protein sequence collections
      - Feature groups: Annotated gene collections
      - Custom gene sets: User-defined gene collections
      
      Analysis Capabilities:
      - Read mapping to gene sequences
      - Gene abundance quantification
      - Functional gene identification
      - Resistance gene detection
      - Virulence factor analysis
      
      Parameters:
      - gene_set_type: Input type - predefined_list/fasta_file/feature_group (required)
      - gene_set_name: Predefined gene set - VFDB/CARD/feature_group/fasta_file (optional)
      - gene_set_fasta: Gene set FASTA data (optional)
      - gene_set_feature_group: Gene set feature group name (optional)
      - paired_end_libs: Paired-end read libraries (optional)
      - single_end_libs: Single-end read libraries (optional)
      - srr_ids: SRA Run IDs (optional)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_metagenomic_read_mapping_app(token: str = None, gene_set_type: str = None, gene_set_name: str = None, gene_set_fasta: str = None, gene_set_feature_group: str = None, paired_end_libs: Dict = None, single_end_libs: Dict = None, srr_ids: str = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_metagenomic_read_mapping_app(api, token=token, user_id=user_id, gene_set_type=gene_set_type, gene_set_name=gene_set_name, gene_set_fasta=gene_set_fasta, gene_set_feature_group=gene_set_feature_group, paired_end_libs=paired_end_libs, single_end_libs=single_end_libs, srr_ids=srr_ids, output_path=output_path, output_file=output_file)

""" Transcriptomics Services """

@tool(name="service_start_rnaseq_app", 
      description="""Align or assemble RNA-seq reads into transcripts with normalized expression levels.
      
      This tool performs comprehensive RNA-seq analysis including read alignment,
      transcript assembly, and differential expression analysis. It supports both
      bacterial and host genome analysis with multiple analysis recipes and
      experimental design support for comparative studies.
      
      Key Features:
      - Multiple analysis recipes: HTSeq-DESeq, Cufflinks, Host-specific
      - Support for bacterial and host genomes
      - Strand-specific and non-strand-specific analysis
      - Multiple input formats: paired-end, single-end, and SRA datasets
      - Experimental condition support for comparative analysis
      - Read trimming and quality control
      - Differential expression analysis
      
      Analysis Recipes:
      - HTSeq-DESeq: HTSeq counting with DESeq2 differential expression
      - Cufflinks: Transcript assembly and quantification
      - Host: Host-specific RNA-seq analysis pipeline
      
      Experimental Design:
      - Experimental conditions: Define sample groups
      - Contrasts: Statistical comparisons between conditions
      - Sample management: Sample ID tracking and organization
      - Condition indexing: Group samples by experimental conditions
      
      Quality Control:
      - Read trimming with TrimGalore
      - Strand-specific analysis support
      - Insert size estimation and validation
      - Read quality assessment
      
      Parameters:
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
      - skip_sampling: Skip sampling step in alignment (optional)""")
def service_start_rnaseq_app(token: str = None, experimental_conditions: List[str] = None, contrasts: str = None, strand_specific: bool = True, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_libs: List[Dict] = None, reference_genome_id: str = None, genome_type: str = None, recipe: str = "HTSeq-DESeq", host_ftp: str = None, output_path: str = None, output_file: str = None, trimming: bool = False, unit_test: str = None, skip_sampling: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_rnaseq_app(api, token=token, user_id=user_id, experimental_conditions=experimental_conditions, contrasts=contrasts, strand_specific=strand_specific, paired_end_libs=paired_end_libs, single_end_libs=single_end_libs, srr_libs=srr_libs, reference_genome_id=reference_genome_id, genome_type=genome_type, recipe=recipe, host_ftp=host_ftp, output_path=output_path, output_file=output_file, trimming=trimming, unit_test=unit_test, skip_sampling=skip_sampling)

@tool(name="service_start_expression_import_app", 
      description="""Parse and transform user differential expression data.
      
      This tool processes and transforms differential expression data from users,
      supporting various expression data formats and metadata integration. It
      provides data validation, format conversion, and integration capabilities
      for downstream analysis workflows.
      
      Key Features:
      - Expression data parsing and validation
      - Metadata integration and template support
      - Data format transformation and standardization
      - User information processing and validation
      - Output format optimization for analysis
      - Data quality assessment and reporting
      
      Data Processing:
      - Expression data file parsing and validation
      - Metadata template processing and integration
      - User information extraction and validation
      - Data format standardization
      - Quality control and error checking
      
      Input Support:
      - Expression data files: Comparison values between samples
      - Metadata files: User-filled metadata templates
      - User information: JSON string with user details
      - Flexible output options
      
      Parameters:
      - xfile: Experiment data file with comparison values (required)
      - mfile: Metadata template file (optional)
      - ustring: User information JSON string (required)
      - output_path: Output directory (optional)
      - output_file: Output basename (optional)""")
def service_start_expression_import_app(token: str = None, xfile: str = None, mfile: str = None, ustring: str = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_expression_import_app(api, token=token, user_id=user_id, xfile=xfile, mfile=mfile, ustring=ustring, output_path=output_path, output_file=output_file)

""" Utilities Services """

@tool(name="service_start_fastqutils_app", 
      description="""Compute multiple sequence alignment and analyze SNP/variance.
      
      This tool performs comprehensive multiple sequence alignment (MSA) analysis
      to identify single nucleotide polymorphisms (SNPs) and other variations
      across multiple sequences. It supports various input formats and alignment
      algorithms for robust sequence comparison and variant analysis.
      
      Key Features:
      - Multiple sequence alignment using MUSCLE, MAFFT, or progressiveMauve
      - Support for DNA and protein sequences
      - Multiple input formats: FASTA files, genome groups, feature groups, ID lists
      - Reference sequence support for alignment optimization
      - SNP and variation analysis
      - Flexible alignment strategies and parameters
      
      Input Options:
      - FASTA files: Direct sequence file input
      - Genome groups: Pre-defined genome collections
      - Feature groups: Annotated feature collections
      - Feature lists: Direct feature ID collections
      - Genome lists: Direct genome ID collections
      - Keyboard input: Direct FASTA sequence input
      
      Alignment Algorithms:
      - MUSCLE: Multiple sequence alignment by log-expectation
      - MAFFT: Multiple alignment using fast Fourier transform
      - progressiveMauve: Progressive multiple alignment
      
      MAFFT Strategies:
      - auto: Automatic strategy selection
      - fftns1, fftns2, fftnsi: FFT-based methods
      - einsi, linsi, ginsi: Iterative refinement methods
      
      Reference Options:
      - none: No reference sequence
      - string: Direct sequence input
      - feature_id: Feature identifier
      - genome_id: Genome identifier
      - first: Use first sequence as reference
      
      Parameters:
      - input_status: Input alignment status - unaligned/aligned (default: unaligned)
      - input_type: Input type - input_group/input_fasta/input_sequence/input_genomegroup/input_featurelist/input_genomelist (default: input_group)
      - fasta_files: FASTA file groups with file and type (optional)
      - select_genomegroup: Genome groups for analysis (optional)
      - feature_groups: Feature groups for analysis (optional)
      - feature_list: Feature ID lists (optional)
      - genome_list: Genome ID lists (optional)
      - aligner: Alignment algorithm - Muscle/Mafft/progressiveMauve (default: Muscle)
      - alphabet: Sequence alphabet - dna/protein (default: dna)
      - fasta_keyboard_input: Direct FASTA sequence input (optional)
      - ref_type: Reference type - none/string/feature_id/genome_id/first (default: none)
      - strategy: MAFFT strategy - auto/fftns1/fftns2/fftnsi/einsi/linsi/ginsi (default: auto)
      - ref_string: Reference sequence identity (optional)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_fastqutils_app(token: str = None, input_status: str = "unaligned", input_type: str = "input_group", fasta_files: List[Dict] = None, select_genomegroup: List[str] = None, feature_groups: List[str] = None, feature_list: List[str] = None, genome_list: List[str] = None, aligner: str = "Muscle", alphabet: str = "dna", fasta_keyboard_input: str = "", ref_type: str = "none", strategy: str = "auto", ref_string: str = "", output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_fastqutils_app(api, token=token, user_id=user_id, input_status=input_status, input_type=input_type, fasta_files=fasta_files, select_genomegroup=select_genomegroup, feature_groups=feature_groups, feature_list=feature_list, genome_list=genome_list, aligner=aligner, alphabet=alphabet, fasta_keyboard_input=fasta_keyboard_input, ref_type=ref_type, strategy=strategy, ref_string=ref_string, output_path=output_path, output_file=output_file)

@tool(name="service_start_id_mapper_app", description="")
def service_start_id_mapper_app(token: str = None) -> str:
    return null

""" Viral Services """

@tool(name="service_start_sars_genome_analysis_app", 
      description="""Assemble SARS-CoV-2 reads into consensus sequences using specialized assembly strategies.
      
      This tool performs comprehensive SARS-CoV-2 genome assembly from sequencing reads
      using specialized assembly strategies optimized for viral genomes. It supports
      multiple sequencing platforms and primer sets for accurate consensus sequence
      generation and variant detection.
      
      Key Features:
      - SARS-CoV-2 specific assembly strategies: auto, onecodex, cdc-illumina, cdc-nanopore, artic-nanopore
      - Multiple primer sets: ARTIC, midnight, qiagen, swift, varskip, varskip-long
      - Support for various sequencing platforms: Illumina, PacBio, Nanopore, IonTorrent
      - Paired-end and single-end read support
      - SRA dataset integration
      - Coverage depth filtering and consensus generation
      - Intermediate output management
      
      Assembly Recipes:
      - auto: Automatic strategy selection based on input data
      - onecodex: OneCodex assembly method
      - cdc-illumina: CDC Illumina-specific assembly
      - cdc-nanopore: CDC Nanopore-specific assembly
      - artic-nanopore: ARTIC Nanopore assembly
      
      Primer Sets:
      - ARTIC: Standard ARTIC primer set
      - midnight: Midnight primer set
      - qiagen: QIAGEN primer set
      - swift: Swift primer set
      - varskip: VarSkip primer set
      - varskip-long: VarSkip long primer set
      
      Input Options:
      - Paired-end reads: Dual read libraries with platform detection
      - Single-end reads: Single read libraries with platform detection
      - SRA datasets: Sequence Read Archive integration
      - Platform inference: Automatic platform detection
      
      Quality Control:
      - Minimum depth filtering: Threshold for consensus inclusion
      - Maximum depth filtering: Upper limit for read consideration
      - Coverage optimization: Balanced depth for consensus quality
      - Intermediate file management: Optional retention of intermediate outputs
      
      Parameters:
      - paired_end_libs: Paired-end read libraries with read1, read2, platform, interleaved, read_orientation_outward (optional)
      - single_end_libs: Single-end read libraries with read, platform (optional)
      - srr_ids: SRA Run IDs (optional)
      - recipe: Assembly recipe - auto/onecodex/cdc-illumina/cdc-nanopore/artic-nanopore (default: auto)
      - primers: Primer set - ARTIC/midnight/qiagen/swift/varskip/varskip-long (default: ARTIC)
      - primer_version: Version number for primer set (optional)
      - min_depth: Minimum coverage depth for consensus (default: 100)
      - max_depth: Maximum coverage depth for consensus (default: 8000)
      - keep_intermediates: Keep intermediate outputs (default: 0)
      - output_path: Output directory (required)
      - output_file: Output basename (required)
      - debug_level: Debugging level (default: 0)""")
def service_start_sars_genome_analysis_app(token: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_ids: List[str] = None, recipe: str = "auto", primers: str = "ARTIC", primer_version: str = None, min_depth: int = 100, max_depth: int = 8000, keep_intermediates: int = 0, output_path: str = None, output_file: str = None, debug_level: int = 0) -> str:
    user_id = extract_userid_from_token(token)
    return start_sars_genome_analysis_app(api, token=token, user_id=user_id, paired_end_libs=paired_end_libs, single_end_libs=single_end_libs, srr_ids=srr_ids, recipe=recipe, primers=primers, primer_version=primer_version, min_depth=min_depth, max_depth=max_depth, keep_intermediates=keep_intermediates, output_path=output_path, output_file=output_file, debug_level=debug_level)

@tool(name="service_start_sars_wastewater_analysis_app", 
      description="""Assemble SARS-CoV-2 reads into consensus sequences for wastewater surveillance.
      
      This tool performs comprehensive SARS-CoV-2 analysis from wastewater samples,
      including read assembly, variant detection, and lineage analysis. It uses
      specialized primer sets and Freyja analysis for accurate variant identification
      and abundance quantification in complex environmental samples.
      
      Key Features:
      - SARS-CoV-2 consensus sequence assembly
      - Multiple primer set support: ARTIC, Midnight, QIAGEN, Swift, VarSkip
      - Variant detection and lineage analysis using Freyja
      - Multiple sequencing platforms: Illumina, PacBio, Nanopore, IonTorrent
      - Time series analysis for wastewater surveillance
      - Coverage and quality control parameters
      - Sample metadata integration
      
      Primer Sets:
      - ARTIC: Standard ARTIC primer set
      - Midnight: Midnight primer set
      - QIAGEN: QIAGEN primer set
      - Swift: Swift primer set
      - VarSkip: VarSkip primer set (standard and long)
      
      Analysis Pipeline:
      - Read preprocessing and quality control
      - Consensus sequence assembly
      - Variant calling and lineage identification
      - Abundance quantification
      - Time series analysis and visualization
      - Coverage assessment and quality metrics
      
      Quality Control:
      - Minimum base quality score filtering
      - Genome coverage requirements
      - Lineage abundance thresholds
      - Coverage depth filtering
      - Confirmed lineage filtering
      
      Parameters:
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
      """)
def service_start_sars_wastewater_analysis_app(token: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_libs: List[Dict] = None, recipe: str = "auto", primers: str = "ARTIC", minimum_base_quality_score: int = 20, minimum_genome_coverage: int = 60, agg_minimum_lineage_abundance: float = 0.01, minimum_coverage_depth: int = 0, confirmedonly: bool = False, minimum_lineage_abundance: float = 0.001, coverage_estimate: int = 10, timeseries_plot_interval: str = "0", primer_version: str = None, barcode_csv: str = None, sample_metadata_csv: str = None, keep_intermediates: bool = True, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_sars_wastewater_analysis_app(api, token=token, user_id=user_id, paired_end_libs=paired_end_libs, single_end_libs=single_end_libs, srr_libs=srr_libs, recipe=recipe, primers=primers, minimum_base_quality_score=minimum_base_quality_score, minimum_genome_coverage=minimum_genome_coverage, agg_minimum_lineage_abundance=agg_minimum_lineage_abundance, minimum_coverage_depth=minimum_coverage_depth, confirmedonly=confirmedonly, minimum_lineage_abundance=minimum_lineage_abundance, coverage_estimate=coverage_estimate, timeseries_plot_interval=timeseries_plot_interval, primer_version=primer_version, barcode_csv=barcode_csv, sample_metadata_csv=sample_metadata_csv, keep_intermediates=keep_intermediates, output_path=output_path, output_file=output_file)

@tool(name="service_start_sequence_submission_app", 
      description="""Submit sequences to public databases with comprehensive metadata.
      
      This tool facilitates the submission of biological sequences to public
      databases with complete metadata integration. It supports multiple input
      formats and comprehensive submitter information for proper attribution
      and database integration.
      
      Key Features:
      - Multiple input formats: FASTA data, files, genome groups, ID lists
      - Comprehensive metadata integration
      - Submitter information management
      - Affiliation and contact details
      - Geographic location tracking
      - Consortium and organizational support
      - Sequence count validation
      
      Input Options:
      - FASTA data: Direct sequence input
      - FASTA files: Workspace file references
      - Genome groups: Pre-defined genome collections
      - ID lists: Sequence identifier collections
      - Metadata files: CSV format metadata
      
      Submitter Information:
      - Personal details: Name, email, phone
      - Affiliation: Organization and consortium
      - Location: Country, state, city, postal code
      - Contact: Street address and phone number
      
      Parameters:
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
      - output_file: Output basename (required)""")
def service_start_sequence_submission_app(token: str = None, input_source: str = None, input_fasta_data: str = None, input_fasta_file: str = None, input_genome_group: str = None, metadata: str = None, affiliation: str = None, first_name: str = None, last_name: str = None, email: str = None, consortium: str = None, country: str = None, phoneNumber: str = None, street: str = None, postal_code: str = None, city: str = None, state: str = None, numberOfSequences: str = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_sequence_submission_app(api, token=token, user_id=user_id, input_source=input_source, input_fasta_data=input_fasta_data, input_fasta_file=input_fasta_file, input_genome_group=input_genome_group, metadata=metadata, affiliation=affiliation, first_name=first_name, last_name=last_name, email=email, consortium=consortium, country=country, phoneNumber=phoneNumber, street=street, postal_code=postal_code, city=city, state=state, numberOfSequences=numberOfSequences, output_path=output_path, output_file=output_file)

@tool(name="service_start_influenza_ha_subtype_conversion_app", 
      description="""Convert HA subtype numbering for influenza hemagglutinin sequences.
      
      This tool performs HA (hemagglutinin) subtype numbering conversion for
      influenza virus sequences, supporting multiple input formats and various
      numbering schemes. It facilitates standardized sequence analysis and
      comparison across different influenza subtypes.
      
      Key Features:
      - Multiple input formats: feature lists, FASTA data, files, feature groups
      - HA subtype numbering conversion
      - Multiple numbering scheme support
      - Influenza sequence standardization
      - Feature group integration
      - Flexible output options
      
      Input Options:
      - Feature lists: Direct feature ID collections
      - FASTA data: Sequence data in FASTA format
      - FASTA files: Workspace file references
      - Feature groups: Annotated feature collections
      
      Conversion Capabilities:
      - HA subtype standardization
      - Numbering scheme conversion
      - Sequence format optimization
      - Feature annotation integration
      
      Parameters:
      - input_source: Input source type - feature_list/fasta_data/fasta_file/feature_group (required)
      - input_fasta_data: Input sequences in FASTA format (optional)
      - input_fasta_file: Input FASTA file (optional)
      - input_feature_group: Input feature group (optional)
      - input_feature_list: Input feature ID list (optional)
      - types: Selected types in the submission (required)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_influenza_ha_subtype_conversion_app(token: str = None, input_source: str = None, input_fasta_data: str = None, input_fasta_file: str = None, input_feature_group: str = None, input_feature_list: str = None, types: str = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_influenza_ha_subtype_conversion_app(api, token=token, user_id=user_id, input_source=input_source, input_fasta_data=input_fasta_data, input_fasta_file=input_fasta_file, input_feature_group=input_feature_group, input_feature_list=input_feature_list, types=types, output_path=output_path, output_file=output_file)

@tool(name="service_start_subspecies_classification_app", 
      description="""Classify viral subspecies using sequence analysis and reference alignments.
      
      This tool performs subspecies classification for viral sequences using
      multiple sequence alignment (MSA) analysis and reference databases. It
      supports various input formats and virus types for accurate taxonomic
      classification and subspecies identification.
      
      Key Features:
      - Multiple input formats: ID lists, FASTA data, files, genome groups
      - Reference MSA support for accurate classification
      - Multiple virus type support
      - Subspecies-level taxonomic classification
      - Sequence alignment and comparison
      - Flexible output options
      
      Input Options:
      - ID lists: Sequence identifier collections
      - FASTA data: Direct sequence input
      - FASTA files: Workspace file references
      - Genome groups: Pre-defined genome collections
      - Reference MSA: Multiple sequence alignment files
      
      Classification Features:
      - Subspecies-level identification
      - Reference sequence alignment
      - Taxonomic classification
      - Virus type-specific analysis
      - Sequence similarity assessment
      
      Parameters:
      - input_source: Input source type - id_list/fasta_data/fasta_file/genome_group (required)
      - input_fasta_data: Input sequences in FASTA format (optional)
      - input_fasta_file: Input FASTA file (optional)
      - input_genome_group: Input genome group (optional)
      - ref_msa_fasta: Reference multiple sequence alignment (optional)
      - virus_type: Virus type for classification (required)
      - output_path: Output directory (required)
      - output_file: Output basename (required)""")
def service_start_subspecies_classification_app(token: str = None, input_source: str = None, input_fasta_data: str = None, input_fasta_file: str = None, input_genome_group: str = None, ref_msa_fasta: str = None, virus_type: str = None, output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_subspecies_classification_app(api, token=token, user_id=user_id, input_source=input_source, input_fasta_data=input_fasta_data, input_fasta_file=input_fasta_file, input_genome_group=input_genome_group, ref_msa_fasta=ref_msa_fasta, virus_type=virus_type, output_path=output_path, output_file=output_file)

@tool(name="service_start_viral_assembly_app", 
      description="""Assemble viral reads into consensus sequences using specialized viral assembly strategies.
      
      This tool performs viral genome assembly from sequencing reads using
      specialized assembly strategies optimized for viral genomes. It supports
      multiple virus types and assembly algorithms designed for viral sequence
      characteristics and genome structures.
      
      Key Features:
      - Specialized viral assembly strategies: auto, IRMA
      - Multiple virus type support: Influenza, Coronavirus, RSV, Ebola, Flu AD
      - Paired-end and single-end read support
      - SRA dataset integration
      - Viral genome size estimation
      - Optimized assembly parameters for viral genomes
      
      Assembly Strategies:
      - Auto: Automatic strategy selection based on input data
      - IRMA: Influenza Reconstruct and Merge Assembler
      
      Virus Modules:
      - FLU: Influenza virus assembly
      - CoV: Coronavirus assembly
      - RSV: Respiratory syncytial virus assembly
      - EBOLA: Ebola virus assembly
      - FLU_AD: Influenza A/D assembly
      
      Input Options:
      - Paired-end reads: Dual read libraries
      - Single-end reads: Single read libraries
      - SRA datasets: Sequence Read Archive integration
      - Multiple file formats: FASTQ, FASTA, H5
      
      Parameters:
      - paired_end_lib: Paired-end read library (optional)
      - single_end_lib: Single-end read library (optional)
      - srr_id: SRA Run ID (optional)
      - recipe: Assembly strategy - auto/irma (default: auto)
      - module: Virus module - FLU/CoV/RSV/EBOLA/FLU_AD (optional)
      - viral_size: Estimated viral genome size (default: 5M)
      - output_path: Output directory (required)
      - output_file: Output basename (required)
      """)
def service_start_viral_assembly_app(token: str = None, paired_end_lib: Dict = None, single_end_lib: Dict = None, srr_id: str = None, recipe: str = "auto", module: str = None, viral_size: str = "5M", output_path: str = None, output_file: str = None) -> str:
    user_id = extract_userid_from_token(token)
    return start_viral_assembly_app(api, token=token, user_id=user_id, paired_end_lib=paired_end_lib, single_end_lib=single_end_lib, srr_id=srr_id, recipe=recipe, module=module, viral_size=viral_size, output_path=output_path, output_file=output_file)

""" Query Tasks """

@tool(name="service_query_tasks", description="Query tasks. Parameters: token (str, optional) - Authentication token for API access, task_ids (List[str]) - List of task IDs to query")
def service_query_tasks(token: str = None, task_ids: List[str] = None) -> str:
    user_id = extract_userid_from_token(token)
    params = {"task_ids": task_ids} if task_ids else None
    return query_tasks(api, token=token, user_id=user_id, params=params)

def main() -> int:
    print(f"Starting BVBRC Service MCP Flask Server on port {port}...", file=sys.stderr)
    try:
        app.run(host="127.0.0.1", port=port, debug=True)
    except KeyboardInterrupt:
        print("Server stopped.", file=sys.stderr)
    except Exception as e:
        print(f"Server error: {e}", file=sys.stderr)
        return 1
    
    return 0


if __name__ == "__main__":
    main()
