"""
BVBRC Service MCP Tools

This module contains all the MCP tool functions for the BVBRC Service MCP Server.
All tools are registered with the FastMCP server instance.
"""

from fastmcp import FastMCP
from json_rpc import JsonRpcCaller
from service_functions import (
    enumerate_apps, start_date_app, start_genome_annotation_app, query_tasks,
    start_genome_assembly_app, start_comprehensive_genome_analysis_app, start_blast_app,
    start_primer_design_app, start_variation_app, start_tnseq_app, start_bacterial_genome_tree_app,
    start_gene_tree_app, start_core_genome_mlst_app, start_whole_genome_snp_app,
    start_taxonomic_classification_app, start_metagenomic_binning_app, start_metagenomic_read_mapping_app,
    start_rnaseq_app, start_expression_import_app, start_sars_wastewater_analysis_app,
    start_sequence_submission_app, start_influenza_ha_subtype_conversion_app,
    start_subspecies_classification_app, start_viral_assembly_app, start_fastq_utils_app,
    start_genome_alignment_app, start_sars_genome_analysis_app, start_msa_snp_analysis_app,
    start_metacats_app, start_proteome_comparison_app, start_comparative_systems_app,
    start_docking_app, start_similar_genome_finder_app
)
from typing import Any, List, Dict, Optional
from service_verifier import (
    verify_date_app_parameters, verify_genome_assembly_app_parameters, verify_genome_annotation_app_parameters,
    verify_comprehensive_genome_analysis_app_parameters, verify_blast_app_parameters, verify_primer_design_app_parameters,
    verify_variation_app_parameters, verify_tnseq_app_parameters, verify_bacterial_genome_tree_app_parameters,
    verify_gene_tree_app_parameters, verify_core_genome_mlst_app_parameters, verify_whole_genome_snp_app_parameters,
    verify_taxonomic_classification_app_parameters, verify_metagenomic_binning_app_parameters, verify_metagenomic_read_mapping_app_parameters,
    verify_rnaseq_app_parameters, verify_expression_import_app_parameters, verify_sars_wastewater_analysis_app_parameters,
    verify_sequence_submission_app_parameters, verify_influenza_ha_subtype_conversion_app_parameters, verify_subspecies_classification_app_parameters,
    verify_viral_assembly_app_parameters, verify_genome_alignment_app_parameters, verify_sars_genome_analysis_app_parameters,
    verify_msa_snp_analysis_app_parameters, verify_metacats_app_parameters, verify_proteome_comparison_app_parameters,
    verify_comparative_systems_app_parameters, verify_docking_app_parameters, verify_similar_genome_finder_app_parameters,
    verify_fastqutils_app_parameters,
)

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


def register_all_tools(mcp: FastMCP, api: JsonRpcCaller, similar_genome_finder_api: JsonRpcCaller, token_provider):
    """
    Register all MCP tools with the FastMCP server instance.
    
    Args:
        mcp: FastMCP server instance
        api: Main service API caller
        similar_genome_finder_api: Similar genome finder API caller
        token_provider: TokenProvider instance for handling authentication tokens
    """
    
    # Basic Services
    # @mcp.tool(name="enumerate_service_apps", description="Enumerate all BV-BRC Service applications.")
    def service_enumerate_apps() -> List[str]:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return ["Error: No authentication token available"]
        
        user_id = extract_userid_from_token(auth_token)
        return enumerate_apps(api, token=auth_token, user_id=user_id)

    @mcp.tool(name="submit_date_application", description="Start the date app. Parameters defined in the date_app_prompt.")
    def service_start_date_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        
        user_id = extract_userid_from_token(auth_token)
        if not verify_date_app_parameters(params):
            return "Error: Invalid parameters"
        return start_date_app(api, token=auth_token, user_id=user_id, **params)

    # Query Tasks
    @mcp.tool(name="query_running_jobs", description="Query running jobs by IDs. Parameters: task_ids (List[str]) - List of task IDs to query")
    def service_query_tasks(token: Optional[str] = None, task_ids: List[str] = None) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        
        user_id = extract_userid_from_token(auth_token)
        params = {"task_ids": task_ids} if task_ids else None
        return query_tasks(api, token=auth_token, user_id=user_id, params=params)

    # Genomics Analysis Services
    @mcp.tool(name="submit_genome_assembly_application", description="Submit a genome assembly application. Parameters defined in the genome_assembly_app_prompt.")
    def service_start_genome_assembly_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_genome_assembly_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_genome_assembly_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_genome_annotation_application", description="Call genes and functionally annotate input contig set using RASTtk pipeline. Parameters defined in the genome_annotation_app_prompt.")
    def service_start_genome_annotation_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_genome_annotation_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_genome_annotation_app(api, token=auth_token, user_id=user_id, **params)


    @mcp.tool(name="submit_comprehensive_genome_analysis_application", description="Perform comprehensive genome analysis from reads, contigs, or GenBank files. Parameters defined in the comprehensive_genome_analysis_app_prompt.")
    def service_start_comprehensive_genome_analysis_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_comprehensive_genome_analysis_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_comprehensive_genome_analysis_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_blast_application", description="Perform homology searches on sequence data using BLAST algorithms. Parameters defined in the blast_app_prompt.")
    def service_start_blast_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_blast_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_blast_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_primer_design_application", description="Use Primer3 to design primers for given DNA sequences. Parameters defined in the primer_design_app_prompt.")
    def service_start_primer_design_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_primer_design_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_primer_design_app(api, token=auth_token, user_id=user_id, **params)

    # @mcp.tool(name="service_start_variation_app", description="Identify and annotate small nucleotide variations (SNVs) relative to a reference genome")
    @mcp.tool(name="submit_variation_application", description="Identify and annotate small nucleotide variations (SNVs) relative to a reference genome. Parameters defined in the variation_app_prompt.")
    def service_start_variation_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_variation_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_variation_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_tnseq_application", description="Analyze TnSeq (Transposon Sequencing) data using TRANSIT pipeline. Parameters defined in the tnseq_app_prompt.")
    def service_start_tnseq_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_tnseq_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_tnseq_app(api, token=auth_token, user_id=user_id, **params)

# Phylogenomics Services

    @mcp.tool(name="submit_bacterial_genome_tree_application", description="Compute phylogenetic tree from PGFam protein and DNA sequences. Parameters defined in the bacterial_genome_tree_app_prompt.")
    def service_start_bacterial_genome_tree_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_bacterial_genome_tree_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_bacterial_genome_tree_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_gene_tree_application", description="Estimate phylogeny of gene or other sequence features. Parameters defined in the gene_tree_app_prompt.")
    def service_start_gene_tree_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_gene_tree_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_gene_tree_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_core_genome_mlst_application", description="Evaluate core genomes from a set of genome groups of the same species. Parameters defined in the core_genome_mlst_app_prompt.")
    def service_start_core_genome_mlst_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_core_genome_mlst_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_core_genome_mlst_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_whole_genome_snp_application", description="Identify SNP differences in a genome group with genomes of the same species. Parameters defined in the whole_genome_snp_app_prompt.")
    def service_start_whole_genome_snp_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_whole_genome_snp_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_whole_genome_snp_app(api, token=auth_token, user_id=user_id, **params)

# Metagenomics Services

    @mcp.tool(name="submit_taxonomic_classification_application", description="Compute taxonomic classification for read data using metagenomic analysis. Parameters defined in the taxonomic_classification_app_prompt.")
    def service_start_taxonomic_classification_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_taxonomic_classification_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_taxonomic_classification_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_metagenomic_binning_application", description="Assemble, bin, and annotate metagenomic sample data. Parameters defined in the metagenomic_binning_app_prompt.")
    def service_start_metagenomic_binning_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_metagenomic_binning_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_metagenomic_binning_app(api, token=auth_token, user_id=user_id, **params)

    # @mcp.tool()
    @mcp.tool(name="submit_metagenomic_read_mapping_application", description="Map metagenomic reads to a defined gene set for functional analysis. Parameters defined in the metagenomic_read_mapping_app_prompt.")
    def service_start_metagenomic_read_mapping_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_metagenomic_read_mapping_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_metagenomic_read_mapping_app(api, token=auth_token, user_id=user_id, **params)

# Transcriptomics Services

    @mcp.tool(name="submit_rnaseq_application", description="Align or assemble RNA-seq reads into transcripts with normalized expression levels. Parameters defined in the rnaseq_app_prompt.")
    def service_start_rnaseq_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_rnaseq_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_rnaseq_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_expression_import_application", description="Parse and transform user differential expression data. Parameters defined in the expression_import_app_prompt.")
    def service_start_expression_import_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_expression_import_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_expression_import_app(api, token=auth_token, user_id=user_id, **params)

# Viral Services

    @mcp.tool(name="submit_sars_wastewater_analysis_application", description="Assemble SARS-CoV-2 reads into consensus sequences for wastewater surveillance. Parameters defined in the sars_wastewater_analysis_app_prompt.")
    def service_start_sars_wastewater_analysis_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_sars_wastewater_analysis_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_sars_wastewater_analysis_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_sequence_submission_application", description="Submit sequences to public databases with comprehensive metadata. Parameters defined in the sequence_submission_app_prompt.")
    def service_start_sequence_submission_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_sequence_submission_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_sequence_submission_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_influenza_ha_subtype_conversion_application", description="Convert HA subtype numbering for influenza hemagglutinin sequences. Parameters defined in the influenza_ha_subtype_conversion_app_prompt.")
    def service_start_influenza_ha_subtype_conversion_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_influenza_ha_subtype_conversion_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_influenza_ha_subtype_conversion_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_subspecies_classification_application", description="Classify viral subspecies using sequence analysis and reference alignments. Parameters defined in the subspecies_classification_app_prompt.")
    def service_start_subspecies_classification_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_subspecies_classification_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_subspecies_classification_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_viral_assembly_application", description="Assemble viral reads into consensus sequences using specialized viral assembly strategies. Parameters defined in the viral_assembly_app_prompt.")
    def service_start_viral_assembly_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_viral_assembly_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_viral_assembly_app(api, token=auth_token, user_id=user_id, **params)

    # Additional Services

    # @mcp.tool()
    @mcp.tool(name="submit_genome_alignment_application", description="Align multiple genomes using progressiveMauve or other alignment algorithms. Parameters defined in the genome_alignment_app_prompt.")
    def service_start_genome_alignment_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_genome_alignment_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_genome_alignment_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_sars_genome_analysis_application", description="Assemble SARS-CoV-2 genome from sequencing reads. Parameters defined in the sars_genome_analysis_app_prompt.")
    def service_start_sars_genome_analysis_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_sars_genome_analysis_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_sars_genome_analysis_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_msa_snp_analysis_application", description="Analyze SNPs in multiple sequence alignments. Parameters defined in the msa_snp_analysis_app_prompt.")
    def service_start_msa_snp_analysis_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_msa_snp_analysis_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_msa_snp_analysis_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_metacats_application", description="Perform metagenomic analysis of community temporal dynamics. Parameters defined in the metacats_app_prompt.")
    def service_start_metacats_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_metacats_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_metacats_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_proteome_comparison_application", description="Compare proteomes across multiple genomes. Parameters defined in the proteome_comparison_app_prompt.")
    def service_start_proteome_comparison_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_proteome_comparison_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_proteome_comparison_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_comparative_systems_application", description="Compare genomic systems across multiple genomes. Parameters defined in the comparative_systems_app_prompt.")
    def service_start_comparative_systems_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_comparative_systems_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_comparative_systems_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_docking_application", description="Perform molecular docking simulations. Parameters defined in the docking_app_prompt.")
    def service_start_docking_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_docking_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_docking_app(api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_similar_genome_finder_application", description="Find similar genomes using sequence similarity. Parameters defined in the similar_genome_finder_app_prompt.")
    def service_start_similar_genome_finder_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_similar_genome_finder_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_similar_genome_finder_app(similar_genome_finder_api, token=auth_token, user_id=user_id, **params)

    @mcp.tool(name="submit_fastqutils_application", description="Perform quality control and processing on FASTQ sequencing data. Parameters defined in the fastqutils_app_prompt.")
    def service_start_fastqutils_app(params: Dict[str, Any]) -> str:
        # Get the appropriate token
        auth_token = token_provider.get_token()
        if not auth_token:
            return "Error: No authentication token available"
        if not verify_fastqutils_app_parameters(params):
            return "Error: Invalid parameters"
        user_id = extract_userid_from_token(auth_token)
        return start_fastq_utils_app(api, token=auth_token, user_id=user_id, **params)

