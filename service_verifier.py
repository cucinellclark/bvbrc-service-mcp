from typing import Dict, Any

def verify_date_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the date app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False
    return True

def verify_genome_assembly_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the genome assembly app.
    """
    if not params.get("paired_end_libs") or not params.get("single_end_libs") or not params.get("srr_ids"):
        print("Error: paired_end_libs, single_end_libs, or srr_ids is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False
    return True

def verify_genome_annotation_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the genome annotation app.
    """
    if not params.get("contigs"):
        print("Error: contigs is required")
        return False
    if not params.get("scientific_name"):
        print("Error: scientific_name is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False
    return True

def verify_comprehensive_genome_analysis_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the comprehensive genome analysis app.
    """
    if not params.get("input_type"):
        print("Error: input_type is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False
    if not params.get("scientific_name"):
        print("Error: scientific_name is required")
        return False
    if not params.get("code"):
        print("Error: code is required")
        return False
    if not params.get("domain"):
        print("Error: domain is required")
        return False
    
    # Check that at least one input source is provided
    input_sources = [
        params.get("paired_end_libs"),
        params.get("single_end_libs"), 
        params.get("srr_ids"),
        params.get("contigs"),
        params.get("genbank_file"),
        params.get("gto")
    ]
    if not any(input_sources):
        print("Error: At least one input source (paired_end_libs, single_end_libs, srr_ids, contigs, genbank_file, or gto) is required")
        return False
    
    return True

def verify_blast_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the blast app.
    """
    if not params.get("input_type"):
        print("Error: input_type is required")
        return False
    if not params.get("input_source"):
        print("Error: input_source is required")
        return False
    if not params.get("db_type"):
        print("Error: db_type is required")
        return False
    if not params.get("db_source"):
        print("Error: db_source is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one input source is provided
    input_sources = [
        params.get("input_fasta_data"),
        params.get("input_id_list"),
        params.get("input_fasta_file"),
        params.get("input_feature_group"),
        params.get("input_genome_group")
    ]
    if not any(input_sources):
        print("Error: At least one input source (input_fasta_data, input_id_list, input_fasta_file, input_feature_group, or input_genome_group) is required")
        return False

    # Check that at least one database source is provided
    db_sources = [
        params.get("db_fasta_data"),
        params.get("db_fasta_file"),
        params.get("db_id_list"),
        params.get("db_feature_group"),
        params.get("db_genome_group"),
        params.get("db_genome_list"),
        params.get("db_taxon_list"),
        params.get("db_precomputed_database")
    ]
    if not any(db_sources):
        print("Error: At least one database source (db_fasta_data, db_fasta_file, db_id_list, db_feature_group, db_genome_group, db_genome_list, db_taxon_list, or db_precomputed_database) is required")
        return False

    return True

def verify_primer_design_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the primer design app.
    """
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("input_type"):
        print("Error: input_type is required")
        return False
    if not params.get("sequence_input"):
        print("Error: sequence_input is required")
        return False

    return True

def verify_variation_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the variation app.
    """
    if not params.get("reference_genome_id"):
        print("Error: reference_genome_id is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one input source is provided
    input_sources = [
        params.get("paired_end_libs"),
        params.get("single_end_libs"),
        params.get("srr_ids")
    ]
    if not any(input_sources):
        print("Error: At least one input source (paired_end_libs, single_end_libs, or srr_ids) is required")
        return False

    return True

def verify_tnseq_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the tnseq app.
    """
    if not params.get("reference_genome_id"):
        print("Error: reference_genome_id is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    return True

def verify_bacterial_genome_tree_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the bacterial genome tree app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one genome input is provided
    genome_sources = [
        params.get("genome_ids"),
        params.get("genome_groups")
    ]
    if not any(genome_sources):
        print("Error: At least one genome input (genome_ids or genome_groups) is required")
        return False

    return True

def verify_gene_tree_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the gene tree app.
    """
    if not params.get("sequences"):
        print("Error: sequences is required")
        return False
    if not params.get("alphabet"):
        print("Error: alphabet is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    return True

def verify_core_genome_mlst_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the core genome MLST app.
    """
    if not params.get("input_schema_selection"):
        print("Error: input_schema_selection is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one genome input is provided
    genome_sources = [
        params.get("input_genome_group"),
        params.get("input_genome_fasta")
    ]
    if not any(genome_sources):
        print("Error: At least one genome input (input_genome_group or input_genome_fasta) is required")
        return False

    return True

def verify_whole_genome_snp_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the whole genome SNP app.
    """
    if not params.get("input_genome_type"):
        print("Error: input_genome_type is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one genome input is provided
    genome_sources = [
        params.get("input_genome_group"),
        params.get("input_genome_fasta")
    ]
    if not any(genome_sources):
        print("Error: At least one genome input (input_genome_group or input_genome_fasta) is required")
        return False

    return True

def verify_taxonomic_classification_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the taxonomic classification app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one input source is provided
    input_sources = [
        params.get("paired_end_libs"),
        params.get("single_end_libs"),
        params.get("srr_libs")
    ]
    if not any(input_sources):
        print("Error: At least one input source (paired_end_libs, single_end_libs, or srr_libs) is required")
        return False

    return True

def verify_metagenomic_binning_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the metagenomic binning app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one input source is provided
    input_sources = [
        params.get("paired_end_libs"),
        params.get("single_end_libs"),
        params.get("srr_ids"),
        params.get("contigs"),
        params.get("genome_group")
    ]
    if not any(input_sources):
        print("Error: At least one input source (paired_end_libs, single_end_libs, srr_ids, contigs, or genome_group) is required")
        return False

    return True

def verify_metagenomic_read_mapping_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the metagenomic read mapping app.
    """
    if not params.get("gene_set_type"):
        print("Error: gene_set_type is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one gene set source is provided
    gene_sources = [
        params.get("gene_set_name"),
        params.get("gene_set_fasta"),
        params.get("gene_set_feature_group")
    ]
    if not any(gene_sources):
        print("Error: At least one gene set source (gene_set_name, gene_set_fasta, or gene_set_feature_group) is required")
        return False

    # Check that at least one input source is provided
    input_sources = [
        params.get("paired_end_libs"),
        params.get("single_end_libs"),
        params.get("srr_ids")
    ]
    if not any(input_sources):
        print("Error: At least one input source (paired_end_libs, single_end_libs, or srr_ids) is required")
        return False

    return True

def verify_rnaseq_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the RNA-seq app.
    """
    if not params.get("reference_genome_id"):
        print("Error: reference_genome_id is required")
        return False
    if not params.get("genome_type"):
        print("Error: genome_type is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one input source is provided
    input_sources = [
        params.get("paired_end_libs"),
        params.get("single_end_libs"),
        params.get("srr_libs")
    ]
    if not any(input_sources):
        print("Error: At least one input source (paired_end_libs, single_end_libs, or srr_libs) is required")
        return False

    return True

def verify_expression_import_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the expression import app.
    """
    if not params.get("xfile"):
        print("Error: xfile is required")
        return False
    if not params.get("ustring"):
        print("Error: ustring is required")
        return False

    return True

def verify_sars_wastewater_analysis_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the SARS wastewater analysis app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one input source is provided
    input_sources = [
        params.get("paired_end_libs"),
        params.get("single_end_libs"),
        params.get("srr_libs")
    ]
    if not any(input_sources):
        print("Error: At least one input source (paired_end_libs, single_end_libs, or srr_libs) is required")
        return False

    return True

def verify_sequence_submission_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the sequence submission app.
    """
    if not params.get("input_source"):
        print("Error: input_source is required")
        return False
    if not params.get("metadata"):
        print("Error: metadata is required")
        return False
    if not params.get("first_name"):
        print("Error: first_name is required")
        return False
    if not params.get("last_name"):
        print("Error: last_name is required")
        return False
    if not params.get("email"):
        print("Error: email is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one sequence input is provided
    sequence_sources = [
        params.get("input_fasta_data"),
        params.get("input_fasta_file"),
        params.get("input_genome_group")
    ]
    if not any(sequence_sources):
        print("Error: At least one sequence input (input_fasta_data, input_fasta_file, or input_genome_group) is required")
        return False

    return True

def verify_influenza_ha_subtype_conversion_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the influenza HA subtype conversion app.
    """
    if not params.get("input_source"):
        print("Error: input_source is required")
        return False
    if not params.get("types"):
        print("Error: types is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one sequence input is provided
    sequence_sources = [
        params.get("input_fasta_data"),
        params.get("input_fasta_file"),
        params.get("input_feature_group"),
        params.get("input_feature_list")
    ]
    if not any(sequence_sources):
        print("Error: At least one sequence input (input_fasta_data, input_fasta_file, input_feature_group, or input_feature_list) is required")
        return False

    return True

def verify_subspecies_classification_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the subspecies classification app.
    """
    if not params.get("input_source"):
        print("Error: input_source is required")
        return False
    if not params.get("virus_type"):
        print("Error: virus_type is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one sequence input is provided
    sequence_sources = [
        params.get("input_fasta_data"),
        params.get("input_fasta_file"),
        params.get("input_genome_group")
    ]
    if not any(sequence_sources):
        print("Error: At least one sequence input (input_fasta_data, input_fasta_file, or input_genome_group) is required")
        return False

    return True

def verify_viral_assembly_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the viral assembly app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one input source is provided
    input_sources = [
        params.get("paired_end_lib"),
        params.get("single_end_lib"),
        params.get("srr_id")
    ]
    if not any(input_sources):
        print("Error: At least one input source (paired_end_lib, single_end_lib, or srr_id) is required")
        return False

    return True

def verify_genome_alignment_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the genome alignment app.
    """
    if not params.get("genome_ids"):
        print("Error: genome_ids is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    return True

def verify_sars_genome_analysis_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the SARS genome analysis app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one input source is provided
    input_sources = [
        params.get("paired_end_libs"),
        params.get("single_end_libs"),
        params.get("srr_ids")
    ]
    if not any(input_sources):
        print("Error: At least one input source (paired_end_libs, single_end_libs, or srr_ids) is required")
        return False

    return True

def verify_msa_snp_analysis_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the MSA SNP analysis app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one sequence input is provided
    sequence_sources = [
        params.get("fasta_files"),
        params.get("select_genomegroup"),
        params.get("feature_groups"),
        params.get("feature_list"),
        params.get("genome_list"),
        params.get("fasta_keyboard_input")
    ]
    if not any(sequence_sources):
        print("Error: At least one sequence input (fasta_files, select_genomegroup, feature_groups, feature_list, genome_list, or fasta_keyboard_input) is required")
        return False

    return True

def verify_metacats_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the metacats app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    return True

def verify_proteome_comparison_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the proteome comparison app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one genome input is provided
    genome_sources = [
        params.get("genome_ids"),
        params.get("user_genomes"),
        params.get("user_feature_groups")
    ]
    if not any(genome_sources):
        print("Error: At least one genome input (genome_ids, user_genomes, or user_feature_groups) is required")
        return False

    return True

def verify_comparative_systems_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the comparative systems app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one genome input is provided
    genome_sources = [
        params.get("genome_ids"),
        params.get("genome_groups")
    ]
    if not any(genome_sources):
        print("Error: At least one genome input (genome_ids or genome_groups) is required")
        return False

    return True

def verify_docking_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the docking app.
    """
    if not params.get("protein_input_type"):
        print("Error: protein_input_type is required")
        return False
    if not params.get("ligand_library_type"):
        print("Error: ligand_library_type is required")
        return False
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one protein input is provided
    protein_sources = [
        params.get("input_pdb"),
        params.get("user_pdb_file")
    ]
    if not any(protein_sources):
        print("Error: At least one protein input (input_pdb or user_pdb_file) is required")
        return False

    # Check that at least one ligand input is provided
    ligand_sources = [
        params.get("ligand_named_library"),
        params.get("ligand_smiles_list"),
        params.get("ligand_ws_file")
    ]
    if not any(ligand_sources):
        print("Error: At least one ligand input (ligand_named_library, ligand_smiles_list, or ligand_ws_file) is required")
        return False

    return True

def verify_similar_genome_finder_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the similar genome finder app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    # Check that at least one query input is provided
    query_sources = [
        params.get("selectedGenomeId"),
        params.get("fasta_file")
    ]
    if not any(query_sources):
        print("Error: At least one query input (selectedGenomeId or fasta_file) is required")
        return False

    return True

def verify_fastqutils_app_parameters(params: Dict[str, Any]) -> bool:
    """
    Verify the parameters for the fastqutils app.
    """
    if not params.get("output_path"):
        print("Error: output_path is required")
        return False
    if not params.get("output_file"):
        print("Error: output_file is required")
        return False

    return True
