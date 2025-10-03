from json_rpc import JsonRpcCaller
from typing import List, Dict, Any
import uuid

def _set_default_output_paths(user_id: str, app_name: str, output_path: str = None, output_file: str = None):
    """Helper function to set default output_path and output_file values."""
    if output_path is None:
        output_path = '/' + user_id + '/CopilotDevWorkflows'
    if output_file is None:
        output_file = app_name + '_' + str(uuid.uuid4())
    return output_path, output_file

def enumerate_apps(api: JsonRpcCaller, token: str = None, user_id: str = None) -> List[str]:
    try:
        result = api.call("AppService.enumerate_apps", {}, 1, token)
        print(result)
        return result
    except Exception as e:
        print(e)
        return []

def start_date_app(api: JsonRpcCaller, token: str = None, user_id: str = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "Date"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "output_path": output_path,
            "output_file": output_file
        }
        
        data = ["Date", params, {}]
        
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_genome_annotation_app(api: JsonRpcCaller, token: str = None, user_id: str = None, genome_id: str = None, contigs: str = None, scientific_name: str = None, tax_id: str = None, my_label: str = None, reference_genome_id: str = None, taxonomy_id: str = None, code: int = 0, domain: str = "auto", public: bool = False, queue_nowait: bool = False, skip_indexing: bool = False, skip_workspace_output: bool = False, output_path: str = None, output_file: str = None, lowvan_min_contig_length: int = 300, lowvan_max_contig_length: int = 35000, reference_virus_name: str = None, fix_errors: bool = None, fix_frameshifts: bool = None, verbose_level: int = None, workflow: str = None, recipe: str = None, disable_replication: bool = None, analyze_quality: bool = None, assembly_output: str = None, custom_pipeline: Dict = None) -> str:
    app_name = "GenomeAnnotation"
    try:
        if genome_id is None:
            genome_id = ""
        if contigs is None:
            contigs = ""
        if scientific_name is None:
            scientific_name = ""
        if tax_id is None:
            tax_id = ""
        if my_label is None:
            my_label = ""
        if reference_genome_id is None:
            reference_genome_id = ""
        if taxonomy_id is None:
            taxonomy_id = ""
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "genome_id": genome_id,
            "contigs": contigs,
            "scientific_name": scientific_name,
            "tax_id": tax_id,
            "my_label": my_label,
            "reference_genome_id": reference_genome_id,
            "taxonomy_id": taxonomy_id,
            "code": code,
            "domain": domain,
            "public": public,
            "queue_nowait": queue_nowait,
            "skip_indexing": skip_indexing,
            "skip_workspace_output": skip_workspace_output,
            "output_path": output_path,
            "output_file": output_file,
            "lowvan_min_contig_length": lowvan_min_contig_length,
            "lowvan_max_contig_length": lowvan_max_contig_length,
            "reference_virus_name": reference_virus_name,
            "fix_errors": fix_errors,
            "fix_frameshifts": fix_frameshifts,
            "verbose_level": verbose_level,
            "workflow": workflow,
            "recipe": recipe,
            "disable_replication": disable_replication,
            "analyze_quality": analyze_quality,
            "assembly_output": assembly_output,
            "custom_pipeline": custom_pipeline
        }
        data = ["GenomeAnnotation", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def query_tasks(api: JsonRpcCaller, token: str = None, user_id: str = None, params: Dict[str, Any] = None) -> str:
    try:
        result = api.call("AppService.query_tasks", [params['task_ids']], 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_genome_assembly_app(api: JsonRpcCaller, token: str = None, user_id: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_ids: List[str] = None, max_bases: int = 10000000000, recipe: str = "auto", racon_iter: int = 2, pilon_iter: int = 2, trim: bool = False, target_depth: int = 200, normalize: bool = False, filtlong: bool = False, genome_size: int = 5000000, min_contig_len: int = 300, min_contig_cov: float = 5.0, output_path: str = None, output_file: str = None, debug: int = 0) -> str:
    app_name = "GenomeAssembly"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "paired_end_libs": paired_end_libs or [],
            "single_end_libs": single_end_libs or [],
            "srr_ids": srr_ids or [],
            "max_bases": max_bases,
            "recipe": recipe,
            "racon_iter": racon_iter,
            "pilon_iter": pilon_iter,
            "trim": trim,
            "target_depth": target_depth,
            "normalize": normalize,
            "filtlong": filtlong,
            "genome_size": genome_size,
            "min_contig_len": min_contig_len,
            "min_contig_cov": min_contig_cov,
            "output_path": output_path,
            "output_file": output_file,
            "debug": debug
        }
        data = ["GenomeAssembly", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_comprehensive_genome_analysis_app(api: JsonRpcCaller, token: str = None, user_id: str = None, input_type: str = None, output_path: str = None, output_file: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_ids: List[str] = None, reference_assembly: str = None, recipe: str = "auto", racon_iter: int = 2, pilon_iter: int = 2, trim: bool = False, normalize: bool = False, filtlong: bool = False, target_depth: int = 200, genome_size: int = 5000000, min_contig_len: int = 300, min_contig_cov: float = 5.0, gto: str = None, genbank_file: str = None, contigs: str = None, scientific_name: str = None, taxonomy_id: int = None, code: int = 0, domain: str = "auto", public: bool = False, queue_nowait: bool = False, skip_indexing: bool = False, reference_genome_id: str = None, analyze_quality: bool = None, debug_level: int = 0) -> str:
    app_name = "ComprehensiveGenomeAnalysis"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "input_type": input_type or "",
            "output_path": output_path,
            "output_file": output_file,
            "paired_end_libs": paired_end_libs or [],
            "single_end_libs": single_end_libs or [],
            "srr_ids": srr_ids or [],
            "reference_assembly": reference_assembly or "",
            "recipe": recipe,
            "racon_iter": racon_iter,
            "pilon_iter": pilon_iter,
            "trim": trim,
            "normalize": normalize,
            "filtlong": filtlong,
            "target_depth": target_depth,
            "genome_size": genome_size,
            "min_contig_len": min_contig_len,
            "min_contig_cov": min_contig_cov,
            "gto": gto or "",
            "genbank_file": genbank_file or "",
            "contigs": contigs or "",
            "scientific_name": scientific_name or "",
            "taxonomy_id": taxonomy_id or "",
            "code": code,
            "domain": domain,
            "public": public,
            "queue_nowait": queue_nowait,
            "skip_indexing": skip_indexing,
            "reference_genome_id": reference_genome_id or "",
            "analyze_quality": analyze_quality,
            "debug_level": debug_level
        }
        data = ["ComprehensiveGenomeAnalysis", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_blast_app(api: JsonRpcCaller, token: str = None, user_id: str = None, input_type: str = None, input_source: str = None, input_fasta_data: str = None, input_id_list: List[str] = None, input_fasta_file: str = None, input_feature_group: str = None, input_genome_group: str = None, db_type: str = None, db_source: str = None, db_fasta_data: str = None, db_fasta_file: str = None, db_id_list: List[str] = None, db_feature_group: str = None, db_genome_group: str = None, db_genome_list: List[str] = None, db_taxon_list: List[str] = None, db_precomputed_database: str = None, blast_program: str = None, blast_evalue_cutoff: float = 1e-5, blast_max_hits: int = 300, blast_min_coverage: int = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "BLAST"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "input_type": input_type or "",
            "input_source": input_source or "",
            "input_fasta_data": input_fasta_data or "",
            "input_id_list": input_id_list or [],
            "input_fasta_file": input_fasta_file or "",
            "input_feature_group": input_feature_group or "",
            "input_genome_group": input_genome_group or "",
            "db_type": db_type or "",
            "db_source": db_source or "",
            "db_fasta_data": db_fasta_data or "",
            "db_fasta_file": db_fasta_file or "",
            "db_id_list": db_id_list or [],
            "db_feature_group": db_feature_group or "",
            "db_genome_group": db_genome_group or "",
            "db_genome_list": db_genome_list or [],
            "db_taxon_list": db_taxon_list or [],
            "db_precomputed_database": db_precomputed_database or "",
            "blast_program": blast_program or "",
            "blast_evalue_cutoff": blast_evalue_cutoff,
            "blast_max_hits": blast_max_hits,
            "blast_min_coverage": blast_min_coverage,
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["BLAST", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_primer_design_app(api: JsonRpcCaller, token: str = None, user_id: str = None, output_file: str = None, output_path: str = None, input_type: str = None, sequence_input: str = None, SEQUENCE_ID: str = None, SEQUENCE_TARGET: List[List[int]] = None, SEQUENCE_INCLUDED_REGION: List[int] = None, SEQUENCE_EXCLUDED_REGION: List[int] = None, SEQUENCE_OVERLAP_JUNCTION_LIST: List[List[int]] = None, PRIMER_PICK_INTERNAL_OLIGO: int = None, PRIMER_PRODUCT_SIZE_RANGE: List[List[int]] = None, PRIMER_NUM_RETURN: int = None, PRIMER_MIN_SIZE: int = None, PRIMER_OPT_SIZE: int = None, PRIMER_MAX_SIZE: int = None, PRIMER_MAX_TM: float = None, PRIMER_MIN_TM: float = None, PRIMER_OPT_TM: float = None, PRIMER_PAIR_MAX_DIFF_TM: float = None, PRIMER_MAX_GC: float = None, PRIMER_MIN_GC: float = None, PRIMER_OPT_GC: float = None, PRIMER_SALT_MONOVALENT: float = None, PRIMER_SALT_DIVALENT: float = None, PRIMER_DNA_CONC: float = None, PRIMER_DNTP_CONC: float = None) -> str:
    app_name = "PrimerDesign"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "output_file": output_file,
            "output_path": output_path,
            "input_type": input_type or "",
            "sequence_input": sequence_input or "",
            "SEQUENCE_ID": SEQUENCE_ID or "",
            "SEQUENCE_TARGET": SEQUENCE_TARGET or [],
            "SEQUENCE_INCLUDED_REGION": SEQUENCE_INCLUDED_REGION or [],
            "SEQUENCE_EXCLUDED_REGION": SEQUENCE_EXCLUDED_REGION or [],
            "SEQUENCE_OVERLAP_JUNCTION_LIST": SEQUENCE_OVERLAP_JUNCTION_LIST or [],
            "PRIMER_PICK_INTERNAL_OLIGO": PRIMER_PICK_INTERNAL_OLIGO,
            "PRIMER_PRODUCT_SIZE_RANGE": PRIMER_PRODUCT_SIZE_RANGE or [],
            "PRIMER_NUM_RETURN": PRIMER_NUM_RETURN,
            "PRIMER_MIN_SIZE": PRIMER_MIN_SIZE,
            "PRIMER_OPT_SIZE": PRIMER_OPT_SIZE,
            "PRIMER_MAX_SIZE": PRIMER_MAX_SIZE,
            "PRIMER_MAX_TM": PRIMER_MAX_TM,
            "PRIMER_MIN_TM": PRIMER_MIN_TM,
            "PRIMER_OPT_TM": PRIMER_OPT_TM,
            "PRIMER_PAIR_MAX_DIFF_TM": PRIMER_PAIR_MAX_DIFF_TM,
            "PRIMER_MAX_GC": PRIMER_MAX_GC,
            "PRIMER_MIN_GC": PRIMER_MIN_GC,
            "PRIMER_OPT_GC": PRIMER_OPT_GC,
            "PRIMER_SALT_MONOVALENT": PRIMER_SALT_MONOVALENT,
            "PRIMER_SALT_DIVALENT": PRIMER_SALT_DIVALENT,
            "PRIMER_DNA_CONC": PRIMER_DNA_CONC,
            "PRIMER_DNTP_CONC": PRIMER_DNTP_CONC
        }
        data = ["PrimerDesign", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_variation_app(api: JsonRpcCaller, token: str = None, user_id: str = None, reference_genome_id: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_ids: List[str] = None, mapper: str = "BWA-mem", caller: str = "FreeBayes", output_path: str = None, output_file: str = None, debug: bool = False) -> str:
    app_name = "Variation"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "reference_genome_id": reference_genome_id or "",
            "paired_end_libs": paired_end_libs or [],
            "single_end_libs": single_end_libs or [],
            "srr_ids": srr_ids or [],
            "mapper": mapper,
            "caller": caller,
            "output_path": output_path,
            "output_file": output_file,
            "debug": debug
        }
        data = ["Variation", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_tnseq_app(api: JsonRpcCaller, token: str = None, user_id: str = None, experimental_conditions: List[str] = None, contrasts: List[str] = None, read_files: List[Dict] = None, reference_genome_id: str = None, recipe: str = "gumbel", protocol: str = "sassetti", primer: str = "", output_path: str = None, output_file: str = None) -> str:
    app_name = "TnSeq"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "experimental_conditions": experimental_conditions or [],
            "contrasts": contrasts or [],
            "read_files": read_files or [],
            "reference_genome_id": reference_genome_id or "",
            "recipe": recipe,
            "protocol": protocol,
            "primer": primer,
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["TnSeq", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_bacterial_genome_tree_app(api: JsonRpcCaller, token: str = None, user_id: str = None, output_path: str = None, output_file: str = None, genome_ids: List[str] = None, genome_groups: List[str] = None, optional_genome_ids: List[str] = None, genome_metadata_fields: str = None, number_of_genes: int = 20, bootstraps: int = 100, max_genomes_missing: int = 0, max_allowed_dups: int = 0) -> str:
    app_name = "BacterialGenomeTree"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "output_path": output_path,
            "output_file": output_file,
            "genome_ids": genome_ids or [],
            "genome_groups": genome_groups or [],
            "optional_genome_ids": optional_genome_ids or [],
            "genome_metadata_fields": genome_metadata_fields or "",
            "number_of_genes": number_of_genes,
            "bootstraps": bootstraps,
            "max_genomes_missing": max_genomes_missing,
            "max_allowed_dups": max_allowed_dups
        }
        data = ["BacterialGenomeTree", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_gene_tree_app(api: JsonRpcCaller, token: str = None, user_id: str = None, sequences: List[str] = None, alignment_program: str = None, trim_threshold: float = None, gap_threshold: float = None, alphabet: str = None, substitution_model: str = None, bootstrap: int = None, recipe: str = "RAxML", tree_type: str = None, feature_metadata_fields: str = None, genome_metadata_fields: str = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "GeneTree"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "sequences": sequences or [],
            "alignment_program": alignment_program or "",
            "trim_threshold": trim_threshold,
            "gap_threshold": gap_threshold,
            "alphabet": alphabet or "",
            "substitution_model": substitution_model or "",
            "bootstrap": bootstrap,
            "recipe": recipe,
            "tree_type": tree_type or "",
            "feature_metadata_fields": feature_metadata_fields or "",
            "genome_metadata_fields": genome_metadata_fields or "",
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["GeneTree", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_core_genome_mlst_app(api: JsonRpcCaller, token: str = None, user_id: str = None, input_genome_type: str = "genome_group", analysis_type: str = "chewbbaca", input_genome_group: str = None, input_genome_fasta: str = None, schema_location: str = None, input_schema_selection: str = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "CoreGenomeMLST"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "input_genome_type": input_genome_type,
            "analysis_type": analysis_type,
            "input_genome_group": input_genome_group or "",
            "input_genome_fasta": input_genome_fasta or "",
            "schema_location": schema_location or "",
            "input_schema_selection": input_schema_selection or "",
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["CoreGenomeMLST", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_whole_genome_snp_app(api: JsonRpcCaller, token: str = None, user_id: str = None, input_genome_type: str = None, majority_threshold: float = 0.5, min_mid_linkage: int = 10, max_mid_linkage: int = 40, analysis_type: str = "Whole Genome SNP Analysis", input_genome_group: str = None, input_genome_fasta: str = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "WholeGenomeSNP"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "input_genome_type": input_genome_type or "",
            "majority_threshold": majority_threshold,
            "min_mid_linkage": min_mid_linkage,
            "max_mid_linkage": max_mid_linkage,
            "analysis_type": analysis_type,
            "input_genome_group": input_genome_group or "",
            "input_genome_fasta": input_genome_fasta or "",
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["WholeGenomeSNP", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_taxonomic_classification_app(api: JsonRpcCaller, token: str = None, user_id: str = None, host_genome: str = "no_host", analysis_type: str = "16S", paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_libs: List[Dict] = None, database: str = "SILVA", save_classified_sequences: bool = False, save_unclassified_sequences: bool = False, confidence_interval: float = 0.1, output_path: str = None, output_file: str = None) -> str:
    app_name = "TaxonomicClassification"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "host_genome": host_genome,
            "analysis_type": analysis_type,
            "paired_end_libs": paired_end_libs or [],
            "single_end_libs": single_end_libs or [],
            "srr_libs": srr_libs or [],
            "database": database,
            "save_classified_sequences": save_classified_sequences,
            "save_unclassified_sequences": save_unclassified_sequences,
            "confidence_interval": confidence_interval,
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["TaxonomicClassification", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_metagenomic_binning_app(api: JsonRpcCaller, token: str = None, user_id: str = None, paired_end_libs: Dict = None, single_end_libs: Dict = None, srr_ids: str = None, contigs: str = None, genome_group: str = None, skip_indexing: bool = False, recipe: str = None, viral_recipe: str = None, output_path: str = None, output_file: str = None, force_local_assembly: bool = False, force_inline_annotation: bool = True, perform_bacterial_binning: bool = True, perform_viral_binning: bool = False, perform_viral_annotation: bool = False, perform_bacterial_annotation: bool = True, assembler: str = "", danglen: str = "50", min_contig_len: int = 400, min_contig_cov: float = 4.0) -> str:
    app_name = "MetagenomicBinning"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "paired_end_libs": paired_end_libs or {},
            "single_end_libs": single_end_libs or {},
            "srr_ids": srr_ids or "",
            "contigs": contigs or "",
            "genome_group": genome_group or "",
            "skip_indexing": skip_indexing,
            "recipe": recipe or "",
            "viral_recipe": viral_recipe or "",
            "output_path": output_path,
            "output_file": output_file,
            "force_local_assembly": force_local_assembly,
            "force_inline_annotation": force_inline_annotation,
            "perform_bacterial_binning": perform_bacterial_binning,
            "perform_viral_binning": perform_viral_binning,
            "perform_viral_annotation": perform_viral_annotation,
            "perform_bacterial_annotation": perform_bacterial_annotation,
            "assembler": assembler,
            "danglen": danglen,
            "min_contig_len": min_contig_len,
            "min_contig_cov": min_contig_cov
        }
        data = ["MetagenomicBinning", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_metagenomic_read_mapping_app(api: JsonRpcCaller, token: str = None, user_id: str = None, gene_set_type: str = None, gene_set_name: str = None, gene_set_fasta: str = None, gene_set_feature_group: str = None, paired_end_libs: Dict = None, single_end_libs: Dict = None, srr_ids: str = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "MetagenomicReadMapping"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "gene_set_type": gene_set_type or "",
            "gene_set_name": gene_set_name or "",
            "gene_set_fasta": gene_set_fasta or "",
            "gene_set_feature_group": gene_set_feature_group or "",
            "paired_end_libs": paired_end_libs or {},
            "single_end_libs": single_end_libs or {},
            "srr_ids": srr_ids or "",
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["MetagenomicReadMapping", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_rnaseq_app(api: JsonRpcCaller, token: str = None, user_id: str = None, experimental_conditions: List[str] = None, contrasts: str = None, strand_specific: bool = True, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_libs: List[Dict] = None, reference_genome_id: str = None, genome_type: str = None, recipe: str = "HTSeq-DESeq", host_ftp: str = None, output_path: str = None, output_file: str = None, trimming: bool = False, unit_test: str = None, skip_sampling: str = None) -> str:
    app_name = "RNAseq"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "experimental_conditions": experimental_conditions or [],
            "contrasts": contrasts or "",
            "strand_specific": strand_specific,
            "paired_end_libs": paired_end_libs or [],
            "single_end_libs": single_end_libs or [],
            "srr_libs": srr_libs or [],
            "reference_genome_id": reference_genome_id or "",
            "genome_type": genome_type or "",
            "recipe": recipe,
            "host_ftp": host_ftp or "",
            "output_path": output_path,
            "output_file": output_file,
            "trimming": trimming,
            "unit_test": unit_test or "",
            "skip_sampling": skip_sampling or ""
        }
        data = ["RNAseq", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_expression_import_app(api: JsonRpcCaller, token: str = None, user_id: str = None, xfile: str = None, mfile: str = None, ustring: str = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "ExpressionImport"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "xfile": xfile or "",
            "mfile": mfile or "",
            "ustring": ustring or "",
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["ExpressionImport", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_sars_wastewater_analysis_app(api: JsonRpcCaller, token: str = None, user_id: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_libs: List[Dict] = None, recipe: str = "auto", primers: str = "ARTIC", minimum_base_quality_score: int = 20, minimum_genome_coverage: int = 60, agg_minimum_lineage_abundance: float = 0.01, minimum_coverage_depth: int = 0, confirmedonly: bool = False, minimum_lineage_abundance: float = 0.001, coverage_estimate: int = 10, timeseries_plot_interval: str = "0", primer_version: str = None, barcode_csv: str = None, sample_metadata_csv: str = None, keep_intermediates: bool = True, output_path: str = None, output_file: str = None, debug_level: int = 0) -> str:
    app_name = "SARSWastewaterAnalysis"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "paired_end_libs": paired_end_libs or [],
            "single_end_libs": single_end_libs or [],
            "srr_libs": srr_libs or [],
            "recipe": recipe,
            "primers": primers,
            "minimum_base_quality_score": minimum_base_quality_score,
            "minimum_genome_coverage": minimum_genome_coverage,
            "agg_minimum_lineage_abundance": agg_minimum_lineage_abundance,
            "minimum_coverage_depth": minimum_coverage_depth,
            "confirmedonly": confirmedonly,
            "minimum_lineage_abundance": minimum_lineage_abundance,
            "coverage_estimate": coverage_estimate,
            "timeseries_plot_interval": timeseries_plot_interval,
            "primer_version": primer_version or "",
            "barcode_csv": barcode_csv or "",
            "sample_metadata_csv": sample_metadata_csv or "",
            "keep_intermediates": keep_intermediates,
            "output_path": output_path,
            "output_file": output_file,
            "debug_level": debug_level
        }
        data = ["SARSWastewaterAnalysis", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_sequence_submission_app(api: JsonRpcCaller, token: str = None, user_id: str = None, input_source: str = None, input_fasta_data: str = None, input_fasta_file: str = None, input_genome_group: str = None, metadata: str = None, affiliation: str = None, first_name: str = None, last_name: str = None, email: str = None, consortium: str = None, country: str = None, phoneNumber: str = None, street: str = None, postal_code: str = None, city: str = None, state: str = None, numberOfSequences: str = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "SequenceSubmission"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "input_source": input_source or "",
            "input_fasta_data": input_fasta_data or "",
            "input_fasta_file": input_fasta_file or "",
            "input_genome_group": input_genome_group or "",
            "metadata": metadata or "",
            "affiliation": affiliation or "",
            "first_name": first_name or "",
            "last_name": last_name or "",
            "email": email or "",
            "consortium": consortium or "",
            "country": country or "",
            "phoneNumber": phoneNumber or "",
            "street": street or "",
            "postal_code": postal_code or "",
            "city": city or "",
            "state": state or "",
            "numberOfSequences": numberOfSequences or "",
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["SequenceSubmission", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_influenza_ha_subtype_conversion_app(api: JsonRpcCaller, token: str = None, user_id: str = None, input_source: str = None, input_fasta_data: str = None, input_fasta_file: str = None, input_feature_group: str = None, input_feature_list: str = None, types: str = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "InfluenzaHASubtypeConversion"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "input_source": input_source or "",
            "input_fasta_data": input_fasta_data or "",
            "input_fasta_file": input_fasta_file or "",
            "input_feature_group": input_feature_group or "",
            "input_feature_list": input_feature_list or "",
            "types": types or "",
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["InfluenzaHASubtypeConversion", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_subspecies_classification_app(api: JsonRpcCaller, token: str = None, user_id: str = None, input_source: str = None, input_fasta_data: str = None, input_fasta_file: str = None, input_genome_group: str = None, ref_msa_fasta: str = None, virus_type: str = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "SubspeciesClassification"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "input_source": input_source or "",
            "input_fasta_data": input_fasta_data or "",
            "input_fasta_file": input_fasta_file or "",
            "input_genome_group": input_genome_group or "",
            "ref_msa_fasta": ref_msa_fasta or "",
            "virus_type": virus_type or "",
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["SubspeciesClassification", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_viral_assembly_app(api: JsonRpcCaller, token: str = None, user_id: str = None, paired_end_lib: Dict = None, single_end_lib: Dict = None, srr_id: str = None, recipe: str = "auto", module: str = None, viral_size: str = "5M", output_path: str = None, output_file: str = None, debug: int = 0) -> str:
    app_name = "ViralAssembly"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "paired_end_lib": paired_end_lib or {},
            "single_end_lib": single_end_lib or {},
            "srr_id": srr_id or "",
            "recipe": recipe,
            "module": module or "",
            "viral_size": viral_size,
            "output_path": output_path,
            "output_file": output_file,
            "debug": debug
        }
        data = ["ViralAssembly", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_fastqutils_app(api: JsonRpcCaller, token: str = None, user_id: str = None, input_status: str = "unaligned", input_type: str = "input_group", fasta_files: List[Dict] = None, select_genomegroup: List[str] = None, feature_groups: List[str] = None, feature_list: List[str] = None, genome_list: List[str] = None, aligner: str = "Muscle", alphabet: str = "dna", fasta_keyboard_input: str = "", ref_type: str = "none", strategy: str = "auto", ref_string: str = "", output_path: str = None, output_file: str = None) -> str:
    app_name = "MSA"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "input_status": input_status,
            "input_type": input_type,
            "fasta_files": fasta_files or [],
            "select_genomegroup": select_genomegroup or [],
            "feature_groups": feature_groups or [],
            "feature_list": feature_list or [],
            "genome_list": genome_list or [],
            "aligner": aligner,
            "alphabet": alphabet,
            "fasta_keyboard_input": fasta_keyboard_input,
            "ref_type": ref_type,
            "strategy": strategy,
            "ref_string": ref_string,
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["MSA", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_genome_alignment_app(api: JsonRpcCaller, token: str = None, user_id: str = None, genome_ids: List[str] = None, recipe: str = "progressiveMauve", seedWeight: float = None, maxGappedAlignerLength: float = None, maxBreakpointDistanceScale: float = None, conservationDistanceScale: float = None, weight: float = None, minScaledPenalty: float = None, hmmPGoHomologous: float = None, hmmPGoUnrelated: float = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "GenomeAlignment"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "genome_ids": genome_ids or [],
            "recipe": recipe,
            "seedWeight": seedWeight,
            "maxGappedAlignerLength": maxGappedAlignerLength,
            "maxBreakpointDistanceScale": maxBreakpointDistanceScale,
            "conservationDistanceScale": conservationDistanceScale,
            "weight": weight,
            "minScaledPenalty": minScaledPenalty,
            "hmmPGoHomologous": hmmPGoHomologous,
            "hmmPGoUnrelated": hmmPGoUnrelated,
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["GenomeAlignment", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_sars_genome_analysis_app(api: JsonRpcCaller, token: str = None, user_id: str = None, paired_end_libs: List[Dict] = None, single_end_libs: List[Dict] = None, srr_ids: List[str] = None, recipe: str = "auto", primers: str = "ARTIC", primer_version: str = None, min_depth: int = 100, max_depth: int = 8000, keep_intermediates: int = 0, output_path: str = None, output_file: str = None, debug_level: int = 0) -> str:
    app_name = "SARS2Assembly"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "paired_end_libs": paired_end_libs or [],
            "single_end_libs": single_end_libs or [],
            "srr_ids": srr_ids or [],
            "recipe": recipe,
            "primers": primers,
            "primer_version": primer_version or "",
            "min_depth": min_depth,
            "max_depth": max_depth,
            "keep_intermediates": keep_intermediates,
            "output_path": output_path,
            "output_file": output_file,
            "debug_level": debug_level
        }
        data = ["SARS2Assembly", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_msa_snp_analysis_app(api: JsonRpcCaller, token: str = None, user_id: str = None, input_status: str = "unaligned", input_type: str = "input_group", fasta_files: List[Dict] = None, select_genomegroup: List[str] = None, feature_groups: List[str] = None, feature_list: List[str] = None, genome_list: List[str] = None, aligner: str = "Muscle", alphabet: str = "dna", fasta_keyboard_input: str = "", ref_type: str = "none", strategy: str = "auto", ref_string: str = "", output_path: str = None, output_file: str = None) -> str:
    app_name = "MSA"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "input_status": input_status,
            "input_type": input_type,
            "fasta_files": fasta_files or [],
            "select_genomegroup": select_genomegroup or [],
            "feature_groups": feature_groups or [],
            "feature_list": feature_list or [],
            "genome_list": genome_list or [],
            "aligner": aligner,
            "alphabet": alphabet,
            "fasta_keyboard_input": fasta_keyboard_input,
            "ref_type": ref_type,
            "strategy": strategy,
            "ref_string": ref_string,
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["MSA", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_metacats_app(api: JsonRpcCaller, token: str = None, user_id: str = None, output_path: str = None, output_file: str = None, p_value: float = 0.05, year_ranges: str = None, metadata_group: str = None, input_type: str = None, alphabet: str = "na", groups: List[str] = None, alignment_file: str = None, group_file: str = None, alignment_type: str = None, auto_groups: List[Dict] = None) -> str:
    app_name = "MetaCATS"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "output_path": output_path,
            "output_file": output_file,
            "p_value": p_value,
            "year_ranges": year_ranges or "",
            "metadata_group": metadata_group or "",
            "input_type": input_type or "",
            "alphabet": alphabet,
            "groups": groups or [],
            "alignment_file": alignment_file or "",
            "group_file": group_file or "",
            "alignment_type": alignment_type or "",
            "auto_groups": auto_groups or []
        }
        data = ["MetaCATS", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_proteome_comparison_app(api: JsonRpcCaller, token: str = None, user_id: str = None, genome_ids: List[str] = None, user_genomes: List[str] = None, user_feature_groups: List[str] = None, reference_genome_index: int = 1, min_seq_cov: float = 0.30, max_e_val: float = 1e-5, min_ident: float = 0.1, min_positives: float = 0.2, output_path: str = None, output_file: str = None) -> str:
    app_name = "GenomeComparison"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "genome_ids": genome_ids or [],
            "user_genomes": user_genomes or [],
            "user_feature_groups": user_feature_groups or [],
            "reference_genome_index": reference_genome_index,
            "min_seq_cov": min_seq_cov,
            "max_e_val": max_e_val,
            "min_ident": min_ident,
            "min_positives": min_positives,
            "output_path": output_path,
            "output_file": output_file
        }
        data = ["GenomeComparison", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_comparative_systems_app(api: JsonRpcCaller, token: str = None, user_id: str = None, output_path: str = None, output_file: str = None, genome_ids: List[str] = None, genome_groups: List[str] = None) -> str:
    app_name = "ComparativeSystems"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "output_path": output_path,
            "output_file": output_file,
            "genome_ids": genome_ids or [],
            "genome_groups": genome_groups or []
        }
        data = ["ComparativeSystems", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_docking_app(api: JsonRpcCaller, token: str = None, user_id: str = None, protein_input_type: str = None, input_pdb: List[str] = None, user_pdb_file: List[str] = None, ligand_library_type: str = None, ligand_named_library: str = None, ligand_smiles_list: List[str] = None, ligand_ws_file: str = None, top_n: int = None, batch_size: int = 10, output_path: str = None, output_file: str = None, enable_debug: bool = None, _tmpdir: str = None) -> str:
    app_name = "Docking"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        params = {
            "protein_input_type": protein_input_type or "",
            "input_pdb": input_pdb or [],
            "user_pdb_file": user_pdb_file or [],
            "ligand_library_type": ligand_library_type or "",
            "ligand_named_library": ligand_named_library or "",
            "ligand_smiles_list": ligand_smiles_list or [],
            "ligand_ws_file": ligand_ws_file or "",
            "top_n": top_n,
            "batch_size": batch_size,
            "output_path": output_path,
            "output_file": output_file,
            "enable_debug": enable_debug,
            "_tmpdir": _tmpdir or ""
        }
        data = ["Docking", params, {}]
        result = api.call("AppService.start_app2", data, 1, token)
        return result
    except Exception as e:
        print(e)
        return []

def start_similar_genome_finder_app(api: JsonRpcCaller, token: str = None, user_id: str = None, selectedGenomeId: str = None, fasta_file: str = None, max_pvalue: float = None, max_distance: float = None, max_hits: int = None, include_reference: bool = None, include_representative: bool = None, include_bacterial: bool = None, include_viral: bool = None, output_path: str = None, output_file: str = None) -> str:
    app_name = "SimilarGenomeFinder"
    try:
        # Set default values if not provided
        output_path, output_file = _set_default_output_paths(user_id, app_name, output_path, output_file)
        
        # Call the Minhash.compute_genome_distance_for_genome2 method
        function_call = ""
        if selectedGenomeId:
            params = [selectedGenomeId, max_pvalue, max_distance, max_hits, include_reference, include_representative, include_bacterial, include_viral]
            function_call = "Minhash.compute_genome_distance_for_genome2"
        elif fasta_file:
            params = [fasta_file, max_pvalue, max_distance, max_hits, include_reference, include_representative, include_bacterial, include_viral]
            function_call = "Minhash.compute_genome_distance_for_fasta2"
        else:
            return "Error: selectedGenomeId or fasta_file is required"
        
        result = api.call(function_call, params, 1, token)
        output_headers = ['genome_id','distance','pvalue','kmers']
        return output_headers, result
    except Exception as e:
        print(e)
        return []