# BV-BRC Service MCP Server

A Model Context Protocol (MCP) server that provides access to BV-BRC (Bacterial and Viral Bioinformatics Resource Center) services using FastMCP stdio.

## Features

- Submit bvbrc tools and services through convenient functions
- FastMCP stdio-based server for better performance and compatibility
- Comprehensive genomics, metagenomics, and transcriptomics analysis tools

## Installation

### Prerequisites

- Python >= 3.10
- pip

### Install Required Dependencies

0. Create and activate a virtual environment:
```bash
python3 -m venv .venv
source .venv/bin/activate
```

### Install This MCP Server

```bash
git clone https://github.com/cucinellclark/bvbrc-service-mcp.git
cd bvbrc-service-mcp
pip install -r requirements.txt
```

## Configuration

The server uses a `config.json` file for configuration:

```json
{
    "service_api_url": "https://p3.theseed.org/services/app_service"
}
```

## Usage

Run the MCP server:

```bash
python server.py
```

The server uses FastMCP stdio for communication, making it compatible with MCP clients that support stdio transport.

## Available Tools

The server provides access to various BV-BRC analysis tools including:

- **Genomics Analysis**: Genome assembly, annotation, comprehensive analysis
- **Metagenomics**: Taxonomic classification, binning, read mapping
- **Transcriptomics**: RNA-seq analysis, expression import
- **Phylogenomics**: Phylogenetic trees, core genome MLST, SNP analysis
- **Viral Analysis**: SARS-CoV-2 analysis, viral assembly, influenza analysis
- **Utilities**: Sequence submission, primer design, variation analysis
