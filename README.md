# BV-BRC Data MCP Server

A Model Context Protocol (MCP) server that provides access to BV-BRC (Bacterial and Viral Bioinformatics Resource Center) data through the bvbrc-solr-python-api module.

## Features

- Access to genome and genome feature data
- Query capabilities for various BV-BRC data types
- Modular tool architecture for different data categories
- Flask-based MCP server implementation

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
cd bvbrc-service
pip install -r requirements.txt
```

## Configuration

The server uses a `config.json` file for configuration:

```json
{
    "service_api_url": "https://p3.theseed.org/services/app_service",
    "port": 8058,
    "token": "<bvbrc_token>"
}
```

## Usage

Run the MCP server:

```bash
python server.py
```

The server will start on port 8059 (configurable in `config.json`).

## Health Check

The server provides a health check endpoint at `/health` that returns the server status.

## Dependencies

- flaskmcp
- httpx >= 0.28
- bvbrc-solr-python-api

## License

This project is part of the BV-BRC MCP Servers collection.
