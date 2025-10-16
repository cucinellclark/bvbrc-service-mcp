from fastmcp import FastMCP
from json_rpc import JsonRpcCaller
from service_tools import register_all_tools
from token_provider import TokenProvider
import json
import sys

with open('config.json', 'r') as f:
    config = json.load(f)

service_api_url = config['service_api_url']
similar_genome_finder_api_url = config.get('similar_genome_finder_api_url', service_api_url)
port = config.get("port", 5001)
mcp_url = config.get("mcp_url", "127.0.0.1")

# Initialize token provider for HTTP mode
token_provider = TokenProvider(mode="http")

api = JsonRpcCaller(service_api_url)
similar_genome_finder_api = JsonRpcCaller(similar_genome_finder_api_url)

mcp = FastMCP("BVBRC Service MCP Server")

# Register all tools with the MCP server and token provider
register_all_tools(mcp, api, similar_genome_finder_api, token_provider)

def main() -> int:
    print(f"Starting BVBRC Service MCP FastMCP Server on port {port}...", file=sys.stderr)
    try:
        mcp.run(transport="http", host=mcp_url, port=port)
    except KeyboardInterrupt:
        print("Server stopped.", file=sys.stderr)
    except Exception as e:
        print(f"Server error: {e}", file=sys.stderr)
        return 1
    
    return 0


if __name__ == "__main__":
    main()
