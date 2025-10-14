from fastmcp import FastMCP
from json_rpc import JsonRpcCaller
from service_tools import register_all_tools
import sys
import os

service_api_url = os.getenv("SERVICE_API_URL")
similar_genome_finder_api_url = os.getenv("SIMILAR_GENOME_FINDER_API_URL")

api = JsonRpcCaller(service_api_url)
similar_genome_finder_api = JsonRpcCaller(similar_genome_finder_api_url)

mcp = FastMCP("BVBRC Service MCP Server")

# Register all tools with the MCP server
register_all_tools(mcp, api, similar_genome_finder_api)

def main() -> int:
    print("Starting BVBRC Service MCP FastMCP Server in stdio mode...", file=sys.stderr)
    try:
        mcp.run(transport="stdio")
    except KeyboardInterrupt:
        print("Server stopped.", file=sys.stderr)
    except Exception as e:
        print(f"Server error: {e}", file=sys.stderr)
        return 1
    
    return 0


if __name__ == "__main__":
    main()
