import requests
import json

with open("config.json", "r") as f:
    config = json.load(f)

token = config["token"]

params = {
    "token": token
}

call_url = "http://127.0.0.1:8058/mcp/tools/call"
list_url = "http://127.0.0.1:8058/mcp/tools/list"


if True:
    print('********* listing tools *********')

    url = list_url

    response = requests.get(url)

    print(response.json())

    print('********* calling tool: service_enumerate_apps *********')

    response = requests.get(url)
    print(response.json())

if False:
    print('********* calling tool: service_start_date_app *********')

    params = {
        "token": token,
        "params": {
            "output_path": "/clark.cucinell@patricbrc.org/home/MCP_Dev",
            "output_file": "test_mcp_date"
        }
    }

    #response = requests.post(url, json={"jsonrpc": "2.0", "id": 1, "name": "service_start_date_app", "params": params})
    #print(response.json())

if False:
    print('********* calling tool: service_start_genome_annotation_app *********')

    params = {
        "token": token,
        "params": {
            "output_path": "/clark.cucinell@patricbrc.org/home/MCP_Dev",
            "output_file": "Pseudomonas aeruginosa PAO1 test_mcp",
            "contigs": "/clark.cucinell@patricbrc.org/home/DevTest/Pao1_Test/.Test_HTSeq/208964.12.fna",
            "scientific_name": "Pseudomonas aeruginosa PAO1 test_mcp",
            "tax_id": "208964",
            "my_label": "test_mcp",
            "reference_genome_id": "",
            "taxonomy_id": "208964"
        }
    }

    response = requests.post(call_url, json={"jsonrpc": "2.0", "id": 1, "name": "service_start_genome_annotation_app", "params": params})
    print(response.json())

if False:
    print('********* calling tool: service_query_tasks *********')

    params = {
        "token": token,
        "params": {
            "task_ids": ["19018027"]
        }
    }

    response = requests.post(call_url, json={"jsonrpc": "2.0", "id": 1, "name": "service_query_tasks", "params": params})
    print(response.json())