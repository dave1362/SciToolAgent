```plaintext    
/SCITOOLAGENT
├── .github                #  GitHub profile
├── app                    # App Directory
│   ├── core                   
│   │   ├── __pycache__
│   │   ├── __init__.py
│   │   └── config.py          # Configuration file
│   └── llms                   # configuration file
│       ├── api_clients        # API client module
│       │   ├── __pycache__
│       │   ├── __init__.py
│       │   └── openai_api.py  # OpenAI API
│       ├── extraction         # Extraction module
│       │   ├── __pycache__
│       │   ├── __init__.py
│       │   ├── extract_answer.py # Answer extraction tool
│       │   └── extract_tool.py   # tool Extraction tool
│       ├── planning           #  Planning module
│       │   ├── __pycache__
│       │   ├── __init__.py
│       │   ├── prompt_templates.py # Prompt template
│       │   └── tool_plan_executor.py # Tool plan executor
│       └── tools              # Tool module
│           ├── __pycache__
│           ├── safety         # Safety Check module
│           │   ├── __pycache__
│           │   ├── __init__.py
│           │   └── toxicity_checker.py # Toxicity checker tools
│           ├── __init__.py
│           ├── utils.py       #  Common method of the utility class
│           └── tool_runner.py # Tool run module
├── data                   # Data directory
│   ├── tool_example.xlsx   # Tool KG example data
│   ├── toxin_compound.csv # Toxic compounds data
│   └── toxin_protein.csv  # Toxic protein data
├── KG                     # Knowledge Graph Storage Directory
│   ├── storage_graph_large
├── scripts                # Script directory
│   ├── __pycache__
│   ├── __init__.py
│   ├── custom_kg_retrievers.py #  Custom Knowledge Graph Retrievers
│   ├── generate_kg_index.py   # Knowledge Graph index generation script
│   ├── retrieve_tool_info.py  # Tool information retrieval script
│   └── utils.py              # script tool method
├── test                   # Test Directory
│   ├── __pycache__
│   ├── __init__.py
│   ├── test_embedding.py    # embedding model test
│   ├── test_openai.py       # OpenAI function test
│   └── test_run_SciToolAgent.py # SciToolAgent Runs the test
├── ToolsAgent             # ToolsAgent service Directory
├── app.py                 # Main application entry
├── Cases.ipynb        # Case study notebook
├── example.env            # example file of environment variables
├── .gitignore             # Git Ignore files
├── requirements_agent.txt # List of project dependencies for agent only
└── requirements.txt       # List of project dependencies
```