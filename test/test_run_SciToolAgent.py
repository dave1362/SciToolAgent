import asyncio
from scripts.retrieve_tool_info import run_SciToolAgent

async def main():
    question = "What is the SMILES of ethanol?"
    file_path_list = []
    is_add_filename = False

    response = await run_SciToolAgent(question, file_path_list, is_add_filename)
    print("Response from run_SciToolAgent:", response)

if __name__ == "__main__":
    asyncio.run(main())
