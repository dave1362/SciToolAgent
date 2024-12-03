import pandas as pd
import json
import os
import openai
import yaml
from tqdm import tqdm

# Set environment variables for OpenAI API
os.environ["OPENAI_API_BASE"] = "Your api base"
os.environ["OPENAI_API_KEY"] = "Your api key"


# Function to get the chat message response from GPT-4 model
def get_chat_messages(prompt_text):
    openai.api_key = os.getenv("OPENAI_API_KEY")
    openai.api_base = os.getenv("OPENAI_API_BASE")
    response = openai.ChatCompletion.create(
        model="gpt-4o",  # Use GPT-4 model (or any model you prefer)
        messages=prompt_text,
        max_tokens=2000,
        n=1,
        stop=None,
        temperature=0.1,
        response_format={"type": "json_object"},

    )
    return response.choices[0].message.content


# Function to get the description of a tool based on its name from tools_data
def get_tool_description(tool_name, tools_data):
    # Check if tool exists in the tools_data dictionary
    if tool_name not in tools_data:
        return f"No description found for tool: {tool_name}"

    tool_info = tools_data[tool_name]
    description = []

    # Construct description excluding certain relationships
    for entry in tool_info:
        if entry[1] not in ["is the input of", "is the output of"]:
            description.append(f"{entry[0]} {entry[1]} {entry[2]}")

    return "\n".join(description)


# Function to calculate the correctness percentage based on the comparison between tool paths
def calculate_correctness_percentage(correct_count, total_count):
    return (correct_count / total_count) * 100


# Main function to evaluate tool paths and calculate correctness
def evaluate_tool_paths(result_file, standard_file,tool_description_file, output_file):
    correct_count = 0
    total_count = 0

    # Load results data
    with open(result_file, 'r', encoding='utf-8') as file:
        lines = file.readlines()

        for line in tqdm(lines, total=len(lines)):
            data = json.loads(line)
            question = data['question']
            tool_path = data['tool_path'] if data['tool_path'] else []

            # Load standard data from jsonl file for comparison
            with open(standard_file, 'r', encoding='utf-8') as f:
                standard_data = f.readlines()
                for item in standard_data:
                    standard_data = json.loads(item)
                    if standard_data['question'] == question:
                        standard_tool = standard_data['Action']
                        break

            if standard_tool is None:
                continue  # Skip if no standard tool path is found for the question

            # Compare tool paths
            if tool_path == standard_tool:
                correct_count += 1
                with open(output_file, 'a', encoding='utf-8') as file:
                    result = {"question": data['question'], "result": "Correct"}
                    file.write(json.dumps(result) + '\n')
            else:
                # Evaluate using GPT-4 if the tool path is incorrect
                prompt = f"""You will be given a scientific question, some scientific functions, a tool path offered by the large language model and a tool path annotated by a human expert. You should evaluate whether the tool path is able to solve the question based on the functions.

                Question: {question}
                Tool path by human: {standard_tool}
                Tool path by LLM: {tool_path}

                Functions description:
                """
                with open(tool_description_file, 'r', encoding='utf-8') as f:
                    tools_data = json.load(f)

                # Add descriptions of the tools in the tool path
                for name in tool_path:
                    prompt += get_tool_description(name, tools_data) + "\n"

                prompt += """Evaluate the tool path using these criteria:

Logicality: The tool path must follow a logical sequence and align with the problem, ensuring the chosen tools lead to a correct solution.
Completeness: The tool path should fully address the problem. Missing key tools that impact the answer's completeness will lower the score.
Sequence: Tools should be used in the correct order when necessary, where one tool’s output serves as the input for the next. In cases where order is not critical, the sequence can be flexible.
Conciseness: The solution should use the minimal number of tools without sacrificing completeness. Redundant or unnecessary tools should be avoided.
Please give your judgment based on the above criteria and respond only with json format like:
{
    "question": "[The question you are evaluating.]",
    "result": "[Correct/Incorrect]",
}
"""

                prompt_text = [{"role": "system", "content": "You are an expert in science."},
                               {"role": "user", "content": f"{prompt}"}]

                response_content = get_chat_messages(prompt_text)

                # Determine the result based on the response
                result_data = json.loads(response_content)
                if result_data.get("result") == "Correct":
                    correct_count += 1

                with open(output_file, 'a', encoding='utf-8') as file:
                    result = {"question": data['question'], "result": "Correct" if result_data.get("result") == "Correct" else "Incorrect"}
                    file.write(json.dumps(result) + '\n')


            total_count += 1  # Increment the total count for each evaluated question

    # Calculate and print correctness percentage
    correctness_percentage = calculate_correctness_percentage(correct_count, total_count)
    print(f"Correctness Percentage: {correctness_percentage:.2f}%")




# Run the evaluation
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Evaluate tool paths for scientific questions.")
    parser.add_argument('--result_file', type=str, required=True, help="Path to the result JSONL file")
    parser.add_argument('--standard_file', type=str, required=True, help="Path to the standard answers JSONL file")
    parser.add_argument('--tool_description_file', type=str, required=True,
                        help="Path to the tool description JSON file")
    parser.add_argument('--output_file', type=str, required=True, help="Path to save the evaluation results")

    args = parser.parse_args()

    evaluate_tool_paths(result_file=args.result_file, standard_file=args.standard_file, tool_description_file=args.tool_description_file, output_file=args.output_file)
