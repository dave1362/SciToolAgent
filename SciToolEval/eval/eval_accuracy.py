import pandas as pd
import json
import os
import openai
import argparse
from tqdm import tqdm

# Set environment variables for OpenAI API
os.environ["OPENAI_API_BASE"] = "Your api base"
os.environ["OPENAI_API_KEY"] = "Your api key"


def get_chat_messages(prompt_text, api_key):
    """
    Call the OpenAI ChatCompletion API to get the model's response.
    """
    openai.api_key = api_key
    try:
        response = openai.ChatCompletion.create(
            model="gpt-4o",
            messages=prompt_text,
            max_tokens=2000,
            temperature=0.1,
            response_format={"type": "json_object"}
        )
        content = response.choices[0].message.content.strip()
        if not content:
            raise ValueError("Empty response from OpenAI API")
        return content
    except Exception as e:
        print(f"Error during API call: {e}")
        return None


def load_jsonl(file_path):
    """
    Load a JSONL file and return a list of data for each line.
    """
    with open(file_path, 'r', encoding='utf-8') as file:
        return [json.loads(line) for line in file]


def evaluate_responses(input_file, standard_file, output_file, api_key):
    """
    Evaluate the correctness of responses and calculate accuracy.
    """
    input_data = load_jsonl(input_file)
    standard_data = {item['question']: item['final_answer'] for item in load_jsonl(standard_file)}
    results = []
    correct_count = 0

    for item in tqdm(input_data, total=len(input_data)):
        question = item.get('question', '')
        response = item.get('final_answer', '')

        # Get the standard answer
        correct_answer = standard_data.get(question, '')

        # Build the prompt
        prompt = f"""You will be given a scientific question, a standard answer, and a response of unknown correctness. Please evaluate the response based on the question and the standard answer, using the following criteria:

If the value is a prediction or a calculation result, it allows for a margin of error (e.g., TPSA), then deviations within the allowable error range (-10% to +10%) can be considered "correct". Significant deviations are considered "incorrect".
As long as the response conveys the same core concepts, regardless of specific wording, it can be considered "correct". If the core content is inconsistent or incorrect, label it as "incorrect."

Question: {question}
Standard answer: {correct_answer}
Response: {response}

The response only needs to contain the correct answer; it does not have to be identical. There may be multiple correct ways to express the same concept like SMILES and name of a molecule, but it must contain the same core information to answer the question correctly.
Please make your judgment based on the above criteria and respond only with JSON format like:
{{
"question": [The question you are evaluating.],
"reason": [The reason why you think the response is correct or incorrect.],
"correctness": [correct/incorrect]
}}"""

        prompt_text = [
            {"role": "system", "content": "You are an expert in science."},
            {"role": "user", "content": prompt}
        ]

        try:
            # Get the model's response
            response_content = get_chat_messages(prompt_text, api_key)

            if not response_content:
                raise ValueError("No response content received.")

            # Try to parse the response as JSON
            try:
                response_json = json.loads(response_content)
            except json.JSONDecodeError as json_error:
                print(f"Invalid JSON response: {response_content}")
                response_json = {"correctness": "incorrect", "reason": "Invalid JSON response"}

            # Record the result
            correctness = response_json.get("correctness", "incorrect")
            results.append({
                "question": question,
                "correct_answer": correct_answer,
                "response": response,
                "reason": response_json.get("reason", "No reason provided"),
                "correctness": correctness
            })

            if correctness == "correct":
                correct_count += 1

        except Exception as e:
            print(f"Error processing question: {question}, Error: {e}")

    # Save the results
    with open(output_file, 'w', encoding='utf-8') as file:
        for result in results:
            file.write(json.dumps(result, ensure_ascii=False) + '\n')

    # Calculate accuracy
    accuracy = correct_count / len(input_data) if input_data else 0
    print(f"Accuracy: {accuracy:.2%}")
    return accuracy


if __name__ == "__main__":
    # Setup argument parser


    parser = argparse.ArgumentParser(description="Evaluate tool paths for scientific questions.")
    parser.add_argument('--input_file', type=str, required=True, help="Path to the result JSONL file")
    parser.add_argument('--standard_file', type=str, required=True, help="Path to the standard answers JSONL file")
    parser.add_argument('--output_file', type=str, required=True, help="Path to save the evaluation results")

    args = parser.parse_args()



    # Run evaluation
    evaluate_responses(args.input_file, args.standard_file, args.output_file)
