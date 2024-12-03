from openai import OpenAI
from config import Config

client = OpenAI(
    api_key=Config().OPENAI_API_KEY,
    base_url=Config().OPENAI_API_BASE
)

response = client.chat.completions.create(

    model = Config().MODEL_NAME,
    messages=[
        {"role": "system", "content": "You are a helpful assistant."},
        {"role": "user", "content": "There are 9 birds in a tree. If you shoot one down, how many are left?"},

    ]
)
print(response.choices[0].message.content)

