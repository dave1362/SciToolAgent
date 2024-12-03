from llama_index.embeddings.openai import OpenAIEmbedding
from app.core.config import Config


text = "This is an example sentence for embedding."

embedding_instance = OpenAIEmbedding(
    model=Config().EMBEDDING_MODEL_NAME,  
    api_key=Config().OPENAI_API_KEY,
    api_base=Config().OPENAI_API_BASE
)

try:
    embedding = embedding_instance._get_text_embedding(text)
    print(embedding)
except Exception as e:
    print(f"Error occurred: {e}")

