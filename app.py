from fastapi import FastAPI
from pydantic import BaseModel
from scripts.retrieve_tool_info import run_SciToolAgent
from app.core.config import Config
import uvicorn

def create_app():
    app = FastAPI()
    
    class ChatQuery(BaseModel):
        question: str
        api_key: str
        api_base: str
        is_add_filename: bool = False
        file_path_list: list

    @app.post("/chat")
    async def chat(query: ChatQuery):
        try:
            question = query.question
            is_add_filename = query.is_add_filename
            file_path_list = query.file_path_list
            print(f"question: {question}\nis_add_filename: {is_add_filename}\nfile_path_list: {file_path_list}\n")
            Config().OPENAI_API_KEY = query.api_key
            Config().OPENAI_API_BASE = query.api_base
            # print(f"OPENAI_API_KEY: {Config().OPENAI_API_KEY}, OPENAI_API_BASE: {Config().OPENAI_API_BASE}")
            response = await run_SciToolAgent(question, file_path_list, is_add_filename=is_add_filename)
            return {"response": response}
        except Exception as e:
            return {
                "error_type": type(e).__name__,
                "message": str(e)
            }
    
    return app

app = create_app()

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=50005)
