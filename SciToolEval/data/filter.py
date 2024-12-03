import json



# 先读取所有问题数据到内存中
with open('level2_question.jsonl', 'r', encoding='utf-8') as f:
    questions_data = {json.loads(line)['question'] for line in f.readlines()}

# 处理答案数据并检查问题是否存在
with open('level2_answer.jsonl', 'r', encoding='utf-8') as f:
    for line in f:
        data = json.loads(line)
        question = data['question']
        # 检查问题是否在 questions_data 集合中
        if question not in questions_data:
            print(question)
