
# SciToolEval Dataset and Evaluation Guide

## 📂 Dataset Structure

The `data/` directory contains the **SciToolEval question data**, categorized into two difficulty levels:

- `level1_question.jsonl`
- `level2_question.jsonl`

### 📄 Example Entry

```json
{
  "tool_path": ["SMILESToInChI", "InChIToInChIKey"],
  "Parameter": "CC(C)OC(=O)C1=CN=C(N=C1C2=CN(C3=CC=CC=C32)C)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC",
  "question": "What is the InChIKey for the molecule represented by the SMILES notation CC(C)OC(=O)C1=CN=C(N=C1C2=CN(C3=CC=CC=C32)C)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC?",
  "answer": "The InChIKey is AZSRSNUQCUDCGG-UHFFFAOYSA-N."
}
````

---

## 🧪 Evaluation Instructions

The `eval/` directory includes evaluation scripts for:

* **Tool planning accuracy**
* **Final answer accuracy**

To use the evaluation scripts, your model predictions should be formatted to match:

* `data/example_standard_toolpath.jsonl`
* `data/example_standard_answer.jsonl`

### 📄 Input Format Examples

**Tool Path Prediction:**

```json
{
  "tool_path": ["NameToSMILES", "GetChi1v", "GetMolFormula", "CalculateTPSA"],
  "question": "Please provide the SMILES string for the molecule with the molecule name \"morphine,\" then calculate and report the Chi^1v valence molecular graph index, the molecular formula, and the topological polar surface area (TPSA) for the molecule."
}
```

**Final Answer Prediction:**

```json
{
  "question": "What is the InChIKey for the molecule represented by the SMILES notation CC(C)OC(=O)C1=CN=C(N=C1C2=CN(C3=CC=CC=C32)C)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC?",
  "final_answer": "The InChIKey for the molecule represented by the SMILES notation CC(C)OC(=O)C1=CN=C(N=C1C2=CN(C3=CC=CC=C32)C)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC is AZSRSNUQCUDCGG-UHFFFAOYSA-N."
}
```
