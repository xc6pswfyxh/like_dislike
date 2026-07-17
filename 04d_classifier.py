# importing libraries
import pandas as pd
import ollama
import json
import re
from tqdm import tqdm


# loading the complete comments dataset
data_cpath = "data/comments/comments.csv"
data_comments = pd.read_csv(data_cpath, encoding="latin1")

# loading manually coded dataset to use manually coded observations for the final dataset 
data_epath = "data/comments/comments_manual_coding.xlsx"
data_eval = pd.read_excel(data_epath, engine="openpyxl")

VARS = ["pers_exp", "emot_exp", "pol_opin", "breadth", "valence", "contr"]

if "coder" not in data_comments.columns:
    data_comments["coder"] = None


# read the system prompt from file
with open("data/comments/system_prompt.txt", "r", encoding="utf-8") as f:
    SYSTEM_PROMPT = f.read()

def classify_comment(text, parent_text=None, model="qwen3:8b"):
    fallback = {"pers_exp": -1, "emot_exp": -1, "pol_opin": -1,
                "breadth": -1, "valence": -1, "contr": -1}
    if not isinstance(text, str) or text.strip() == "":
        return {"pers_exp": 0, "emot_exp": 0, "pol_opin": 0,
                "breadth": 0, "valence": 4, "contr": 0}
    user_msg = f"Comment: {text}"
    if parent_text:
        user_msg += f"\n\n(This is a reply to the following parent comment: {parent_text})"
    response = ollama.chat(
        model=model,
        messages=[{"role": "system", "content": SYSTEM_PROMPT},
                   {"role": "user", "content": user_msg}],
        options={"temperature": 0}
    )
    out = response["message"]["content"].strip()
    match = re.search(r"\{.*\}", out, re.DOTALL)
    if not match:
        return fallback
    try:
        parsed = json.loads(match.group(0))
    except json.JSONDecodeError:
        return fallback
    result = {}
    bounds = {"pers_exp": (0,1), "emot_exp": (0,1), "pol_opin": (0,1),
              "breadth": (0,3), "valence": (1,4), "contr": (0,1)}
    for key, (lo, hi) in bounds.items():
        val = parsed.get(key, -1)
        try:
            val = int(val)
        except (TypeError, ValueError):
            val = -1
        result[key] = val if lo <= val <= hi else -1
    return result


# composite id + parent lookup 
def make_id(topic, post_number):
    return f"{topic}_{int(post_number)}"

data_comments["row_id"] = data_comments.apply(lambda r: make_id(r["topic"], r["post_number"]), axis=1)
data_eval["row_id"] = data_eval.apply(lambda r: make_id(r["topic"], r["post_number"]), axis=1)

reply_col = "reply_to_post_number" if "reply_to_post_number" in data_comments.columns else None
text_by_id = dict(zip(data_comments["row_id"], data_comments["raw"].fillna("")))

def get_parent_text(row):
    if reply_col and pd.notna(row.get(reply_col)):
        parent_id = make_id(row["topic"], row[reply_col])
        return text_by_id.get(parent_id)
    return None

# step 1: fill in the manual codes wherever they exist 
manual_lookup = data_eval.set_index("row_id")[VARS].to_dict(orient="index")

for i, row in data_comments.iterrows():
    if row["row_id"] in manual_lookup:
        for var in VARS:
            data_comments.at[i, var] = manual_lookup[row["row_id"]][var]
        data_comments.at[i, "coder"] = "manual"

# step 2: run llm only on blank obs
FINAL_MODEL = "qwen3:8b"

texts = data_comments["raw"].fillna("")
to_classify = data_comments[data_comments["coder"] != "manual"]
print(f"{len(to_classify)} comments need LLM classification "
      f"({len(data_comments) - len(to_classify)} already manually coded)")

for i, row in tqdm(to_classify.iterrows(), total=len(to_classify), desc=f"Classifying ({FINAL_MODEL})"):
    parent_text = get_parent_text(row)
    codes = classify_comment(texts.loc[i], parent_text=parent_text, model=FINAL_MODEL)
    for var in VARS:
        data_comments.at[i, var] = codes[var]
    data_comments.at[i, "coder"] = "llm"


# saving
data_comments.to_excel("data/comments/data_comments_classified.xlsx", index=False)
print("Saved data_comments_classified.xlsx")