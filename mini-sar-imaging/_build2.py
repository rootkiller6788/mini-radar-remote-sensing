import os
BASE = os.path.dirname(os.path.abspath(__file__))
def w(rel, content):
    fp = os.path.join(BASE, rel)
    os.makedirs(os.path.dirname(fp), exist_ok=True)
    with open(fp, "w", encoding="utf-8") as f: f.write(content)
    print(f"  {rel}: {len(content.splitlines())} lines")
