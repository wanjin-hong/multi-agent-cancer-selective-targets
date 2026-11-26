import pandas as pd, numpy as np, sys
f = "/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv"
df = pd.read_csv(f, index_col=0)
# columns are "GENE (Entrez)"; strip suffix
df.columns = [c.split(" (")[0] for c in df.columns]
print("lines:", df.shape[0], "genes:", df.shape[1])
for anchor in ["WWTR1","KLF5"]:
    a = df[anchor]
    # Pearson vs every gene, pairwise NaN-aware
    corrs = df.corrwith(a, method="pearson")
    n = df.notna().sum()
    out = pd.DataFrame({"gene": corrs.index, "r": corrs.values, "n_lines": n.values})
    out = out.dropna().sort_values("r", ascending=False)
    out.to_csv(f"codepend_{anchor}.csv", index=False)
    print(anchor, "top 25 positive:")
    print(out.head(25).to_string(index=False))
    print(anchor, "top 10 negative:")
    print(out.tail(10).iloc[::-1].to_string(index=False))
    # cross-correlation between the two anchors
    print("WWTR1-KLF5 r:", df["WWTR1"].corr(df["KLF5"]))
