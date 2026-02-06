import pandas as pd

t1 = pd.read_csv("inputs/LINCS.csv")
t2 = pd.read_csv("inputs/JUMP-CP.csv")

s1 = set(t1["Pubchem CID"].dropna().astype(float).astype(int))
s2 = set(t2["Pubchem CID"].dropna().astype(float).astype(int))

common_cids = s1 | s2

df_out = pd.DataFrame(sorted(common_cids), columns=["cid"])
df_out.to_csv("outputs/cid_unions.csv", index=False)
