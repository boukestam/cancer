import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

df = pd.read_json("results/signatures.json")

plt.figure(figsize=(10, 6))
sns.barplot(x=df.index, y=df["contribution"])
plt.xlabel("Mutation Signature")
plt.ylabel("Contribution")
plt.title("Mutation Signature Contribution in Breast Cancer Samples")
plt.xticks(rotation=45)
plt.savefig("results/signature_distribution.png")
plt.show()
