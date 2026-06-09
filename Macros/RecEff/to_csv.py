import pandas as pd

name = "events_f0"

file = f"./Inputs/{name}.dat"

df = pd.read_csv(file, delimiter=" ", names=["run", "entry"])
# Create empty columns for efficiency
df["type"] = ""
df["status"] = False
print(df)
# Save
df.to_csv(f"{name}.csv", index=False)
