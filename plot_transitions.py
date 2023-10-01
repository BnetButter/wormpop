import pandas as pd
import matplotlib.pyplot as plt

# Read the TSV file
df = pd.read_csv('stage_transitions.tsv', sep='\t')

# Plotting the data
plt.figure(figsize=(10,6))

for column in df.columns:
    if column != 'Timestep':
        plt.plot(df['Timestep'], df[column], label=column)

plt.xlabel('Timestep')
plt.ylabel('Value')
plt.legend()
plt.title('Timestep vs Other Fields')
plt.grid(True)
plt.tight_layout()
plt.savefig("stage_transitions.png")
