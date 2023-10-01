import pandas as pd
import matplotlib.pyplot as plt

# Read in the TSV file
data = pd.read_csv('summary_4.tsv', sep='\t')

# Plot the specified columns against the "Timestep" column
plt.figure(figsize=(10,6))

cols_to_plot = ["Number Worms", "Number Eggs", "Number Larvae", "Number Dauer", "Number Adults", "Number Parlads"]

for col in cols_to_plot:
    plt.plot(data["Timestep"], data[col], label=col)

# Setting labels, title and legend
plt.xlabel('Timestep')
plt.ylabel('Count')
plt.title('Counts as a function of Timestep')
plt.legend()

# Display the plot
plt.tight_layout()
plt.grid(True)
plt.savefig("result.png")
