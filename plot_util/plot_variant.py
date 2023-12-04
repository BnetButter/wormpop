import pandas as pd
import matplotlib.pyplot as plt
import sys

input_file, output_file = sys.argv[1:]

# Reading the TSV data into a DataFrame
data = pd.read_csv(input_file, sep='\t')

# Plotting the data
plt.figure(figsize=(10, 6))
for column in data.columns:
    plt.plot(data[column], label=column)

# Labeling the axes
plt.xlabel("Timestep")
plt.ylabel("Count")

# Adding a legend
plt.legend()

# Saving the plot as a PNG file'
plt.savefig(output_file)

output_file
