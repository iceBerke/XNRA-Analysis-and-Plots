import matplotlib.pyplot as plt
import numpy as np

# Define your data
ranges = ['[0-50]', '[50-100]', '[100-150]', '[150-200]']
treatments = ['20°C', '30°C', '40°C']

# Data: percentages for each treatment at each range
# Each row = one treatment, each column = one range
data = {
    '20°C': [15, 35, 30, 20],
    '30°C': [25, 40, 25, 10],
    '40°C': [30, 25, 35, 10]
}

# Set up the bar positions
x = np.arange(len(ranges))  # the label locations
width = 0.25  # width of each bar
multiplier = 0

fig, ax = plt.subplots(figsize=(10, 6))

# Create bars for each treatment
for treatment, percentages in data.items():
    offset = width * multiplier
    ax.bar(x + offset, percentages, width, label=treatment)
    multiplier += 1

# Customize the plot
ax.set_xlabel('Value Ranges', fontsize=12)
ax.set_ylabel('Percentage (%)', fontsize=12)
ax.set_title('Data Distribution by Temperature Treatment', fontsize=14)
ax.set_xticks(x + width)
ax.set_xticklabels(ranges)
ax.legend(title='Temperature', loc='upper right')
ax.set_ylim(0, 100)  # for percentages
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.show()
