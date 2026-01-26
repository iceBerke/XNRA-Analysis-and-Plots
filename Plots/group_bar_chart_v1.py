import matplotlib.pyplot as plt
import numpy as np

# Define your data
ranges = ['[0-50]', '[50-100]', '[100-150]', '[150-200]']
treatments = ['20°C', '30°C', '40°C']

data = {
    '20°C': [15, 35, 30, 20],
    '30°C': [25, 40, 25, 10],
    '40°C': [30, 25, 35, 10]
}

# Key parameters to control spacing
n_treatments = len(treatments)
x = np.arange(len(ranges))  # the label locations
width = 0.1  # width of EACH individual bar (make this smaller if bars overlap)
group_spacing = 1.0  # total width allocated to each group (default is 0.8)

# Calculate positions so bars are centered within each group
fig, ax = plt.subplots(figsize=(20, 12))

for i, (treatment, percentages) in enumerate(data.items()):
    # Center the bars around each x position
    offset = (i - n_treatments/2 + 0.5) * width
    ax.bar(x + offset, percentages, width, label=treatment)

# Customize the plot
ax.set_xlabel('Value Ranges', fontsize=12)
ax.set_ylabel('Percentage (%)', fontsize=12)
ax.set_title('Data Distribution by Temperature Treatment', fontsize=14)
ax.set_xticks(x)
ax.set_xticklabels(ranges)
ax.legend(title='Temperature', loc='upper right')
ax.set_ylim(0, 100)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.show()