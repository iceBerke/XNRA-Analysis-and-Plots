import matplotlib.pyplot as plt
import numpy as np

# Define your data
ranges = ['[0-7]', '[7-35]', '[35-62]', '[62-90]', '[90-118]','[118-145]','[145-173]','[173-201]','[201-229]','[229-256]','[256-700]','[700-3470]']
treatments = ['NT', '250°C', '300°C','400°C','500°C','550°C','600°C']

data = {
  'NT': [8.23, 8.13, 8.43, 8.43, 8.33, 8.38, 8.38, 8.43, 8.38, 8.43, 8.33, 8.33],
'250°C': [4.03, 2.73, 7.53, 8.43, 8.33, 8.38, 8.38, 8.43, 8.38, 8.43, 8.33, 8.33],
'300°C': [3.00, 1.52, 3.80, 4.03, 6.33, 8.33, 8.33, 8.33, 8.33, 8.33, 8.33, 8.33],
'400°C': [2.50, 2.02, 3.03, 5.63, 6.23, 7.23, 7.83, 7.53, 8.33, 8.33, 8.33, 8.33],
'500°C': [2.00, 1.52, 2.23, 4.03, 4.53, 5.53, 6.33, 6.63, 6.83, 7.23, 8.33, 8.33],
'550°C': [2.50, 2.02, 3.03, 5.23, 5.53, 6.43, 6.83, 6.93, 7.23, 7.63, 8.33, 8.33],
'600°C': [1.20, 2.00, 3.03, 4.43, 4.43, 5.83, 6.03, 6.03, 6.53, 6.93, 8.33, 8.33],
}

# Key parameters to control spacing
n_treatments = len(treatments)
x = np.arange(len(ranges))  # the label locations
width = 0.1  # width of EACH individual bar (make this smaller if bars overlap)
group_spacing = 1.0  # total width allocated to each group (default is 0.8)

# Calculate positions so bars are centered within each group
fig, ax = plt.subplots(figsize=(19, 11.5))

for i, (treatment, percentages) in enumerate(data.items()):
    # Center the bars around each x position
    offset = (i - n_treatments/2 + 0.5) * width
    ax.bar(x + offset, percentages, width, label=treatment)

# Customize the plot
ax.set_xlabel('Penetration Depth (nm)', fontsize=10)
ax.set_ylabel('Percentage Sodium (%)', fontsize=10)
#ax.set_title('Data Distribution by Temperature Treatment', fontsize=14)
ax.set_xticks(x)
ax.set_xticklabels(ranges)
#ax.legend(title='Temperature', loc='upper right')
ax.legend(title='Temperature', bbox_to_anchor=(0.5, 1.0),loc='upper center', ncol=7) 

ax.set_ylim(0, 10)
ax.grid(axis='y', alpha=0.3)

plt.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.1)

#plt.tight_layout()
plt.show()