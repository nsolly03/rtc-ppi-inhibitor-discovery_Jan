import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
import seaborn as sns

# Set style for R-like appearance
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Data
data = {
    'Week': list(range(1, 53)),
    '2024': [17749, 20114, 29197, 35215, 31774, 26793, 24121, 21478, 20826, 19519,
             18302, 17546, 15360, 14116, 14537, 13506, 14207, 15188, 15700, 16061,
             14569, 14308, 13924, 13816, 13523, 13135, 10476, 10909, 10261, 10517,
             10613, 10622, 10143, 10730, 10299, 10284, 10393, 11208, 12386, 12773,
             12779, 12986, 13053, 13434, 13595, 15336, 16618, 16267, 16991, 16038,
             15134, 13383],
    '2025': [12901, 12829, 14210, 13857, 13815, 14080, 14111, 14062, 14144, 13345,
             13274, 13804, 14014, 13714, 13687, 16032, 11446, 12173, 12238, 12838,
             13187, 12550, 12299, 13245, 12998, 12557, 10184, 10310, 10425, 10214,
             10317, 10375, 9236, 9358, 9671, 9759, 9381, 9318, 10675, 9971,
             9833, 10572, 11014, 11364, 12395, 17987, 22278, 24314, 22218, 23902,
             22079, 19319],
    '2026': [18745, 19037, 22677, 23821] + [np.nan] * 48
}

df = pd.DataFrame(data)

# Create figure with R-like aesthetics
fig, ax = plt.subplots(figsize=(14, 8), dpi=300)

# Define colors (R-like palette)
colors = {'2024': '#E41A1C', '2025': '#377EB8', '2026': '#4DAF4A'}
line_styles = {'2024': '-', '2025': '-', '2026': '--'}
markers = {'2024': 'o', '2025': 's', '2026': '^'}

# Plot each year
for year in ['2024', '2025', '2026']:
    mask = ~df[year].isna()
    ax.plot(df.loc[mask, 'Week'], df.loc[mask, year], 
            color=colors[year], 
            linewidth=2.5, 
            linestyle=line_styles[year],
            marker=markers[year],
            markersize=4,
            markevery=4,  # Show markers every 4th point
            label=year,
            alpha=0.9)

# Formatting
ax.set_xlabel('Epidemiological Week', fontsize=14, fontweight='bold', family='serif')
ax.set_ylabel('Number of Cases', fontsize=14, fontweight='bold', family='serif')
ax.set_title('Weekly Case Trends by Year (2024-2026)', 
             fontsize=16, fontweight='bold', family='serif', pad=20)

# Format y-axis with comma separator
def thousands(x, pos):
    return f'{int(x):,}'

ax.yaxis.set_major_formatter(FuncFormatter(thousands))

# Set axis limits and ticks
ax.set_xlim(0, 53)
ax.set_xticks(range(0, 53, 4))
ax.set_ylim(0, df[['2024', '2025', '2026']].max().max() * 1.1)

# Grid styling
ax.grid(True, linestyle='--', alpha=0.3, linewidth=0.8)
ax.set_axisbelow(True)

# Legend
legend = ax.legend(title='Year', 
                   loc='upper right', 
                   frameon=True, 
                   fancybox=True,
                   shadow=True,
                   fontsize=11,
                   title_fontsize=12)
legend.get_frame().set_alpha(0.95)
legend.get_frame().set_facecolor('white')
legend.get_frame().set_edgecolor('gray')

# Spine styling
for spine in ax.spines.values():
    spine.set_linewidth(1.2)
    spine.set_edgecolor('gray')

# Tick parameters
ax.tick_params(axis='both', which='major', labelsize=11, length=6, width=1.2)

# Add subtle background
ax.set_facecolor('#FAFAFA')
fig.patch.set_facecolor('white')

# Tight layout
plt.tight_layout()

# Save the figure
plt.savefig('/home/claude/epidemic_trend.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()

print("Professional epidemic trend chart created successfully!")
print("\nKey observations:")
print(f"- 2024: Peak at Week 4 ({df['2024'].max():,} cases)")
print(f"- 2025: Peak at Week 48 ({df['2025'].max():,} cases)")
print(f"- 2026: Data available through Week 4 (current: {df['2026'].dropna().iloc[-1]:,} cases)")