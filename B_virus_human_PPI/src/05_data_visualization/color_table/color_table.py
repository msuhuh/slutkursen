import matplotlib.pyplot as plt

# Reverse the dictionary to have Enveloped/Non-enveloped at the bottom
reversed_color_table = {
    # Enveloped/Non-enveloped
    "Non-enveloped": "#FFE1DF",
    "Enveloped": "#AACCB6",
    # Gap
    # DNA/RNA
    "DNA": "#FFCCB6",
    "RNA": "#A2E1DB",
    # Gap
    # Virus types
    "dsDNA non-enveloped": "#B22222",
    "(-) ssRNA enveloped": "#228B22",
    "dsDNA enveloped": "#000000",
    "dsRNA non-enveloped": "#9932CC",
    "(+) ssRNA enveloped": "#8B4513",
    "(+) ssRNA non-enveloped": "#FF4500",
    "ssRNA circular enveloped": "#DAA520",
    "ssDNA non-enveloped": "#4682B4"
}

# Create the color table visualization with reversed order
fig, ax = plt.subplots(figsize=(8, 6))
y_offset = 0
for i, (label, color) in enumerate(reversed_color_table.items()):
    ax.add_patch(plt.Rectangle((0, y_offset), 1, 1, color=color))
    ax.text(1.2, y_offset + 0.5, label, va='center', ha='left', fontsize=10)
    y_offset += 1
    # Add small gap after Enveloped/Non-enveloped and DNA/RNA groups
    if label in {"Enveloped", "RNA"}:
        y_offset += 0.5

ax.set_xlim(0, 2)
ax.set_ylim(0, y_offset)
ax.set_aspect('equal')  # Ensure squares are perfectly square
ax.axis('off')
plt.show()



 