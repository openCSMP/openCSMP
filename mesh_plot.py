import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

# Node coordinates
px = [0, 2, 5, 6, 0, 2, 4.0, 6, 1.8, 2.8, 3.5, 4.7, 6, 0, 2, 4.2, 6, 0, 2.5, 4.3, 4.8, 5.0]
py = [5, 5, 5, 5, 4, 4, 3.5, 3, 3.2, 3.0, 2.5, 2.5, 2, 2.5, 1.5, 1.0, 0, 0, 0, 0, 3.0, 4.5]

# Triangular elements (node lists)
tri_elements = [
    [1,5,6], [0,4,5], [0,5,1], [5,9,6],
    [2,21,3], [3,21,7],
    [4,8,5], [5,8,9], [6,9,10],
    [4,13,8], [8,13,14], [8,14,10], [8,10,9],
    [10,15,11], [11,15,12], [12,15,16],
    [13,17,14], [14,17,18], [15,18,19], [15,19,16],
]

# Quadrilateral elements split into 2 triangles each for plotting
quad_elements = [
    [1,6,21,2], [6,10,11,20], [6,20,7,21],
    [7,20,11,12], [10,19,18,15],
]

fig, ax = plt.subplots(figsize=(8, 8))

# Plot triangular elements
for tri in tri_elements:
    xs = [px[n] for n in tri] + [px[tri[0]]]
    ys = [py[n] for n in tri] + [py[tri[0]]]
    ax.plot(xs, ys, 'k-', linewidth=0.8)

# Plot quad elements
for quad in quad_elements:
    xs = [px[n] for n in quad] + [px[quad[0]]]
    ys = [py[n] for n in quad] + [py[quad[0]]]
    ax.plot(xs, ys, 'k-', linewidth=0.8)

# Label nodes
for i, (x, y) in enumerate(zip(px, py)):
    ax.plot(x, y, 'ko', markersize=4)
    ax.annotate(str(i), (x, y), textcoords="offset points",
                xytext=(5, 5), fontsize=8, color='blue')

# Label element centroids
all_elements = tri_elements + [[q[0],q[1],q[2]] for q in quad_elements]
elem_labels = list(range(len(all_elements)))
for i, el in enumerate(all_elements):
    cx = np.mean([px[n] for n in el])
    cy = np.mean([py[n] for n in el])
    ax.annotate(f'({i})', (cx, cy), fontsize=7, color='red',
                ha='center', va='center')

# Boundary box
ax.set_xlim(-0.3, 6.5)
ax.set_ylim(-0.3, 5.5)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Mesh — Nodes (blue) and Elements (red)')
ax.set_aspect('equal')
ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig('mesh.png', dpi=150)
plt.show()

