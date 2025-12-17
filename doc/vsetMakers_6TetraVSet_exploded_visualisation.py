import rhinoscriptsyntax as rs
import random
import scriptcontext as sc

def create_tetrahedron_mesh(nodes, element_nodes, color, offset=(0,0,0)):
    """Create a mesh for a tetrahedron from node indices and coordinates, with offset."""

    # Extract node coordinates + apply offset
    pts = [(nodes[i][0]+offset[0], nodes[i][1]+offset[1], nodes[i][2]+offset[2]) for i in element_nodes]

    # Build vertex list
    vertices = [(pt[0], pt[1], pt[2]) for pt in pts]

    # Faces (zero-based indices into vertices list)
    faces = [
        (1, 2, 3),  # Face0
        (0, 3, 2),  # Face1
        (0, 1, 3),  # Face2
        (0, 2, 1)   # Face3
    ]

    # Add mesh
    mesh_id = rs.AddMesh(vertices, faces)

    # Assign color
    if mesh_id:
        rs.ObjectColor(mesh_id, color)
    return mesh_id, vertices


def visualize_tetra_vset():
    # Clear existing objects
    objs = rs.AllObjects()
    if objs: rs.DeleteObjects(objs)

    # Define node coordinates (from VSet)
    px = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0]
    py = [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0]
    pz = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0]
    nodes = [(px[i], py[i], pz[i]) for i in range(len(px))]

    # Tetrahedral elements
    elements = [
        [0, 1, 3, 4],  # Element 0
        [4, 5, 1, 3],  # Element 1
        [4, 5, 3, 7],  # Element 2
        [1, 2, 3, 6],  # Element 3
        [3, 1, 6, 5],  # Element 4
        [5, 6, 3, 7]   # Element 5
    ]

    # Boundary flags (for node labeling)
    bflags = {
        0: "CNR1", 1: "CNR2", 2: "CNR3", 3: "CNR4",
        4: "CNR5", 5: "CNR6", 6: "CNR7", 7: "CNR8"
    }

    # Colors for tetrahedra
    colors = [(random.randint(0,255), random.randint(0,255), random.randint(0,255)) for _ in elements]

    # Spacing for separating tets
    spacing = 3.0

    for i, elem in enumerate(elements):
        # Arrange in 2D grid
        row = i // 3
        col = i % 3
        offset = (col*spacing, row*spacing, 0)

        mesh_id, verts = create_tetrahedron_mesh(nodes, elem, colors[i], offset)

        # Centroid
        centroid = [sum(v[j] for v in verts)/4.0 for j in range(3)]
        rs.AddTextDot(f"Tet {i}", centroid)

        # Add nodes + labels for each tetrahedron (local copy)
        for j, v in enumerate(verts):
            rs.AddPoint(v)
            orig_node_idx = elem[j]
            label = f"Node {orig_node_idx} {bflags[orig_node_idx]}"
            label_pos = (v[0]+0.1, v[1]+0.1, v[2]+0.1)
            rs.AddTextDot(label, label_pos)

    # Set view to perspective + zoom extents
    view = sc.doc.Views.ActiveView
    if view:
        view.ActiveViewport.ChangeToPerspectiveProjection(True, 50.0)
        rs.ZoomExtents()


# Run
visualize_tetra_vset()