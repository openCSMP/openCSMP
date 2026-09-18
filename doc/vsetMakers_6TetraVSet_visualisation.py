import rhinoscriptsyntax as rs
import Rhino
import random
import scriptcontext as sc

def create_tetrahedron_mesh(nodes, element_nodes, color):
    """Create a mesh for a tetrahedron from node indices and coordinates."""

    # Extract the four node coordinates for the tetrahedron
    pts = [nodes[i] for i in element_nodes]

    # Build vertex list
    vertices = [(pt[0], pt[1], pt[2]) for pt in pts]

    # Faces (zero-based indices into vertices list)
    faces = [
        (1, 2, 3),  # Face0
        (0, 3, 2),  # Face1
        (0, 1, 3),  # Face2
        (0, 2, 1)   # Face3
    ]

    # Add mesh directly with rs
    mesh_id = rs.AddMesh(vertices, faces)

    # Assign color
    if mesh_id:
        rs.ObjectColor(mesh_id, color)
    return mesh_id


def visualize_tetra_vset():
    # Clear existing objects
    objs = rs.AllObjects()
    if objs: rs.DeleteObjects(objs)

    # Define node coordinates (from VSet)
    px = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0]
    py = [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0]
    pz = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0]

    nodes = [(px[i], py[i], pz[i]) for i in range(len(px))]

    # Define tetrahedral elements (from VSet 'plist')
    elements = [
        [0, 1, 3, 4],  # Element 0
        [4, 5, 1, 3],  # Element 1
        [4, 5, 3, 7],  # Element 2
        [1, 2, 3, 6],  # Element 3
        [3, 1, 6, 5],  # Element 4
        [5, 6, 3, 7]   # Element 5
    ]

    # Boundary flags (node labels)
    bflags = {
        0: "CNR1", 1: "CNR2", 2: "CNR3", 3: "CNR4",
        4: "CNR5", 5: "CNR6", 6: "CNR7", 7: "CNR8"
    }

    # Colors for tetrahedra
    colors = [(random.randint(0,255), random.randint(0,255), random.randint(0,255)) for _ in elements]

    # Add tetrahedra meshes
    for i, elem in enumerate(elements):
        create_tetrahedron_mesh(nodes, elem, colors[i])

        # Centroid for labeling
        centroid = [sum(nodes[idx][j] for idx in elem)/4.0 for j in range(3)]
        rs.AddTextDot(f"Tet {i}", centroid)

    # Add node points + labels
    for i, node in enumerate(nodes):
        rs.AddPoint(node)
        label = f"Node {i}: ({node[0]:.1f}, {node[1]:.1f}, {node[2]:.1f}) {bflags[i]}"
        label_pos = (node[0]+0.1, node[1]+0.1, node[2]+0.1)
        rs.AddTextDot(label, label_pos)

    # Set view to perspective + zoom extents
    view = sc.doc.Views.ActiveView
    if view:
        view.ActiveViewport.ChangeToPerspectiveProjection(True, 50.0)
        rs.ZoomExtents()


# Run
visualize_tetra_vset()