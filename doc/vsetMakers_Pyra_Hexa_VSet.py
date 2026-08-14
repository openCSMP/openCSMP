import rhinoscriptsyntax as rs
import Rhino
import random

def create_hexahedron_mesh(nodes, element_nodes, color):
    """Create a mesh for a hexahedron from node indices and coordinates."""
    pts = [nodes[i] for i in element_nodes]
    
    # Convert vertices to plain tuples
    vertices = [(pt[0], pt[1], pt[2]) for pt in pts]
    
    # Define the six quadrilateral faces of the hexahedron
    faces = [
        (0, 3, 2, 1),  # BOTTOM
        (0, 1, 5, 4),  # FRONT
        (1, 2, 6, 5),  # RIGHT
        (2, 3, 7, 6),  # BACK
        (0, 4, 7, 3),  # LEFT
        (4, 5, 6, 7)   # TOP
    ]
    
    mesh_id = rs.AddMesh(vertices, faces)
    rs.ObjectColor(mesh_id, color)
    return mesh_id

def create_pyramid_mesh(nodes, element_nodes, color):
    """Create a mesh for a pyramid from node indices and coordinates."""
    pts = [nodes[i] for i in element_nodes]
    
    # Convert vertices to plain tuples
    vertices = [(pt[0], pt[1], pt[2]) for pt in pts]
    
    # Faces: base is quad, others are triangles
    faces = [
        (0, 1, 2, 3),  # Base
        (0, 1, 4),     # Front
        (1, 2, 4),     # Right
        (2, 3, 4),     # Back
        (0, 4, 3)      # Left
    ]
    
    mesh_id = rs.AddMesh(vertices, faces)
    rs.ObjectColor(mesh_id, color)
    return mesh_id

def visualize_pyramid_hexa_vset():
    rs.DeleteObjects(rs.AllObjects())
    
    # Define node coordinates
    nodes = []
    iDim_i, iDim_j, iDim_k = 4, 4, 4
    iDim_k2 = iDim_k * iDim_j
    for k in range(iDim_k):
        for j in range(iDim_j):
            for i in range(iDim_i):
                nodes.append((float(i), float(j), float(k)))
    nodes.append((1.5, 1.5, 1.5))  # barycenter
    
    # Boundary flags (keep your original logic here if needed)
    bflags = ["NOT"] * 65
    
    # Elements (hexahedrons)
    elements = []
    iDim_km1_2 = (iDim_k - 1) * (iDim_j - 1)
    for k in range(iDim_k - 1):
        for j in range(iDim_j - 1):
            for i in range(iDim_i - 1):
                if i == 1 and j == 1 and k == 1:
                    continue
                element = [
                    k * iDim_k2 + j * iDim_j + i,
                    k * iDim_k2 + j * iDim_j + i + 1,
                    k * iDim_k2 + (j + 1) * iDim_j + i + 1,
                    k * iDim_k2 + (j + 1) * iDim_j + i,
                    (k + 1) * iDim_k2 + j * iDim_j + i,
                    (k + 1) * iDim_k2 + j * iDim_j + i + 1,
                    (k + 1) * iDim_k2 + (j + 1) * iDim_j + i + 1,
                    (k + 1) * iDim_k2 + (j + 1) * iDim_j + i
                ]
                elements.append(element)
    
    # Pyramid elements
    elements.append([21, 22, 26, 25, 64])
    elements.append([37, 38, 22, 21, 64])
    elements.append([22, 38, 42, 26, 64])
    elements.append([25, 26, 42, 41, 64])
    elements.append([37, 21, 25, 41, 64])
    elements.append([38, 37, 41, 42, 64])
    
    # Colors
    colors = [(random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
              for _ in range(len(elements))]
    
    # Add meshes
    for i, elem in enumerate(elements):
        if len(elem) == 8:  # Hexahedron
            create_hexahedron_mesh(nodes, elem, colors[i])
        else:               # Pyramid
            create_pyramid_mesh(nodes, elem, colors[i])
        
        # Centroid for labeling
        centroid = [0.0, 0.0, 0.0]
        for node_idx in elem:
            centroid[0] += nodes[node_idx][0] / len(elem)
            centroid[1] += nodes[node_idx][1] / len(elem)
            centroid[2] += nodes[node_idx][2] / len(elem)
        rs.AddTextDot(f"{'Hex' if len(elem)==8 else 'Pyr'} {i}", centroid)
    
    # Add nodes with labels
    for i, node in enumerate(nodes):
        rs.AddPoint(node)
        label = f"Node {i}: ({node[0]:.1f}, {node[1]:.1f}, {node[2]:.1f}) {bflags[i]}"
        label_pos = (node[0] + 0.1, node[1] + 0.1, node[2] + 0.1)
        rs.AddTextDot(label, label_pos)
    
    rs.ViewDisplayMode("Shaded")
    rs.ZoomExtents()

visualize_pyramid_hexa_vset()