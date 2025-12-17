import rhinoscriptsyntax as rs
import random

def create_hexahedron_mesh(nodes, element_nodes, color, offset=(0,0,0)):
    dx, dy, dz = offset
    vertices = [(nodes[i][0]+dx, nodes[i][1]+dy, nodes[i][2]+dz) for i in element_nodes]
    faces = [
        [0,3,2,1],  # bottom
        [0,1,5,4],  # front
        [1,2,6,5],  # right
        [2,3,7,6],  # back
        [0,4,7,3],  # left
        [4,5,6,7]   # top
    ]
    mesh_id = rs.AddMesh(vertices, faces)
    rs.ObjectColor(mesh_id, color)
    return mesh_id, vertices

def create_prism_mesh(nodes, element_nodes, color, offset=(0,0,0)):
    dx, dy, dz = offset
    vertices = [(nodes[i][0]+dx, nodes[i][1]+dy, nodes[i][2]+dz) for i in element_nodes]
    faces = [
        [0,2,1],    # bottom triangle
        [0,1,4,3],  # quad
        [1,2,5,4],  # quad
        [0,3,5,2],  # quad
        [3,4,5]     # top triangle
    ]
    mesh_id = rs.AddMesh(vertices, faces)
    rs.ObjectColor(mesh_id, color)
    return mesh_id, vertices

def visualize_prism_hexa_vset():
    rs.DeleteObjects(rs.AllObjects())

    # ---- Node coordinates ----
    iDim_i, iDim_j, iDim_k = 4,4,4
    iDim_k2 = iDim_i * iDim_j
    nodes = []
    for k in range(iDim_k):
        for j in range(iDim_j):
            for i in range(iDim_i):
                nodes.append((float(i), float(j), float(k)))
    nodes.append((1.5, 1.5, 1.5))  # barycenter node index 64

    # ---- Boundary flags ----
    bflags = ["NOT"] * len(nodes)
    corners = [
        (0,0,0,"CNR1"), (iDim_i-1,0,0,"CNR2"), (iDim_i-1,iDim_j-1,0,"CNR3"), (0,iDim_j-1,0,"CNR4"),
        (0,0,iDim_k-1,"CNR5"), (iDim_i-1,0,iDim_k-1,"CNR6"), (iDim_i-1,iDim_j-1,iDim_k-1,"CNR7"), (0,iDim_j-1,iDim_k-1,"CNR8")
    ]
    for i,j,k,label in corners:
        idx = k*iDim_k2 + j*iDim_j + i
        bflags[idx] = label

    # ---- Elements ----
    elements = []
    element_types = []

    # Hexahedrons & Prisms
    for k in range(iDim_k-1):
        for j in range(iDim_j-1):
            for i in range(iDim_i-1):
                if i==1 and j==1:
                    # Prism 1
                    n1 = iDim_k2*k + iDim_j*j + i
                    n2 = n1+1
                    n3 = iDim_k2*k + iDim_j*(j+1) + i + 1
                    n4 = iDim_k2*k + iDim_j*(j+1) + i
                    n5 = iDim_k2*(k+1) + iDim_j*j + i
                    n6 = n5+1
                    n7 = iDim_k2*(k+1) + iDim_j*(j+1) + i + 1
                    n8 = iDim_k2*(k+1) + iDim_j*(j+1) + i
                    elements.append([n5,n6,n8,n1,n2,n4])
                    element_types.append("Prism")
                    elements.append([n3,n4,n2,n7,n8,n6])
                    element_types.append("Prism")
                else:
                    # Hexahedron
                    n0 = iDim_k2*k + iDim_j*j + i
                    element = [
                        n0,
                        n0+1,
                        n0+iDim_j+1,
                        n0+iDim_j,
                        n0+iDim_k2,
                        n0+iDim_k2+1,
                        n0+iDim_k2+iDim_j+1,
                        n0+iDim_k2+iDim_j
                    ]
                    elements.append(element)
                    element_types.append("Hexahedron")

    # ---- Colors ----
    colors = [(random.randint(0,255), random.randint(0,255), random.randint(0,255)) for _ in elements]

    # ---- Add meshes with spacing ----
    spacing = 1.5
    for i, elem in enumerate(elements):
        # Compute grid offset
        offset = ( (i % 4)*spacing, (i//4 % 4)*spacing, (i//16)*spacing )

        if element_types[i]=="Hexahedron":
            mesh_id, verts = create_hexahedron_mesh(nodes, elem, colors[i], offset)
        else:
            mesh_id, verts = create_prism_mesh(nodes, elem, colors[i], offset)

        # Centroid label
        cx = sum(v[0] for v in verts)/len(verts)
        cy = sum(v[1] for v in verts)/len(verts)
        cz = sum(v[2] for v in verts)/len(verts)
        rs.AddTextDot(f"{element_types[i][:3]} {i}", (cx, cy, cz))

        # Node labels
        for local_idx, node_idx in enumerate(elem):
            vx, vy, vz = verts[local_idx]
            label = f"Node {node_idx}: ({nodes[node_idx][0]:.1f},{nodes[node_idx][1]:.1f},{nodes[node_idx][2]:.1f}) {bflags[node_idx]}"
            rs.AddTextDot(label, (vx+0.1, vy+0.1, vz+0.1))
            rs.AddPoint((vx,vy,vz))

    rs.ViewDisplayMode("Shaded")
    rs.ZoomExtents()


# Run
if __name__=="__main__":
    visualize_prism_hexa_vset()