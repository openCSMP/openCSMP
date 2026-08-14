/**
    FEM Conventions 
    
     1. Outward Normal (Standard FEA)
       In computational geometry and finite element analysis (FEA), 
       the goal is to ensure the normal vector (the vector perpendicular 
       to the face) points outward from the closed volume.
       To achieve this outward normal, you must use the Right-Hand Rule (RHR):
       Curl your right-hand fingers in the order the nodes 
       are listed (e.g., $N_1 \to N_2 \to N_3$).
       Your right thumb points in the direction of the surface normal.
       For a closed volume, if you follow the node order around the face's perimeter, 
       the nodes must appear to be ordered Counter-Clockwise (CCW) when viewed along 
       the direction of the outward normal!
*/

// TODO: some face-node numbers are inconsistent in that the edges of adjacent faces are not correctly aligned. The normal orientations are all correct but edge-based calcs may play up


// Standard hexahedron element definition for computational analysis.
// Adjacency is consistent: shared edges between faces are traversed in opposite directions.
const vector<vector<size_t>> hexaFaceNodes_corrected = {
    // 1. Bottom: Edge winding (0->1, 1->2, 2->3, 3->0)
    {0, 3, 2, 1},

    // 2. Front: Must start with 1 or 0, traverse shared edge 1 -> 0
    {1, 5, 4, 0},

    // 3. Right: Must start with 2 or 1, traverse shared edge 2 -> 1
    {2, 6, 5, 1},

    // 4. Back: Must start with 3 or 2, traverse shared edge 3 -> 2
    {3, 7, 6, 2},

    // 5. Left: Must start with 0 or 3, traverse shared edge 0 -> 3
    {0, 4, 7, 3},

    // 6. Top: Node winding 4 -> 5 -> 6 -> 7
    {4, 5, 6, 7}
};


// Corrected triangular prism (wedge) face node definition.
// All faces are consistently oriented (e.g., outward-facing normals).
const vector<vector<size_t>> prismFaceNodes_corrected = {
    // 1. Face 0: Bottom Triangle (CCW)
    {0, 2, 1},

    // 2. Face 1: Side Quad (Reversed from original to ensure 0->1 bottom edge and 4->3 top edge)
    // Original: {0, 3, 4, 1}
    {0, 1, 4, 3},

    // 3. Face 2: Side Quad (Reversed from original to ensure 1->2 bottom edge and 5->4 top edge)
    // Original: {1, 4, 5, 2}
    {1, 2, 5, 4},

    // 4. Face 3: Side Quad (This one was already correct)
    {0, 3, 5, 2},

    // 5. Face 4: Top Triangle (CCW)
    {3, 4, 5}
};

// PYRAMID
   switch (face_id) {
        
        // Base Face 4: {0, 1, 2, 3} (CW, Normal -Z)
        case 4: return vector<uint32_t>{0, 1, 2, 3}; 
        
        // Face 0: Connects base edge 3->0 (opposite of base 0->3). Vertical edges 0->4, 4->3.
        case 0: return vector<uint32_t>{0, 3, 4}; 

        // Face 1: Connects base edge 0->1 (opposite of base 1->0). Vertical edges 1->4, 4->0.
        case 1: return vector<uint32_t>{1, 0, 4}; 

        // Face 2: Connects base edge 1->2 (opposite of base 2->1). Vertical edges 2->4, 4->1.
        case 2: return vector<uint32_t>{2, 1, 4}; 

        // Face 3: Connects base edge 2->3 (opposite of base 3->2). Vertical edges 3->4, 4->2.
        // TODO: this currently wrong in CSMP!
        case 3: return vector<uint32_t>{3, 2, 4}; 
