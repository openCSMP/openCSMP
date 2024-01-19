#ifndef PARTITION_VSET_H
#define PARTITION_VSET_H

#include "CSP_VSet.h"
#include "VSetConnectivity.h"

namespace csp {

template<typename fT, stl_index dim>
bool partitionVSet( VSet<fT,dim>&,                                // input VSet (not constant because might be rectified)
                    const std::vector<std::vector<stl_index> >&,  // target regions 
                    std::vector<std::set<std::pair<stl_index,int> > >&,
                    std::vector<VSet<csp_float,dim> >&,	// to store sub-VSets
                    VSetConnectivity<fT,dim>&,
                    bool );

// uses information from file to initialize element types (see below) 
template<typename fT, stl_index dim>
void createPfverts( VSet<fT,dim>& vset, bool isoparametric, const char* etypes_file="ICEM_element_types.txt" ); 

} // end namespace csp

#endif

/*
'ICEM_element_types.txt' 'plist' specs. for all CSP elements (SKM15/7/2005)
element-type  nodes  segments  faces, nodes-per-face followed by node numbers for each face
TRI_3  3 3 3
2  1 2
2  2 0
2  0 1
TRI_3_X 4 3 3
2  1 2
2  2 0
2  0 1
TRI_6  6 3 3  
3  1 4 2
3  2 5 0
3  0 3 1
TRI_6_X  7 3 3
3  1 4 2
3  2 5 0
3  0 3 1
TETRA_4  4 6 4
3  1 2 3
3  0 3 2
3  0 1 3
3  0 2 1
TETRA_10  10 6 4
6  1 5 2 9 3 8
6  0 7 3 9 2 6
6  0 4 1 8 3 7
6  0 6 2 5 1 4 
PYRA_5  5 8 5
4  0 3 2 1
3  0 1 4
3  1 2 4
3  2 3 4
3  0 4 3
PYRA_13  13 8 5
7  0 8 3 7 2 6 1 5
6  0 5 1 10 4 9
6  1 6 2 11 4 10
6  2 7 3 12 4 11
6  0 9 4 12 3 8 
PENTA_6  6 9 5
3  0 2 1 
4  0 1 4 3
4  1 2 5 4
4  0 3 5 2
3  3 4 5
PENTA_15  15 9 5
6  0 8 2 7 1 6
8  0 6 1 10 4 12 3 9
8  1 7 2 11 5 13 4 10
8  0 9 3 14 5 11 2 8
6  3 12 4 13 5 14
QUAD_4  4 4 4
2  0 1
2  1 2
2  2 3
2  3 0
QUAD_4_X  5 4 4
2  0 1
2  1 2
2  2 3
2  3 0
QUAD_8  8 4 4
3  0 4 1
3  1 5 2
3  2 6 3
3  3 7 0
QUAD_8_X  9 4 4 
3  0 4 1
3  1 5 2
3  2 6 3
3  3 7 0
HEXA_8  8 12 6
4  0 3 2 1
4  0 1 5 4
4  1 2 6 5
4  2 3 7 6
4  0 4 7 3
4  5 6 7 4
HEXA_20  20 12 6
8  0 11 3 13 2 9 1 8
8  0 8 1 10 5 16 4 12
8  1 9 2 14 6 17 5 10
8  2 13 3 15 7 18 6 14
8  0 12 4 19 7 15 3 11
8  4 16 5 17 6 18 7 19
BAR_2  2 1 1
2  0 1
BAR_3  3 1 1
3  0 1 2
END
*/
