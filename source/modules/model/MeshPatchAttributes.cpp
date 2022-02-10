//
//  MeshPatchAttributes.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 23/5/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "MeshPatchAttributes.h"

using namespace std;

namespace csmp {

MeshPatchAttributes::MeshPatchAttributes( size_t cells, ELEMENT_DIMENSION dim )
 : cells_(cells), cell_dimension_(dim)
 {
 }
 

void MeshPatchAttributes::Cells( size_t n )
 {
    cells_ = n;
 }
 
 
size_t MeshPatchAttributes::Cells() const
 {
    return cells_;
 }
 
 
void MeshPatchAttributes::Geometry( ELEMENT_DIMENSION dim )
 {
    cell_dimension_ = dim;
 }
 
 
ELEMENT_DIMENSION MeshPatchAttributes::Geometry() const
 {
    return cell_dimension_;
 }


} // end csmp
