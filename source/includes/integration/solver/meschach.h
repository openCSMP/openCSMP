#ifndef ANSI_PORT_OF_MESCHACH_H
#define ANSI_PORT_OF_MESCHACH_H

#include "CSMP_definitions.h"

#ifdef CSMP_PC_WINDOWS_NT_CODE_WARRIOR
  #define CSMP_WITH_MESCHACH
  #define _MESCHACH
#endif

#ifdef CSMP_WITH_MESCHACH

#include "iter.h"
#include "sparse.h"
#include "sparse2.h"
#include "SparseMatrix.h"

void convert_CSP_SparseMatrix_to_SPMAT( const csmp::SparseMatrix& mcsp, SPMAT* spmat );

#endif//CSMP_WITH_MESCHACH

#endif
