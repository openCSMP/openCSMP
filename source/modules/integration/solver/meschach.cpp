#include "meschach.h"
#include "CSMP_definitions.h"
#include <map>

#ifdef CSMP_WITH_MESCHACH
using namespace csmp;

void convert_CSP_SparseMatrix_to_SPMAT( const csmp::SparseMatrix& mcsp, SPMAT* spmat )
{
    //  giving the sparse matrix the correct size
    // this is done before:  spmat = sp_resize( spmat, static_cast<int>(mcsp.Rows()), static_cast<int>(mcsp.Cols()) );

    // iterator to row elements
    std::map<size_t,double>::const_iterator  it;

    // assigning the elements to the SPMAT structure
    std::cout <<"\nconvert_CSP_SparseMatrix_to_SPMAT: Converting the matrix..."<< std::endl;
    for ( unsigned int i=0; i<mcsp.Rows(); i++ )
        for ( it = mcsp.RowBegin(i); it != mcsp.RowEnd(i); it++ )
            // spmat->row[ static_cast<int>(i) ].elt[ (*it).first ].val = (*it).second;
            sp_set_val( spmat, static_cast<int>(i), static_cast<int>((*it).first), (*it).second );

    std::cout <<"\ndone."<< std::endl;
}

#endif
