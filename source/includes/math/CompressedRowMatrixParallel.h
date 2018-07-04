#ifndef COMPRESSED_ROW_MATRIX_PARALLEL_H
#define COMPRESSED_ROW_MATRIX_PARALLEL_H

#include "CompressedRowMatrix.h"
#include "SparseMatrix.h"

namespace csmp {

class CompressedRowMatrixParallel : public CompressedRowMatrix
{
public:
    CompressedRowMatrixParallel();
private:
    void RemoveHalo( int32& );
    void CreatePartitions( uint32 n_blocks, std::vector<std::pair<std::pair<uint32,uint32>,std::vector<bool> > >&  ) const;

    // vector mu st be initialized with the partition ranges
    void ScanHaloOfPartition( uint32, std::vector<std::pair<std::pair<uint32,uint32>,std::vector<bool> > >& ) const;

    // returns those rows as true in the corresponding boolean vector
    uint32 RowsWithHaloElements( int32, std::vector<std::pair<std::pair<uint32,uint32>,std::vector<bool> > >&  ) const;

    void ComputeVariableIndicesForBroadcast( uint32, const std::vector<std::pair<std::pair<uint32,uint32>,std::vector<bool> > >&,
                                             std::vector<int32>&, std::vector<int32>&, std::vector<int32>&  ) const;

    void CreateReceiveList( uint32, const std  ::vector<std::pair<std::pair<uint32,uint32>,std::vector<bool> > >&,
                            std::vector<int32>&, std::vector<int32>&, std::vector<int32>&  ) const;

    void CreateIpts( int32, const std::vector<std::pair<std::pair<uint32,uint32>,std::vector<bool> > >&,
                     std::vector<int32>&, std::vector<int32>&, std::vector<int32>&  ) const;

    void CreateIptr(int32, const std::vector<std::pair<std::pair<uint32,uint32>,std::vector<bool> > >&,
                          std::vector<int32>&, std::vector<int32>&, std::vector<int32>&  ) const;


    void Renumber( uint32, const std::vector<std::pair<std::pair<uint32,uint32>,std::vector<bool> > >&,
                        std::vector<int32>&, std::vector<int32>& ) const;
};
}
#endif // COMPRESSEDROWMATRIXPARALLEL_H
