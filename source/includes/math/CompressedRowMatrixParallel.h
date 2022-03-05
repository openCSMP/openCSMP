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
    void RemoveHalo( int32_t& ); // TODO: check whether we do need int64_t here
    void CreatePartitions( uint32_t n_blocks, std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >&  ) const;

    // vector mu st be initialized with the partition ranges
    void ScanHaloOfPartition( uint32_t, std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >& ) const;

    // returns those rows as true in the corresponding boolean vector
    uint32_t RowsWithHaloElements( int32_t, std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >&  ) const;

    void ComputeVariableIndicesForBroadcast( uint32_t, const std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >&,
                                             std::vector<int32_t>&, std::vector<int32_t>&, std::vector<int32_t>&  ) const;

    void CreateReceiveList( uint32_t, const std  ::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >&,
                            std::vector<int32_t>&, std::vector<int32_t>&, std::vector<int32_t>&  ) const;

    void CreateIpts( int32_t, const std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >&,
                     std::vector<int32_t>&, std::vector<int32_t>&, std::vector<int32_t>&  ) const;

    void CreateIptr(int32_t, const std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >&,
                          std::vector<int32_t>&, std::vector<int32_t>&, std::vector<int32_t>&  ) const;


    void Renumber( uint32_t, const std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >&,
                        std::vector<int32_t>&, std::vector<int32_t>& ) const;
};
}
#endif // COMPRESSEDROWMATRIXPARALLEL_H
