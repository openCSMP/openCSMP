#include "CompressedRowMatrixParallel.h"
#include "SparseMatrix.h"
#include "STL_utilities.h"
#include "Exception.h"

using namespace std;

namespace csmp {

CompressedRowMatrixParallel::CompressedRowMatrixParallel()
{
}

/** Removes outerhalo node entries from the matrix

@param  nrhalo = total number of rows which need to be removed
*/
void CompressedRowMatrixParallel::RemoveHalo(int32_t& nrhalo)
 {
   // This method should be used when doing a mesh partitioning
   // rather than a matrix partitioning. In the case of a mesh partitioning
   // the rows of the outerhalo nodes need to be removed from the global
   // solution matrix. This needs to be done in SAMGp_Solver:SolveMatrixEquation()
   // Outer halo nodes have the highest node ID's, so rows (Rows()-nrhalo)
   // till Rows() need to be removed. The columns should be kept as is.
   int32_t new_nna = ia[ia.size()-nrhalo-1]-1;

   ja.resize(new_nna);            vector<int32_t>( ja ).swap( ja );
   a.resize (new_nna);            vector<double>( a ).swap( a );
   ia.resize(ia.size()-nrhalo);   vector<int32_t>( ia ).swap( ia );
 }

/**

Break matrix into even row ranges and return them into argument vector.

@param ranges The ranges are stored in the first pair of the vector.

*/
void CompressedRowMatrixParallel::CreatePartitions( uint32_t n_parts, vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >&  ranges ) const
 {
      if  ( !ranges.empty() )  ranges.erase( ranges.begin(), ranges.end() );

       uint32_t nrow= static_cast<uint32_t>(floor(static_cast<double>((ia.size()-1)/n_parts)));
       // N1 partitions which store (nrow) rows
       // (n_parts - N1) partitions which store (nrow+1) rows
       uint32_t N1 = n_parts - (ia.size() - nrow * n_parts) + 1U; // Number of partitions with (nrow) rows

       ranges.resize(n_parts);
       vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >( ranges ).swap( ranges );
       int32_t p(0);

       for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::iterator
             it=ranges.begin(); it!=ranges.end(); it++, p++)
         {
              if (p < N1) { 			// Partition gets size nrow
                 (*it).first.first = p*nrow;
                 (*it).first.second = (p+1U)*nrow - 1U;
                 }
              else if (p >= N1) { 			// Partition gets size nrow+1
                 (*it).first.first = N1*nrow+(p-N1)*(nrow+1U);
                 (*it).first.second = N1*nrow+(p-N1)*(nrow+1U)+nrow;
                 }
               else cout <<"\nCompressedRowMatrix::CreatePartitions: partition setting failed.\n";
          }
 }


/**

Looks along the rows of the current partition for non-zero elements. Where
such elements are found in column ranges corresponding to the other
processors, the row entry of the boolean vector (one for each processor)
is set to true else false. The method also takes care of the initialization
of the boolean vectors.

@section implementation Implementation

While all boolean vectors are set to the same size = nrows of current
bock, the column ranges of the partitions are used in the scanning.
*/
void CompressedRowMatrixParallel::ScanHaloOfPartition( uint32_t n_block, vector<pair<pair<uint32_t,uint32_t>,std::vector<bool> > >&  partition_vector  ) const
 {
    // for the current block the boolean vectors (one for each connected processor) are initialized
    for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::iterator
            it=partition_vector.begin(); it!=partition_vector.end(); it++ ) {
            // make zero-initialized storage for flags that signify the rows in the blocks that have elements
            (*it).second.resize(  partition_vector[n_block].first.second - partition_vector[n_block].first.first + 1U );
            // set vector to false
            fill( (*it).second.begin(), (*it).second.end(), false );
       }

    size_t i,j;
    // loop over the rows of interest initializing boolean vectors
    for ( i=partition_vector[n_block].first.first; i<= partition_vector[n_block].first.second; i++ ){
      if ( i < ia.size()-1U ) 							// if not last row
        for ( j=ia[i]-1U; j<ia[i+1U]-1U; j++ ){       // loop over colum
          for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::iterator
              it=partition_vector.begin(); it!=partition_vector.end(); it++ ) {
               // for those row entries that fall into the column range of the iterated slice set boolean vector to true
               if ( ja[j]-1U >= (*it).first.first  and  ja[j]-1U <= (*it).first.second )
                 (*it).second[i-partition_vector[n_block].first.first] = true;
               }
     }
      else
        for ( j=ia[i]-1U; j<ja.size(); j++ ){
          for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::iterator
              it=partition_vector.begin(); it!=partition_vector.end(); it++ ) {
               if ( ja[j]-1U >= (*it).first.first  and  ja[j]-1U <= (*it).first.second )
                 (*it).second[i-partition_vector[n_block].first.first] = true;
               }
     }
   }

 } // end ScanHaloOfPartition




/**

The SAMGp_CommunicationData dealing with receiving data are determined for the curren processor
*/
void CompressedRowMatrixParallel::CreateReceiveList( uint32_t n_block,
                                             const vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >&  partition_vector,
                                             vector<int32_t>& irankrec,
                                             vector<int32_t>& ireclists,
                                             vector<int32_t>& iptr ) const
{
    assert( n_block < partition_vector.size() );

    set<int32_t> receivelist;
    size_t quantity(0U);
    iptr.resize(1U);
    vector<int32_t>( iptr ).swap( iptr );
    iptr[0]=1U;

    for( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::const_iterator
           it=partition_vector.begin(); it!=partition_vector.end(); it++)
        // block is not the partition
        if((*it)!=partition_vector[n_block])
         {
           // check if this block contains non-zeroes
           for ( vector<bool>::const_iterator z=(*it).second.begin(); z!=(*it).second.end(); z++ )
             // Loop over rows of partition
             for( size_t i=partition_vector[n_block].first.first; i<=partition_vector[n_block].first.second; i++ )
               // Loop over elements of row
               for ( size_t j=ia[i]-1U; j<ia[i+1U]-1U;j++ )
                 // If within block
                 if(ja[j]-1U >= (*it).first.first and ja[j]-1U <= (*it).first.second)
                   // Store column number in receivelist
                   receivelist.insert( ja[j] );
           // If new entries have been added to list
           if (receivelist.size()>quantity)
            {
              // Store number of entries + 1 in iptr
              iptr.push_back( receivelist.size()+1 );
              quantity = receivelist.size();
            }

          }
    ireclists.assign(receivelist.begin(), receivelist.end());

    set<int32_t>  processors;
    int32_t       n_proc(0);

    // calculate the size of the receivelist vector that will be created
    for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::const_iterator
          it=partition_vector.begin(); it!=partition_vector.end(); it++, n_proc++ )
       for ( vector<bool>::const_iterator  bt=(*it).second.begin(); bt!=(*it).second.end(); bt++ )
         if ( *bt == true ) {
              // register the processors which deals with the columns in question as the receiving processors
              if ( (*it) != partition_vector[n_block] ) processors.insert( n_proc );
              // count true elements in this column range
           }

    irankrec.assign( processors.begin(), processors.end() );
    assert(iptr.size() == (irankrec.size()+1U));

}


/*
   cout <<"Irec size = "<< (*rit.first).second.size() << endl;
   for (vector<int32_t>::const_iterator i=(*rit.first).second.begin(); i!=(*rit.first).second.end();i++)
      cout <<"Reclist: "<< *i << endl;
   cout <<"\nLooping over block "<< dummy << endl;
   cout <<"Min, Max: "<< (*it).first.first <<", "<< (*it).first.second << endl;
   cout <<"Looping over row " << i << endl;
   cout <<"Found column "<< ja[j] << endl;



*/



void CompressedRowMatrixParallel::ComputeVariableIndicesForBroadcast( uint32_t n_block,
                                                              const vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >&  partition_vector,
                                                              vector<int32_t>& iranksnd, vector<int32_t>& isndlists, vector<int32_t>& ipts ) const
{
    assert( n_block < partition_vector.size() );
    uint32_t n_proc(0U);
    ipts.resize(1U);   vector<int32_t>( ipts ).swap( ipts );
    ipts[0]=1U;
    set<int32_t> proc;
    size_t counter(1U);
    int32_t  temp2(0);


    // Loop over partition_vector
    for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::const_iterator
          it=partition_vector.begin(); it!=partition_vector.end(); it++, n_proc++)
        if ((*it)!=partition_vector[n_block])
          {
            set<int32_t> sendlist;
            vector<int32_t> temp;
            // loop over entries in matrix
            if ( (*it).first.second+1U != ia.size() ) 	// if not last column
              for( int32_t j=ia[(*it).first.first]; j!=ia[(*it).first.second+1U]; j++ )
                {
                  // if column is within block
                  if( (ja[j-1U]-1U)>=partition_vector[n_block].first.first and (ja[j-1U]-1U)<=partition_vector[n_block].first.second)
                    {
                      // insert processor number in proc
                      proc.insert(n_proc);
                      // store column number in set
                      sendlist.insert(ja[j-1U]);
                    }
                }
            // if entries in this part of matrix, increase ipts
            if (!sendlist.empty())
              {
                counter += sendlist.size();
                ipts.push_back(counter);
              }
            temp2=sendlist.size();
            temp.assign(sendlist.begin(), sendlist.end());
            // loop over set and store values in isndlist
            for( uint32_t i=0U; i<temp.size(); i++)
              isndlists.push_back(temp[i]);
          }

    iranksnd.assign( proc.begin(), proc.end() );
    assert(ipts.size() == (iranksnd.size()+1U));

}
/*
cout<<"\Partition"<<n_block<<" ranging from "<<partition_vector[n_block].first.first<<" to "<<partition_vector[n_block].first.second<<  endl;
cout<<"\nBLOCK"<<n_proc<<" ranging from "<<(*it).first.first<<" to "<<(*it).first.second<<  endl;
cout <<"Entry # j = "<<j-1U<<" wih corresponding column = "<<ja[j-1U]-1U<< endl;
cout <<"Column in partition area...storing "<<ja[j-1U]<<" in set "<< endl;
cout <<"sendlist size = "<< sendlist.size() << endl;
cout <<"STORING NEW IPTS ENTRY" << endl;

*/




void CompressedRowMatrixParallel::CreateIpts( int32_t n_block,
										const std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >& partition_vector,
										vector<int32_t>& iranksnd, vector<int32_t>& isndlists, vector<int32_t>& ipts) const

{

	// resize of the ipts vector
	ipts.resize(iranksnd.size()+1U);  vector<int32_t>( ipts ).swap( ipts );
	uint32_t proc(0), number(0);
	ipts[0]=1U;

	//loop over partition vector
	for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::const_iterator
						 it=partition_vector.begin(); it!=partition_vector.end(); it++ )
			// if columns outside the partition
			if( (*it) != partition_vector[n_block] ){
						 proc++;
						 // loop over boolean vector
						 for ( int32_t i=0; i<static_cast<int32_t>((*it).second.size()); i++ )
							 // if there is non-zero element in the row in the column range of the looped over block
							 if ( (*it).second[i] == true ) number++;
						 ipts[proc]=number+1;
				 }

		ipts[iranksnd.size()] = isndlists.size()+1;

}


void CompressedRowMatrixParallel::CreateIptr( int32_t n_block,
										const std::vector<std::pair<std::pair<uint32_t,uint32_t>,std::vector<bool> > >& partition_vector,
										vector<int32_t>& irankrec, vector<int32_t>& ireclists, vector<int32_t>& iptr) const

{

	// resize of the iptr vector
	iptr.resize(irankrec.size()+1U);  vector<int32_t>( iptr ).swap( iptr );
		uint32_t proc(0), number(1U);
		iptr[0]=number;

	//loop over partition vector
	for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::const_iterator
						 it=partition_vector.begin(); it!=partition_vector.end(); it++ )
			// if columns outside the partition
			if( (*it) != partition_vector[n_block] ){
					 proc++;
					 // loop over boolean vector
					 for ( int32_t i=0; i<static_cast<int32_t>((*it).second.size()); i++ )
						 // if there is non-zero element in the row in the column range of the looped over block
						 if ( (*it).second[i] == true ) number++;
					 iptr[proc]=number;
				 }

		iptr[irankrec.size()] = ireclists.size()+1;
}
/**

Put into boolean vector of the current partition true values where there
are halo elements in the corresponding row.

@param partition_vector the number of true values in this vector.

*/
void CompressedRowMatrixParallel::Renumber( uint32_t n_block, const vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >& partition_vector,
																		vector<int32_t>& ireclists, vector<int32_t>& isndlists) const
{
	// isndlists is straightforward:
	for (vector<int32_t>::iterator it=isndlists.begin();it!=isndlists.end();it++)
		 (*it) -= partition_vector[n_block].first.first;

	uint32_t counter=partition_vector[n_block].first.second - partition_vector[n_block].first.first+1U;
	// loop over ireclists
	for (vector<int32_t>::iterator it=ireclists.begin();it!=ireclists.end();it++)
		{
			 // ireclists should never have entries with the same value
			 counter++;
			 (*it) = counter;
		}
}

/**

Put into boolean vector of the current partition true values where there
are halo elements in the corresponding row.

@return the number of true values in this vector.
*/
uint32_t CompressedRowMatrixParallel::RowsWithHaloElements( int32_t n_block,
                                                  vector<pair<pair<uint32_t,uint32_t>,vector<bool> > > & partition_vector  ) const
 {
    // falsify boolean vector of current partition
    fill( partition_vector[n_block].second.begin(), partition_vector[n_block].second.end(), false );

    size_t i,j;
    // loop over the rows of interest initializing boolean vectors
    for ( i=partition_vector[n_block].first.first; i<= partition_vector[n_block].first.second; i++ ){
      if ( i < ia.size()-1U ){ 							// if not last row
        for ( j=ia[i]-1U; j<ia[i+1U]-1U; j++ )          // loop over colum
          for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::iterator
              it=partition_vector.begin(); it!=partition_vector.end(); it++ )
               // for those row entries that fall into the column range of the iterated slice set boolean vector to true
               if ( ja[j]-1U < partition_vector[n_block].first.first  or  ja[j]-1U > partition_vector[n_block].first.second ) // if outside current partition
                 if ( (*it).second[i-partition_vector[n_block].first.first] ) {
                    partition_vector[n_block].second[i-partition_vector[n_block].first.first] = true;
                  }
      }
      else {
        for ( j=ia[i]-1U; j<ja.size(); j++ )
          for ( vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >::iterator
              it=partition_vector.begin(); it!=partition_vector.end(); it++ )
               // for those row entries that fall into the column range of the iterated slice set boolean vector to true
               if ( ja[j]-1U < partition_vector[n_block].first.first  or  ja[j]-1U > partition_vector[n_block].first.second ) // if outside current partition
                 if ( (*it).second[i-partition_vector[n_block].first.first] ) {
                    partition_vector[n_block].second[i-partition_vector[n_block].first.first] = true;
                  }
      }

    }

      return count(partition_vector[n_block].second.begin(), partition_vector[n_block].second.end(),true);
   
   }// end RowsWithHaloElements

/*
cout<<"RANGE: "<< partition_vector[n_block].first.first <<": "<< partition_vector[n_block].first.second << endl;
cout <<"\nia["<<i<<"] = "<<ia[i]<< endl;
cout <<"ja["<<j<<"]-1 = "<<ja[j]-1U<< endl;
cout<<"Outside current block"<<endl;
cout<<"Non-zero-element"<< endl;
*/

} //end namespace csmp
