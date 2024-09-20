#include <cmath>
#include <limits>
#include <cassert>
#include <fstream>
#include <utility>
#include "CompressedRowMatrix.h"
#include "SparseMatrix.h"
#include "Exception.h"
#include "CSMP_global_enumerations.h"
//#include "Region.h"
#include "Element.h"
#include "Face.h"
#include "NimbleRegion.h"
#include "NodeManifold.h"

//#define debug_sparsity_pattern

using namespace std;

namespace csmp {

/*! \file CSMP_mathUtilities.cpp */

/**
@addtogroup CSMPglobalFunctions
@{
*/

void print(  vector<pair<pair<uint32_t,uint32_t>,vector<bool> > >&  v )
 {
       cout <<"\nvector of off-diagonal elements:\n";
       uint32_t n(0);
       
         for (auto it=v.begin(); it != v.end(); it++ )
            {
                 cout <<"\nBlock "<< n++ <<" range: "<< (*it).first.first <<" - "<< (*it).first.second << endl;
                 cout <<"boolean vector (size="<< (*it).second.size() <<"):\n";
                 for ( auto i=(*it).second.begin(); i!=(*it).second.end(); i++ )
                   if ( *i ) cout <<" true  ";
                   else  cout <<"false ";
                   
                cout << endl << endl;
            }
  }
/**
@}
*/


/**
  generates sparsity pattern based on mesh connectivity information (i.e., neighbouring nodes)
  manifold nodes are also taken into account
  ia and ja are initialised, and values in a are initialised to zeros
  elimination of essential conditions is taken into account
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE>
void generateSparsityPatternEliminatingEssentialConditions( CompressedRowMatrix& G,
                                                            const std::map<Parameter,size_t>& test_operands,
                                                            std::vector<size_t>& DOF_indexes,
                                                            const ModelSubDomain<dim,CELLTYPE>& gref)

{
    assert(!test_operands.empty());
    if(test_operands.empty()) throw logic_error ("generateSparsityPatternEliminatingEssentialConditions: there is no test operand, nothing can be done");

    G.ia.push_back(0);
    int32_t index=0;
    double initial_value(0.);
    size_t num_nodes = gref.Nodes();

    for (const auto & test_operand : test_operands) {
      if (test_operand.first.key.place != NODE)
        throw csmp::Exception(ERROR, "generateSparsityPatternEliminatingEssentialConditions()", "only supporting nodal variables");

      csmp::Index prop_key = test_operand.first.key;
      uint32_t offset = static_cast<uint32_t>(test_operand.second);

      uint32_t variable_size{1U};
      if (prop_key.type == SCALAR) variable_size = 1;
      else if (prop_key.type == VECTOR) variable_size = dim;
      else if (prop_key.type == TENSOR) variable_size = dim*dim;
      else if (prop_key.type == ARRAY) variable_size = prop_key.dataDepth;
      else if (prop_key.type == FLAGGEDARRAY) variable_size = prop_key.dataDepth;
      else throw csmp::Exception(ERROR, "generateSparsityPatternEliminatingEssentialConditions()", "variable type not supported");

      //looping over interior nodes (no need to check manifold nodes)
      for (auto nit = gref.NodesBegin(); nit != gref.PerimeterNodesBegin(); nit++) {
        if((*nit)->Status(prop_key)==DIRICH) continue; //ignoring dirichlet nodes
        for (auto i{0U}; i < variable_size; i++) {
          set<uint32_t> node_indexes;
          //current node
          size_t idx = (*nit)->Idx();
          size_t pos = DOF_indexes[idx];
          if(pos != NULL_IDX)
            node_indexes.insert(static_cast<uint32_t>(pos) * variable_size + i + offset);

          //neighboring nodes of current node
          for (auto n{0U}; n < (*nit)->Neighbors(); n++) {
              auto nd = (*nit)->Neighbor(n);
              if(nd->Status(prop_key)==DIRICH) continue; //ignoring dirichlet nodes
              if(!gref.Contains(nd)) continue; //ignoring node outside domain
              idx = nd->Idx();
              pos = DOF_indexes[idx];
              if(pos != NULL_IDX)
                node_indexes.insert(static_cast<uint32_t>(pos) * variable_size + i + offset);
          }

          index += node_indexes.size();
          G.ia.push_back(index);
          for ( const auto& node_index : node_indexes ) {
          // TODO: the CRS only takes 'int32_t' for ja as required by samg
            assert( node_index < numeric_limits<int32_t>::max() );
            G.ja.push_back( static_cast<int32_t>(node_index) );
            G.a.push_back(initial_value);
          }
        }
      } //end looping over interior nodes

      //looping over perimeter nodes (need to check manifold nodes)
      for (auto nit = gref.PerimeterNodesBegin(); nit != gref.NodesEnd(); nit++) {
        if ((*nit)->Status(prop_key) == DIRICH) continue; //ignoring dirichlet nodes
        for (auto i{0U}; i < variable_size; i++) {
          set<size_t> node1_pos;
          //current node
          auto node1 = (*nit);
          size_t idx1 = node1->Idx();
          size_t pos1 = DOF_indexes[idx1];
          if (pos1 != NULL_IDX)
            node1_pos.insert(pos1 * variable_size + i + offset);

          //neighboring nodes of current node
          for (auto nb{0U}; nb < node1->Neighbors(); nb++) {
            auto nd = node1->Neighbor(nb);
            if(nd->Status(prop_key)==DIRICH) continue; //ignoring dirichlet nodes
            if(!gref.Contains(nd)) continue; //ignoring node outside domain
            auto idx = nd->Idx();
            auto pos = DOF_indexes[idx];
            if (pos != NULL_IDX)
              node1_pos.insert(pos * variable_size + i + offset);
          }

          //manifold nodes of current node
          if(node1->IsManifold()) {
            auto md = node1->Manifold();
            uint32_t branches = md->Branches();
            for (auto mn{0U}; mn < branches; mn++) {
              auto node = md->N(mn);
              if (node->Status(prop_key) == DIRICH) continue; //ignoring dirichlet nodes
              if (node == node1) continue; //ignoring current node itself
              if(!gref.Contains(node)) continue; //ignoring node outside domain
              //current manifold node
              auto idx = node->Idx();
              auto pos = DOF_indexes[idx];
              if (pos != NULL_IDX) {
                node1_pos.insert(pos * variable_size + i + offset);
              }
              //neighboring nodes of manifold node
              for (auto nn{0U}; nn < node->Neighbors(); nn++) {
                auto nd = node->Neighbor(nn);
                if (nd->Status(prop_key) == DIRICH) continue; //ignoring dirichlet nodes
                if(!gref.Contains(nd)) continue; //ignoring node outside domain
                idx = nd->Idx();
                pos = DOF_indexes[idx];
                if (pos != NULL_IDX)
                  node1_pos.insert(pos * variable_size + i + offset);
              }
            }
          }

          index += node1_pos.size();
          G.ia.push_back(index);
          for ( const auto& node_index : node1_pos ) {
            assert( node_index < numeric_limits<int32_t>::max() );
            G.ja.push_back(static_cast<uint32_t>(node_index));
            G.a.push_back(initial_value);
          }

        }
      } //end looping over perimeter nodes

      offset += variable_size * num_nodes;

    } //end looping over test_operands

    G.ia.pop_back();
    G.ia.push_back( static_cast<uint32_t>(G.ja.size()) );


#ifdef debug_sparsity_pattern
    cout<<"node number = "<<gref.Nodes()<<endl;
    cout<<"ia size = "<<G.ia.size()<<endl;
    cout<<"ja size = "<<G.ja.size()<<endl;
    cout<<"a size = "<<G.a.size()<<endl;

    //check whether all nodes in a manifold have the same row size
    set<NodeManifold<dim>*> manifolds;
    for (auto nit = gref.PerimeterNodesBegin(); nit != gref.NodesEnd(); nit++)
      if ( (*nit)->IsManifold() ) manifolds.insert((*nit)->Manifold());

    for (const auto & test_operand : test_operands) {
      csmp::Index prop_key = test_operand.first.key;
      for(auto md : manifolds) {
        auto branches = md->Branches();
        auto node1 = md->N(0);
        if(node1->Status(prop_key)==DIRICH) continue; //ignoring dirichlet nodes
        if(!gref.Contains(node1)) continue; //ignoring node outside of domain
        size_t idx1 = node1->Idx();
        size_t pos1 = DOF_indexes[idx1];
        if(pos1 != NULL_IDX) {
          uint32_t num1(0);
          for(auto in_dex = G.ia[pos1]; in_dex < G.ia[pos1+1]; in_dex++) num1++;
          for(auto n{1}; n < branches; n++) {
            auto node2 = md->N(n);
            if(node2->Status(prop_key)==DIRICH) continue; //ignoring dirichlet nodes
            if(!gref.Contains(node2)) continue; //ignoring node outside of domain
            auto idx2 = node2->Idx();
            auto pos2 = DOF_indexes[idx2];
            if(pos2 != NULL_IDX) {
              uint32_t num2(0);
              for(auto in_dex = G.ia[pos2]; in_dex < G.ia[pos2+1]; in_dex++) num2++;
              if(num1 != num2) {
                cout<<"manifold nodes row sizes do not match: "<<num1<<" vs "<<num2<<endl;
                cout<<"row of node 1 contains ";
                for(auto in_dex = G.ia[pos1]; in_dex < G.ia[pos1+1]; in_dex++) cout<<G.ja[in_dex]<<" ";
                cout<<endl;
                cout<<"row of node 2 contains ";
                for(auto in_dex = G.ia[pos2]; in_dex < G.ia[pos2+1]; in_dex++) cout<<G.ja[in_dex]<<" ";
                cout<<endl;

                //node 1
                set<uint32_t> node1_pos;
                cout<<"node 1:"<<endl;
                cout<<"  self = "<<pos1<<", ";
                node1_pos.insert(pos1);
                for (auto nb{0U}; nb < node1->Neighbors(); nb++) {
                  auto nd = node1->Neighbor(nb);
                  if(nd->Status(prop_key)==DIRICH) continue;
                  if(!gref.Contains(nd)) continue;
                  auto idx = nd->Idx();
                  auto pos = DOF_indexes[idx];
                  if (pos != NULL_IDX) {
                    cout << pos << " ";
                    node1_pos.insert(pos);
                  }
                }
                cout<<endl;
                for (auto mn{0U}; mn < branches; mn++) {
                  auto node = md->N(mn);
                  if (node->Status(prop_key) == DIRICH) continue;
                  if (node == node1) continue;
                  if(!gref.Contains(node)) continue;
                  //node itself
                  auto idx = node->Idx();
                  auto pos = DOF_indexes[idx];
                  if (pos != NULL_IDX) {
                    cout << "  manifold = " << pos << ", ";
                    node1_pos.insert(pos);
                  }
                  //neighbours
                  for (auto nn{0U}; nn < node->Neighbors(); nn++) {
                    auto nd = node->Neighbor(nn);
                    if (nd->Status(prop_key) == DIRICH) continue;
                    if(!gref.Contains(nd)) continue;
                    idx = nd->Idx();
                    pos = DOF_indexes[idx];
                    if (pos != NULL_IDX) {
                      cout << pos << " ";
                      node1_pos.insert(pos);
                    }
                  }
                  cout<<endl;
                }
                cout<<endl;

                //node 2
                set<uint32_t> node2_pos;
                cout<<"\nnode 2:"<<endl;
                cout<<"  self = "<<pos2<<", ";
                node2_pos.insert(pos2);
                for (auto nb{0U}; nb < node2->Neighbors(); nb++) {
                  auto nd = node2->Neighbor(nb);
                  if(nd->Status(prop_key)==DIRICH) continue;
                  if(!gref.Contains(nd)) continue;
                  auto idx = nd->Idx();
                  auto pos = DOF_indexes[idx];
                  if (pos != NULL_IDX) {
                    cout << pos << " ";
                    node2_pos.insert(pos);
                  }
                }
                cout<<endl;
                for (auto mn{0U}; mn < branches; mn++) {
                  auto node = md->N(mn);
                  if (node->Status(prop_key) == DIRICH) continue;
                  if (node == node2) continue;
                  if(!gref.Contains(node)) continue;
                  //node itself
                  auto idx = node->Idx();
                  auto pos = DOF_indexes[idx];
                  if (pos != NULL_IDX) {
                    cout << "  manifold = " << pos << ", ";
                    node2_pos.insert(pos);
                  }
                  //neighbours
                  for (auto nn{0U}; nn < node->Neighbors(); nn++) {
                    auto nd = node->Neighbor(nn);
                    if (nd->Status(prop_key) == DIRICH) continue;
                    if(!gref.Contains(nd)) continue;
                    idx = nd->Idx();
                    pos = DOF_indexes[idx];
                    if (pos != NULL_IDX) {
                      cout << pos << " ";
                      node2_pos.insert(pos);
                    }
                  }
                  cout<<endl;
                }

                cout<<endl;
                cout<<"\nnode1_pos size = "<<node1_pos.size()<<", contains:"<<endl;
                for(auto position : node1_pos) cout<<position<<" ";
                cout<<endl;
                cout<<"node2_pos size = "<<node2_pos.size()<<", contains:"<<endl;
                for(auto position : node2_pos) cout<<position<<" ";
                cout<<endl;

                throw csmp::Exception(ERROR, "generateSparsityPatternEliminatingEssentialConditions()",
                                      "manifold nodes row sizes do not match");
              }
            }
          }
        }
      }
    }
#endif

  }

  template void generateSparsityPatternEliminatingEssentialConditions<1U, Element> (CompressedRowMatrix&, const std::map<Parameter,size_t>&, std::vector<size_t>&, const  ModelSubDomain<1U,Element>&);
  template void generateSparsityPatternEliminatingEssentialConditions<2U, Element> (CompressedRowMatrix&, const std::map<Parameter,size_t>&, std::vector<size_t>&, const  ModelSubDomain<2U,Element>&);
  template void generateSparsityPatternEliminatingEssentialConditions<3U, Element> (CompressedRowMatrix&, const std::map<Parameter,size_t>&, std::vector<size_t>&, const  ModelSubDomain<3U,Element>&);
  template void generateSparsityPatternEliminatingEssentialConditions<1U, Face> (CompressedRowMatrix&, const std::map<Parameter,size_t>&, std::vector<size_t>&, const  ModelSubDomain<1U,Face>&);
  template void generateSparsityPatternEliminatingEssentialConditions<2U, Face> (CompressedRowMatrix&, const std::map<Parameter,size_t>&, std::vector<size_t>&, const  ModelSubDomain<2U,Face>&);
  template void generateSparsityPatternEliminatingEssentialConditions<3U, Face> (CompressedRowMatrix&, const std::map<Parameter,size_t>&, std::vector<size_t>&, const  ModelSubDomain<3U,Face>&);



CompressedRowMatrix::CompressedRowMatrix(std::vector<int32_t>&& input_ia,
                                         std::vector<int32_t>&& input_ja,
                                         std::vector<double>&& input_a)
  : ia(input_ia),
    ja(input_ja),
    a(input_a)
{
  if(Rows()>=numeric_limits<int32_t>::max()) throw out_of_range("CompressedRowMatrix(ctor): matrix size too large for SAMG solver");
  if(verbose_) cout<<"CompressedRowMatrix: called custom move constructor"<<endl;
}



CompressedRowMatrix::CompressedRowMatrix( csmp::SparseMatrix& spmat )
 {
    Initialize( spmat );
 }


/*
CompressedRowMatrix::CompressedRowMatrix( const CompressedRowMatrix& crm )
: ia(crm.ia),
  ja(crm.ja),
  a(crm.a)
 {
   if(verbose_) cout<<"CompressedRowMatrix: called copy constructor"<<endl;
 }

CompressedRowMatrix::CompressedRowMatrix( CompressedRowMatrix&& crm ) noexcept
  : ia(std::move(crm.ia)),
    ja(std::move(crm.ja)),
    a(std::move(crm.a))
  {
    if(verbose_) cout<<"CompressedRowMatrix: called move constructor"<<endl;
  }

CompressedRowMatrix&  CompressedRowMatrix::operator=( const CompressedRowMatrix& crm )
  {
       if ( &crm != this ) {
             ia = crm.ia;
             ja = crm.ja;
             a  = crm.a;
         }
       if(verbose_) cout<<"CompressedRowMatrix: called assignment operator"<<endl;
       return *this;
  }

  CompressedRowMatrix&  CompressedRowMatrix::operator=( CompressedRowMatrix&& crm ) noexcept
  {
    if ( &crm != this ) {
      ia = std::move(crm.ia);
      ja = std::move(crm.ja);
      a  = std::move(crm.a);
    }
    if(verbose_) cout<<"CompressedRowMatrix: called move assignment operator"<<endl;
    return *this;
  }
*/


/*
 Julian Mindel:  I proceeded to comment out the old code which contained the version of the () operator used before
 I have left it here below in the comment section for legacy purposes.

double  CompressedRowMatrix::operator()( uint32_t i, uint32_t j ) const
{
  // if the diagonal element is requested
  if ( i == j ) return  a[ static_cast<uint32_t>(ia[i]) ];

  // now all row elements are stored to the right of the diagonal (by convention)
  if (  a.size() - static_cast<uint32_t>(ia[i]) == ia.size() or  ia[i+1] - ia[i] == static_cast<int32_t>(ia.size()) )
    return  a[ static_cast<uint32_t>(ia[i]) + j - 1U ];

  // i is the diagonal element of the matrix
  for ( uint32_t  index=static_cast<uint32_t>(ia[i]); index <= ia.size(); index++ )
    if ( ja[index]  ==  static_cast<int32_t>(j)  ) return a[ index  ];

  return 0.;
}
*/


double CompressedRowMatrix::operator()( uint32_t i, uint32_t j ) const
  {
    assert( i < ja.size()-1U );
    assert( j < ja.size()-1U );
    assert( i < numeric_limits<int32_t>::max() );
    assert( j < numeric_limits<int32_t>::max() );

    return this->At(i,j);
  }


/**
 * range checked random access function that reports out of bound and missing element errors
 * @param i row index
 * @param j column index
 * @return double matrix element
 */
double CompressedRowMatrix::At( uint32_t i, uint32_t j ) const
{
  if ( i >= Rows() ) {
    cerr <<"\nCompressedRowMatrix::At("<< i <<","<< j <<"): ";
    cerr <<"Row access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::At");
  }
  if ( i >= Cols() ) {
    cerr <<"\nCompressedRowMatrix::At("<< i <<","<< j <<"): ";
    cerr <<"Column access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::At");
  }

  if(IsFormattedForSAMG()) {
    for (auto index = ia[i]; index < ia[i + 1]; index++) {
      if (ja[index - 1]-1 == j) return a[index - 1];
    }
  } else {
    if(j <= i) {
      for (auto index = ia[i]; index < ia[i + 1]; index++) {
        if (ja[index] == j) return a[index];
      }
    } else {
      for (auto index = (ia[i+1]-1); index >= ia[i]; index--) {
        if (ja[index] == j) return a[index];
      }
    }
  }

  return 0.;
}


void CompressedRowMatrix::Assign( uint32_t i, uint32_t j, double val )
{
  assert( i < ja.size()-1U );
  assert( j < ja.size()-1U );

  if ( i >= Rows() ) {
    cerr <<"\nCompressedRowMatrix::Assign("<< i <<","<< j <<","<< val <<"): ";
    cerr <<"Row access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::Assign");
  }
  if ( i >= Cols() ) {
    cout <<"\nCompressedRowMatrix::Assign("<< i <<","<< j <<","<< val <<"): ";
    cout <<"Column access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::Assign");
  }

  //Check if the compressed row matrix is converted into SAMG format.
  //This operation only applies to the matrix in its original form.
  if(IsFormattedForSAMG()) {
    cerr <<"\nCompressedRowMatrix::Assign: Error: this operation needs to be performed before the matrix is turned into SAMG format."<<endl;
    throw logic_error("CompressedRowMatrix::Assign: Error: this operation needs to be performed before the matrix is turned into SAMG format.");
  }

  if(j <= i) {
    for ( uint32_t  index = ia[i]; index < ia[i+1]; index++ )
      if (ja[index] == j) {a[index] = val; return;}
  } else {
    for ( uint32_t  index = ia[i+1]-1; index >= ia[i]; index-- )
      if (ja[index] == j) {a[index] = val; return;}
  }

  cerr <<"\nCompressedRowMatrix::Assign: Error: Cannot find target element in the compressed row matrix, i ="<<i<<", j = "<<j<<endl;
  cerr <<"candidate col IDs in the row are:"<<endl;
  for ( uint32_t index=ia[i]; index <ia[i+1]; index++ ) cerr<<ja[index]<<", ";
  cerr<<endl;
  throw runtime_error("CompressedRowMatrix::Assign: Error: Cannot find target element in the compressed row matrix.");

}


void CompressedRowMatrix::MultiplyEntryWith( uint32_t i, uint32_t j, double val )
{
  assert( i < ja.size()-1U );
  assert( j < ja.size()-1U );

  if ( i >= Rows() ) {
    cerr <<"\nCompressedRowMatrix::MultiplyEntryWith("<< i <<","<< j <<","<< val <<"): ";
    cerr <<"Row access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::MultiplyEntryWith");
  }
  if ( i >= Cols() ) {
    cout <<"\nCompressedRowMatrix::MultiplyEntryWith("<< i <<","<< j <<","<< val <<"): ";
    cout <<"Column access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::MultiplyEntryWith");
  }

  //Check if the compressed row matrix is converted into SAMG format.
  //This operation only applies to the matrix in its original form.
  if(IsFormattedForSAMG()) {
    cerr <<"\nCompressedRowMatrix::MultiplyEntryWith: Error: this operation needs to be performed before the matrix is turned into SAMG format."<<endl;
    throw runtime_error("CompressedRowMatrix::MultiplyEntryWith: Error: this operation needs to be performed before the matrix is turned into SAMG format.");
  }

  if(j <= i) {
    for (uint32_t index = ia[i]; index < ia[i + 1]; index++) {
      if (ja[index] == j) {a[index] *= val; return;}
    }
  } else {
    for (uint32_t index = ia[i+1]-1; index >= ia[i]; index--) {
      if (ja[index] == j) {a[index] *= val; return;}
    }
  }

  cerr <<"\nCompressedRowMatrix::MultiplyEntryWith: Error: Cannot find target element in the compressed row matrix."<<endl;
  throw runtime_error("CompressedRowMatrix::MultiplyEntryWith: Error: Cannot find target element in the compressed row matrix.");

}


void CompressedRowMatrix::Add( uint32_t i, uint32_t j, double val )
{
  // zero elements are not stored
  if ( !( val > 0. || val < 0. ) ) return;

  assert( i < ja.size()-1U );
  assert( j < ja.size()-1U );

  if ( i >= Rows() ) {
    cerr <<"\nCompressedRowMatrix::Add("<< i <<","<< j <<","<< val <<"): ";
    cerr <<"Row access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::Add");
  }
  if ( i >= Cols() ) {
    cout <<"\nCompressedRowMatrix::Add("<< i <<","<< j <<","<< val <<"): ";
    cout <<"Column access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::Add");
  }

  //Check if the compressed row matrix is converted into SAMG format.
  //This operation only applies to the matrix in its original form.
  if(IsFormattedForSAMG()) {
    cerr <<"\nCompressedRowMatrix::Add: Error: this operation needs to be performed before the matrix is turned into SAMG format."<<endl;
    throw runtime_error("CompressedRowMatrix::Add: Error: this operation needs to be performed before the matrix is turned into SAMG format.");
  }

  if(j <= i) {
    for (uint32_t index = ia[i]; index < ia[i + 1]; index++) {
      if (ja[index] == j) {a[index] += val; return;}
    }
  }else {
    for (uint32_t index = ia[i+1]-1; index >= ia[i]; index--) {
      if (ja[index] == j) {a[index] += val; return;}
    }
  }

  cerr <<"\nCompressedRowMatrix::Add: Error: Cannot find target element in the compressed row matrix, i ="<<i<<", j = "<<j<<endl;
  cerr <<"candidate col IDs in the row are:"<<endl;
  for ( uint32_t index=ia[i]; index <ia[i+1]; index++ ) cerr<<ja[index]<<", ";
  cerr<<endl;

  throw runtime_error("CompressedRowMatrix::Add: Error: Cannot find target element in the compressed row matrix.");

}

/**
Check whether the compressed row matrix has been converted into SAMG format.
*/
bool CompressedRowMatrix::IsFormattedForSAMG() const {
  auto it_end = ia.end()-1;
  if(*it_end > ja.size()) return true;

  return false;
}


void CompressedRowMatrix::ZeroRow( size_t row )
{
  assert( row < Rows() );
  assert( row >= 0);

  if ( row >= Rows() || row < 0) {
    cerr <<"\nCompressedRowMatrix::ZeroRow("<< row <<"): ";
    cerr <<"Row access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::ZeroRow");
  }

  //Check if the compressed row matrix is converted into SAMG format.
  //This operation should only apply to the matrix in its original form.
  if(IsFormattedForSAMG()) {
    cerr <<"\nCompressedRowMatrix::ZeroRow: Error: this operation should perform before the matrix is turned into SAMG format."<<endl;
    throw runtime_error("CompressedRowMatrix::ZeroRow:: Error: this operation should perform before the matrix is turned into SAMG format.");
  }

  for ( uint32_t index=ia[row]; index <ia[row+1]; index++ ) {
    a[index] = 0.;
  }

}


//this function only sets all values in vector a to zeros, but keep ia and ja unchanged.
void CompressedRowMatrix::Zero()
{
  std::fill(a.begin(), a.end(), 0.);
}


//this function erase all elements in ia, ja and a
void CompressedRowMatrix::Erase()
{
  ia.erase( ia.begin(), ia.end() );
  ia.erase( ja.begin(), ja.end() );
  a.erase( a.begin(), a.end() );
}


//this function add elements in one row by corresponding elements in another row, but excluding the diagonal elements
void CompressedRowMatrix::AddRowByAnotherRow(uint32_t i, uint32_t j)
{
    assert( i < ja.size()-1U );
    assert( j < ja.size()-1U );

    if ( i >= Rows() ) {
      cerr <<"\nCompressedRowMatrix::AddRowByAnotherRow("<< i <<","<< j <<"): ";
      cerr <<"Row access index i = "<< i << " out of range."<< endl;
      throw range_error("CompressedRowMatrix::AddRowByAnotherRow");
    }
    if ( i >= Rows() ) {
      cerr <<"\nCompressedRowMatrix::AddRowByAnotherRow("<< i <<","<< j <<"): ";
      cerr <<"Row access index i = "<< j << " out of range."<< endl;
      throw range_error("CompressedRowMatrix::AddRowByAnotherRow");
    }

    //Check if the compressed row matrix is converted into SAMG format.
    //This operation only applies to the matrix in its original form.
    if(IsFormattedForSAMG()) {
      cerr <<"\nCompressedRowMatrix::AddRowByAnotherRow: Error: this operation needs to be performed before the matrix is turned into SAMG format."<<endl;
      throw runtime_error("CompressedRowMatrix::AddRowByAnotherRow: Error: this operation needs to be performed before the matrix is turned into SAMG format.");
    }

    for (uint32_t index = ia[j]; index < ia[j + 1]; index++) {
      auto col_id = ja[index];
      auto value = a[index];
      if(col_id!=i && col_id!=j) //do not add diagonal elements
        Add(i, col_id, value);
    }

}


//this function assign elements in one row by corresponding elements in another row, but excluding the diagonal elements
void CompressedRowMatrix::AssignRowByAnotherRow(uint32_t i, uint32_t j)
{
    assert( i < ja.size()-1U );
    assert( j < ja.size()-1U );

    if ( i >= Rows() ) {
      cerr <<"\nCompressedRowMatrix::AssignRowByAnotherRow("<< i <<","<< j <<"): ";
      cerr <<"Row access index i = "<< i << " out of range."<< endl;
      throw range_error("CompressedRowMatrix::AssignRowByAnotherRow");
    }
    if ( i >= Rows() ) {
      cerr <<"\nCompressedRowMatrix::AssignRowByAnotherRow("<< i <<","<< j <<"): ";
      cerr <<"Row access index i = "<< j << " out of range."<< endl;
      throw range_error("CompressedRowMatrix::AssignRowByAnotherRow");
    }

    //Check if the compressed row matrix is converted into SAMG format.
    //This operation only applies to the matrix in its original form.
    if(IsFormattedForSAMG()) {
      cerr <<"\nCompressedRowMatrix::AssignRowByAnotherRow: Error: this operation needs to be performed before the matrix is turned into SAMG format."<<endl;
      throw runtime_error("CompressedRowMatrix::AssignRowByAnotherRow: Error: this operation needs to be performed before the matrix is turned into SAMG format.");
    }

    for (uint32_t index = ia[j]; index < ia[j + 1]; index++) {
      auto col_id = ja[index];
      auto value = a[index];
      if(col_id!=i && col_id!=j) //do not assign diagonal elements
        Assign(i, col_id, value);
    }
}




/**
 
Initialises the public CompressedRowMatrix vectors ia, ja, a 
from the given csmp::SparseMatrix.  

@param A The fully accumulated global solution matrix.

@section implementation Implementation

If nnu denotes the number of rows (variables), the non-zero entries 
of the i-th row (1 ≤ i ≤ nnu) are stored in a(j) where

ia(i) ≤ j ≤ ia(i+1)-1.

In particular, according to the above-mentioned convention about the location of the diagonal element,
a(ia(i)) contains the diagonal entry of row i. Note that ia(1) = 1 and  ia(nnu+1) = nna+1 where nna denotes
the total number of matrix entries stored.

ia, ja define the order of the off-diagonal elements in the solution
matrix storage vector a (that holds them row by row). 

The pointer vector ja has to be defined so that ja(j) (1 ≤ j ≤ nna) equals the original matrix’ column index of
a(j), ie, a(j) corresponds to the variable u(ja(j)). In particular, since a(ia(i)) contains the diagonal entry of row i,
we have ja(ia(i))=i.

While, in the compressed row format, all diagonal entries must be 
stored, off-diagonal entries are only stored if they are non-zero.  

@section application Application

Transfer of the global solution matrix to conventional solvers.  

@section messages Messages

Reports if the solution matrix contains zero diagonal entries.  
*/

void CompressedRowMatrix::Initialize( const SparseMatrix& A ) 
 {
      ia.resize( (A.Rows() + 1U) ); ia.shrink_to_fit();
      // ja is constructed with zero diagonal entries
      ja.resize( A.Entries(), 0 );  ja.shrink_to_fit();
      // 'a' stores the non-zero entries of the sparse matrix, row after row
      a.resize( ja.size() );        a.shrink_to_fit();

      // looping over all rows intializing ja and testing for diagonal entries which are zero
      // here n counts from 0 to j=nnu, i.e. all non-zero elements in the matrix
      uint32_t n(0U);
      ia[0] = 0;

      for ( auto i{0}; i < A.Rows(); i++ )
       {
          uint32_t  diag(UNSPECIFIED);
          // looping over the non-zero elements row i
          for ( auto rit=A.RowBegin(i); rit!=A.RowEnd(i); rit++ ) {
               // copying A's entry row(i) into the compressed row storage vector 'a'
               a[n]  = (*rit).second;
               // recording the corresponding column index in 'ja'
               // (NB: rit.first points to matrix column index from 0..rows-1)
               ja[n] = static_cast<uint32_t>((*rit).first);
               // if i=j, i.e., if this is a diagonal elemnt, its position is recorded by 'diag'
               // if the diagonal element is zero, however, it will not have been stored in 'ja'
               // so that this situation is never encountered and diag remains UNSPECIFIED
               if ( ja[n] == static_cast<uint32_t>(i) ) diag = static_cast<uint32_t>(n);
               n++;
	          }
          if ( diag == UNSPECIFIED ) {
               cout <<"\nCompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal: ";
               cout <<"\nSparseMatrix (rows=columns="<< A.Rows() <<") Zero entries (i=j): "<< endl;
               cout.setf(ios::scientific);
               long prec = cout.precision(15U);
               for ( auto i2{0}; i2 < A.Rows(); i2++ )
                 if ( std::fabs(A(i2,i2)) < std::numeric_limits<double>::epsilon() )
                   cout <<"\n\t"<< i2 <<": "<< A(i2,i2);
               cout << endl;
               cout.unsetf( ios::scientific );
               cout.precision(prec);
               A.Out();
               throw underflow_error("CompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal.");
            }
        
          // setting matrix such that diagonal element is at the beginning of next row 
          ia[i+1U]       = static_cast<uint32_t>(n);
          // inserting the diagonal elements at the beginning of each row
          const uint32_t istart = static_cast<uint32_t>(ia[i]);
          uint32_t jatemp = ja[ istart ];
          double  atemp  = a[ istart ];
          uint32_t dindex = diag;
          ja[ istart ]   = ja[ dindex ];
          a[ istart]     = a[ dindex ];
          a[ dindex ]    = atemp;
          ja[ dindex ]   = jatemp;
       }

     // converting C++ array indices (0..n-1) into Fortran indices (1..n) 
     //for ( vector<int32_t>::iterator it=ia.begin(); it!=ia.end(); it++ ) (*it)++;
     for (auto& it : ia) it++;
     //for ( vector<int32_t>::iterator it=ja.begin(); it!=ja.end(); it++ ) (*it)++;
     for (auto& it : ja) it++;

}  // end Initialize



//this function is similar to the previous one, but does not convert compressed row matrix to the SAMG format
void CompressedRowMatrix::ConvertFromSparseMatrix( const SparseMatrix& A ) 
 {
      ia.resize( (A.Rows() + 1U) ); ia.shrink_to_fit();
      // ja is constructed with zero diagonal entries
      ja.resize( A.Entries(), 0 );  ja.shrink_to_fit();
      // 'a' stores the non-zero entries of the sparse matrix, row after row
      a.resize( ja.size() );        a.shrink_to_fit();

      // looping over all rows intializing ja and testing for diagonal entries which are zero
      // here n counts from 0 to j=nnu, i.e. all non-zero elements in the matrix
      uint32_t n(0U);
      ia[0] = 0;

      for ( auto i{0}; i < A.Rows(); i++ )
       {
          uint32_t  diag(UNSPECIFIED);
          // looping over the non-zero elements row i
          for ( auto rit=A.RowBegin(i); rit!=A.RowEnd(i); rit++ ) {
               // copying A's entry row(i) into the compressed row storage vector 'a'
               a[n]  = (*rit).second;
               // recording the corresponding column index in 'ja'
               // (NB: rit.first points to matrix column index from 0..rows-1)
               ja[n] = static_cast<uint32_t>((*rit).first);
               // if i=j, i.e., if this is a diagonal elemnt, its position is recorded by 'diag'
               // if the diagonal element is zero, however, it will not have been stored in 'ja'
               // so that this situation is never encountered and diag remains UNSPECIFIED
               if ( ja[n] == static_cast<uint32_t>(i) ) diag = static_cast<uint32_t>(n);
               n++;
	          }
          if ( diag == UNSPECIFIED ) {
               cout <<"\nCompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal: ";
               cout <<"\nSparseMatrix (rows=columns="<< A.Rows() <<") Zero entries (i=j): "<< endl;
               cout.setf(ios::scientific);
               long prec = cout.precision(15U);
               for ( auto i2{0}; i2 < A.Rows(); i2++ )
                 if ( std::fabs(A(i2,i2)) < std::numeric_limits<double>::epsilon() )
                   cout <<"\n\t"<< i2 <<": "<< A(i2,i2);
               cout << endl;
               cout.unsetf( ios::scientific );
               cout.precision(prec);
               A.Out();
               throw underflow_error("CompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal.");
            }
            ia[i+1U]       = static_cast<uint32_t>(n);
       }

}  // end ConvertFromSparseMatrix


//this function converts compressed row matrix to the SAMG format
void CompressedRowMatrix::ConvertToSAMGFormat()
{
  if(IsFormattedForSAMG()) {
    cout<<"Already in SAMG format, nothing was done"<<endl;
    return;
  }

  for ( auto i{0}; i < ia.size()-1; i++ ) {
    uint32_t  diag(UNSPECIFIED);
    for(auto n=ia[i];n<ia[i+1];n++){
      if ( ja[n] == i ) {diag = n; break;}
    }
    if ( diag == UNSPECIFIED ) {
      cout <<"\nCompressedRowMatrix::ConvertToSAMGFormat: Error: no diagonal element can be found in row: "<<i<<endl;
      Out();
      throw runtime_error("CompressedRowMatrix::ConvertToSAMGFormat: Error: no diagonal element can be found.");
    }

    // setting matrix such that diagonal element is at the beginning of next row
    // inserting the diagonal elements at the beginning of each row
    const auto istart = ia[i];
    uint32_t jatemp = ja[ istart ];
    double  atemp  = a[ istart ];
    uint32_t dindex = diag;
    ja[ istart ]   = ja[ dindex ];
    a[ istart]     = a[ dindex ];
    a[ dindex ]    = atemp;
    ja[ dindex ]   = jatemp;
  }

  // converting C++ array indices (0..n-1) into Fortran indices (1..n)
  for (auto& it : ia) (it)++;
  for (auto& it : ja) (it)++;

}




/**
 
Initialises the public CompressedRowMatrix vectors ia, ja, a for given 
SparseMatrix in case the Point-based approach is selected.
*/
void CompressedRowMatrix::InitializePointBased( const SparseMatrix& A, size_t nsys ) 
 {
      // resize internal storage
      ia.resize( (A.Rows() + 1U) );
      ja.resize( A.Entries() );
      a.resize( ja.size() );
   
      map<size_t,double>::const_iterator rit;
      long      i, j, k, row;
      uint32_t   diag;
      bool      zero_diag_element(false);
	  
      std::vector<uint32_t>  temp( ja.size() ); // auxilary vector

      const size_t nnu_(A.Rows());
	    for ( k = 0U; k < nnu_; k++)
	      temp[k] = static_cast<uint32_t>(k%(nnu_/nsys)*nsys+k/(nnu_/nsys));
	  		
      for ( i=j=0U, ia[0]=0; i < nnu_; i++ )
       {
	        row = i%nsys*(nnu_/nsys)+i/nsys; // amending the order rows will be written in a[]
          for ( diag=-1, rit=A.RowBegin(row); rit!=A.RowEnd(row); rit++ )
            {
               a[j]  = (*rit).second;
               ja[j] = temp[(*rit).first];
               // rit.first points to matrix entries indexed from 0..rows-1
               if ( ja[j] == temp[row] ) diag = static_cast<uint32_t>(j);
               j++;
            }
          if ( diag == -1 ) zero_diag_element = true;
          
          ia[i+1] = static_cast<uint32_t>(j);
			
		      // inserting the diagonal elements at the beginning of each row
          const uint32_t istart = ia[i];
          ja[static_cast<uint32_t>(istart)] = ja[ static_cast<uint32_t>(diag) ];
          a[static_cast<uint32_t>(istart)]  = a[ static_cast<uint32_t>(diag) ];
          a[static_cast<uint32_t>(diag)]    = a[ static_cast<uint32_t>(istart) ];
          ja[static_cast<uint32_t>(diag)]   = ja[ static_cast<uint32_t>(istart) ];
       }

      if ( zero_diag_element ) {
          cout <<"\nCompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal: ";
          throw underflow_error("SparseM atrix::OutCompressedRowFormat");
       }
	   
     // converting C array indices (0..n-1) into Fortran indices (1..n) 
     for ( vector<int32_t>::iterator l=ia.begin(); l!=ia.end(); ++l )  (*l)++;
     for ( vector<int32_t>::iterator l=ja.begin(); l!=ja.end(); ++l )  (*l)++;

}  // end InitializePointBased



uint32_t CompressedRowMatrix::Rows() const
{
    return static_cast<uint32_t>(ia.size()-1U);
}


uint32_t CompressedRowMatrix::Cols() const
{
    return static_cast<uint32_t>(ia.size()-1U);
}


size_t CompressedRowMatrix::NonZeroEntries() const
{
    return ja.size();
}


/** Outputs matrix to screen.
*/
void CompressedRowMatrix::Out() const
 {
    bool SAMG_format = IsFormattedForSAMG();

    cout << flush <<"\nCompressedRowMatrix::Out: "<< endl;
    if(SAMG_format) cout<<"Matrix has been converted SAMG format"<<endl;
    cout <<"\nrow index vector 'ia' with size = "<<ia.size()<<"\n";
    for (auto it : ia) cout << it <<" ";
    cout <<"\ncolumn index vector 'ja' with size = "<<ja.size()<<"\n";
    for (auto it : ja) cout << it <<" ";
    cout <<"\nmatrix elements 'a' with size = "<<a.size()<<"\n";

    for(auto i{0U};i<ia.size()-1;i++) {
      for(size_t index=ia[i];index<ia[i+1];index++){
        if(!SAMG_format) cout<<ja[index]<<":"<<a[index]<<" ";
        else cout<<ja[index-1]<<":"<<a[index-1]<<" ";
      }
      cout<<endl;
    }
    cout << endl;
    cout.flush();   
 }



/** Outputs matrix to text file.
*/
void CompressedRowMatrix::Out( const string& outfile ) const
 {
    ofstream ofs(outfile);
    assert( ofs.is_open() );

    bool SAMG_format = IsFormattedForSAMG();

    ofs << flush <<"\nCompressedRowMatrix::Out: "<< endl;
    if(SAMG_format) cout<<"Matrix has been converted SAMG format"<<endl;
    ofs <<"\nrow index vector 'ia' with size = "<<ia.size()<<"\n";
    for (auto it : ia) ofs << it <<" ";
    ofs <<"\ncolumn index vector 'ja' with size = "<<ja.size()<<"\n";
    for (auto it : ja) ofs << it <<" ";
    ofs <<"\nmatrix elements 'a' with size = "<<a.size()<<"\n";

    const long precision = ofs.precision();
    ofs.precision(15);
    for(auto i{0U};i<ia.size()-1;i++) {
      for(size_t index=ia[i];index<ia[i+1];index++){
        if(!SAMG_format) ofs<<ja[index]<<":"<<a[index]<<" ";
        else ofs<<ja[index-1]<<":"<<a[index-1]<<" ";
      }
      ofs<<endl;
    }
    ofs.precision(precision);
    ofs << endl;
    ofs.flush();

 }


} // end csmp

