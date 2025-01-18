//
//  parentElementStorageOptions_Test.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 19/8/2024.
//

#include "parentElementStorageOptions_Test.h"
#include "Element.h"
#include "NodeParentElementVector.h"
#include "IsoparametricLinearQuadrilateral.h"

using namespace std;

namespace csmp {


void parentElementStorageOptions_Test::run()
 {
    Test_AlternativeContainersForSizeAndSpeed();
 }

struct NodeParentStorage {
    NodeParentStorage( size_t size ) : pe_ptrs_(size,nullptr), pn_idx_(size,ZERO) {}
    size_t Size() const { return pn_idx_.size(); }
    size_t SizeOf() const { return sizeof(*this) + pe_ptrs_.size() * sizeof(Element<3>*) + pn_idx_.size() * sizeof(ONE_BYTE_NUMBER); }
    void Resize( size_t size ) { pe_ptrs_.resize(size); pn_idx_.resize(size); }
    size_t Padding() const { return SizeOf() - (sizeof(pe_ptrs_) + pe_ptrs_.size() * sizeof(Element<3>*) + sizeof(pn_idx_) + pn_idx_.size() * sizeof(ONE_BYTE_NUMBER)); };
    // parent storage as in Node
    vector<Element<3>*>     pe_ptrs_;  ///< parent element pointers
    vector<ONE_BYTE_NUMBER> pn_idx_;   ///< local parent node number (0...nodes-1)
 };





/**
     Compares:
      - current (19/8/2024) implementation in Node using 2 vectors
      - map
      - unordered map
      - vector of pairs
*/
void parentElementStorageOptions_Test::Test_AlternativeContainersForSizeAndSpeed()
 {
    bool running_in_debug_mode = false;
#ifndef NDEBUG
    running_in_debug_mode = true;
#endif
   if ( running_in_debug_mode )
     cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: running in DEBUG mode. Switch on NDEBUG for meaningful comparisons"<< endl;
   
   // 0. build random integer generator to simulate varying parent element sizes
   // --------------------------------------------------------------------------
   // size generation engine
   // generating point locations with random number generator
   random_device rd; // obtain a random number from hardware
//   mt19937 gen(rd()); // seed the generator
   default_random_engine generator;
   // define the range of potential node parents
   const uint32_t min_vals{3}, max_vals{15};
   uniform_int_distribution<uint32_t> int_distr( min_vals, max_vals );
   
   
    // 1. comparisons for a million objects
    // ------------------------------------------------
    const long n_objects(1e6); // modify test here
    IsoparametricLinearQuadrilateral FE(2U);
    // generating a bunch of finite elements to have meaningful pointers
    vector<Element<3>> barebone_elmts( max_vals+1, Element<3>(&FE) );
    uint32_t counter{0u};
    for ( auto& it : barebone_elmts ) it.Idx( counter++ );
    uniform_int_distribution<uint32_t> nodes_distr( 0, FE.Nodes()-1 ); // node nums for quad
    
  cout <<"\n\n"<<"========================================================================================"<< endl << endl;


    // Node<>: current storage mechanism (139 Mb)
    // ==========================================
    {
      // 1.1 construction of container
      // -----------------------------
      vector<NodeParentStorage> storage1;
      storage1.reserve( n_objects );
      size_t estimated_size1{0ul};

      auto t0 = chrono::high_resolution_clock::now();
      for ( long i{0}; i<n_objects; ++i ) {
           storage1.emplace_back( NodeParentStorage( int_distr(generator) ) );
           estimated_size1 += storage1.back().SizeOf();
        }
      auto t1 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to build 'NodeParentStorage': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
      cout <<" milliseconds."<< endl;

      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: estimated size of 'NodeParentStorage': ";
      cout << static_cast<double>( estimated_size1 / 1e6 ) <<" Mb";
      cout << endl;
      
      // 1.2 resizing of container
      // -------------------------
      auto t2 = chrono::high_resolution_clock::now();
      // ___________________________________________________________________________
      for ( long i{0}; i<n_objects; ++i ) storage1[i].Resize( int_distr(generator) );
      //
      auto t3 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to resize 'NodeParentStorage': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t3-t2).count();
      cout <<" milliseconds."<< endl;
        
      // 1.3 iterating over container
      // ----------------------------
      auto t4 = chrono::high_resolution_clock::now();
      // ________________________________
      for ( long i{0}; i<n_objects; ++i ) {
        for ( uint32_t j{0u}; j<storage1[i].Size(); ++j ) {
             storage1[i].pn_idx_[j] = ONE;
             storage1[i].pe_ptrs_[j] = nullptr;
          }
        }
      //
      auto t5 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to iterate over 'NodeParentStorage': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t5-t4).count();
      cout <<" milliseconds."<< endl;
    }

  cout <<"\n\n"<<"========================================================================================"<< endl << endl;


    // NodeParentElementVector<>: proposed storage mechanism
    // =====================================================
    {
      // 1.1 construction of container
      // -----------------------------
      vector<NodeParentElementVector<3>> storage; // using default constructor
      storage.reserve( n_objects );
      size_t estimated_size{0ul};

      auto t0 = chrono::high_resolution_clock::now();
      for ( long i{0}; i<n_objects; ++i ) {
           storage.emplace_back( NodeParentElementVector<3>( int_distr(generator) ) );
           for ( uint32_t j{0u}; j<storage.back().Capacity(); ++j )
             storage.back().Assign( &barebone_elmts[ int_distr(generator) ], nodes_distr(generator) );
           // estimating the size of the storage
           estimated_size += storage.back().SizeOf();
        }
      auto t1 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to build 'NodeParentElementVector': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
      cout <<" milliseconds."<< endl;

      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: estimated size of 'NodeParentElementVector': ";
      cout << static_cast<double>( estimated_size / 1e6 ) <<" Mb";
      cout << endl;
      
      // 1.2 resizing of container
      // -------------------------
      // since there is no such operator it is just grown and shrink to fit is applied
      auto t2 = chrono::high_resolution_clock::now();
      // ___________________________________________________________________________
      for ( long i{0}; i<n_objects; ++i ) {
           storage[i].Assign( &barebone_elmts[ int_distr(generator) ], nodes_distr(generator) );
           storage[i].ShrinkToFit();
        }
      //
      auto t3 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to resize 'NodeParentElementVector': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t3-t2).count();
      cout <<" milliseconds."<< endl;
        
      // 1.3 iterating over container
      // ----------------------------
      Element<3>* eptr{ nullptr };
      size_t    n_node{0ul};
      
      auto t4 = chrono::high_resolution_clock::now();
      // ________________________________
      for ( long i{0}; i<n_objects; ++i ) {
        for ( uint32_t j{0u}; j<storage[i].Parents(); ++j ) {
             eptr    = storage[i].ParentElement(j);
             n_node += storage[i].LocalNodeNumber(j);
          }
        }
      auto t5 = chrono::high_resolution_clock::now();
      // avoiding that the optimiser removes code the results of which are not used
      cout <<"\n"<<"confirming use of values: "<< eptr <<" "<< n_node << endl;
      //
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to iterate over 'NodeParentElementVector': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t5-t4).count();
      cout <<" milliseconds."<< endl;

      // 1.4 testing individual functions
      // --------------------------------
      n_node  = 0ul;
      auto t6 = chrono::high_resolution_clock::now();
      // ________________________________
      for ( long i{0}; i<n_objects; ++i )
        for ( const auto& eit : barebone_elmts )
          n_node += storage[i].LocalNodeNumber( &eit );

      auto t7 = chrono::high_resolution_clock::now();
      cout <<"\n\t\t"<<"final value: "<< n_node << endl;
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time for 'NodeParentElementVector::LocalNodeNumber1': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t7-t6).count();
      cout <<" milliseconds."<< endl;
/*
      n_node  = 0ul;
      auto t8 = chrono::high_resolution_clock::now();
      // ________________________________
      for ( long i{0}; i<n_objects; ++i )
        for ( const auto& eit : barebone_elmts )
          n_node += storage[i].LocalNodeNumber2( &eit );

      auto t9 = chrono::high_resolution_clock::now();
      cout <<"\n\t\t"<<"final value: "<< n_node << endl;
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time for 'NodeParentElementVector::LocalNodeNumber2': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t9-t8).count();
      cout <<" milliseconds."<< endl;
*/
    }

  cout <<"\n\n"<<"========================================================================================"<< endl << endl;


    // map<uint32_t,Element<3>*>: standard map (430 Mb)
    // ================================================
    {
      // 1.1 construction of container
      // -----------------------------
      vector<map<uint32_t,Element<3>*>> storage2;
      storage2.reserve( n_objects );
      auto t0 = chrono::high_resolution_clock::now();
      for ( long i{0}; i<n_objects; ++i ) {
           map<uint32_t,Element<3>*> parent_map;
           const uint32_t desired_parents{ int_distr(generator) };
           for ( uint32_t j{0u}; j < desired_parents; ++j )
             parent_map.insert( make_pair(j,nullptr) );
           // storing
           storage2.emplace_back( parent_map );
        }
        
      auto t1 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to build 'map': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
      cout <<" milliseconds."<< endl;

      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: stop to read memory required by 'map': ";
      cout << endl;

      // 1.2 resizing of container
      // -------------------------
      auto t2 = chrono::high_resolution_clock::now();
      // ________________________________
      for ( long i{0}; i<n_objects; ++i ) {
           // additions or subtractions
           auto n_modifications = int_distr(generator);
           for ( uint32_t j{0u}; j<n_modifications; ++j )
             storage2[i].insert( make_pair( int_distr(generator), nullptr ) );
        }
      //
      auto t3 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to resize 'map': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t3-t2).count();
      cout <<" milliseconds."<< endl;
        
      // 1.3 iterating over container
      // ----------------------------
      auto t4 = chrono::high_resolution_clock::now();
      // ________________________________
      for ( long i{0}; i<n_objects; ++i ) {
        for ( auto& it : storage2[i] ) {
             it.second = static_cast<Element<3>*>(nullptr);
          }
        }
      //
      auto t5 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to iterate over 'map': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t5-t4).count();
      cout <<" milliseconds."<< endl;
    }

  cout <<"\n\n"<<"========================================================================================"<< endl << endl;


    // unordered_map<uint32_t,Element<3>*>: unordered (hash) map (350 Mb) NB: faster build than map
    // ============================================================================================
    {
      // 1.1 construction of container
      // -----------------------------
      vector<unordered_map<uint32_t,Element<3>*>> storage3;
      storage3.reserve( n_objects );
      auto t0 = chrono::high_resolution_clock::now();
      for ( long i{0}; i<n_objects; ++i ) {
           unordered_map<uint32_t,Element<3>*> parent_map;
           const uint32_t desired_parents{ int_distr(generator) };
           for ( uint32_t j{0u}; j < desired_parents; ++j )
             parent_map.insert( make_pair(j,nullptr) );
           // storing
           storage3.emplace_back( parent_map );
        }

      auto t1 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to build 'unordered_map': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
      cout <<" milliseconds."<< endl;

      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: stop to read memory required by 'unordered_map': ";
      cout << endl;

      // 1.2 resizing of container
      // -------------------------
      auto t2 = chrono::high_resolution_clock::now();
      // ________________________________
      for ( long i{0}; i<n_objects; ++i ) {
           // additions or subtractions
           auto n_modifications = int_distr(generator);
           for ( uint32_t j{0u}; j<n_modifications; ++j )
             storage3[i].insert( make_pair( int_distr(generator), nullptr ) );
        }
      //
      auto t3 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to resize 'unordered_map': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t3-t2).count();
      cout <<" milliseconds."<< endl;
        
      // 1.3 iterating over container
      // ----------------------------
      auto t4 = chrono::high_resolution_clock::now();
      // ________________________________
      for ( long i{0}; i<n_objects; ++i ) {
        for ( auto& it : storage3[i] ) {
             it.second = static_cast<Element<3>*>(nullptr);
          }
        }
      //
      auto t5 = chrono::high_resolution_clock::now();
      cout<<"\n"<<"Test_AlternativeContainersForSizeAndSpeed: time to iterate over 'unordered_map': ";
      cout << chrono::duration_cast<chrono::milliseconds>(t5-t4).count();
      cout <<" milliseconds."<< endl;
    }
    
    cout <<"\n"<<"finished 'Test_AlternativeContainersForSizeAndSpeed' test"<< endl;


 } // end Test_AlternativeContainersForSizeAndSpeed
  
/*
    short print_five{5}; // two-byte (-32,768 to 32,767)
    cout <<"\nsize of 'short': "<< sizeof(short) <<", value: "<< print_five << endl;
    cout << endl;

    or
    
    byte print_five{5}; // two-byte (-32,768 to 32,767)
    cout <<"\nsize of 'short': "<< sizeof(byte) <<", value: "<< (int)print_five << endl;
    cout << endl;

*/


} // end csmp
