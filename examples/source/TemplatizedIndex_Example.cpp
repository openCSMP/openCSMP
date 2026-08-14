#include "TemplatizedIndex_Example.h"
#include "Index.h"

#include <iostream>
#include <set>
#include <vector>
#include <functional>
#include <numeric>

using namespace std;

namespace csmp {

void TemplatizedIndex_Example::Specifications()
{
  SetTitle( "How to templatize Index class to permit static polymorphic handling of variables" );
  SetDifficulty( 1 );
  SetCategory( "Software Functionality" );
  AddAuthor( "SKM" );
  AddDescription( "source in: TemplatedIndex_Example.cpp" );
  AddDescription( "logical extension of Index for variables that are known already at compile time" );
  AddRequirement( "Index.h");
  
} // Specifications


// type build from enums, index can still be assigned at runtime
/*
template<VARIABLE_TYPE t,PLACEMENT p>
struct INDEX : public Index {
     // use singleton to count created variables (like compile time map<name,index>)
     explicit INDEX( unsigned int i ) : index(i) {};
     INDEX( const Index& i ) : index(i.index) {};
     Index& operator=( const Index& i ) { index = i.index; return *this; };
     ~INDEX() {};
     enum { type=t };
     enum { place=p };
     unsigned int  index;
 };
*/ 



 
//template<typename A, typename B>
//void add( const A&, const B& );
 
static void add( const Index&, const Index& ) {
    std::cout <<"\nadd: standard version, no specialisation found."<< std::endl;
 } 

static void add( const INDEX<SCALAR,NODE>&, const Index& ) {
    std::cout <<"\nadd: standard version, SCALAR/NODE partial version."<< std::endl;
 } 

static void add( const INDEX<SCALAR,NODE>&, const INDEX<SCALAR,NODE>& ) {
    std::cout <<"\nadd: standard version, SCALAR/NODE full version."<< std::endl;
 } 
 
 
// PROCESSING 

/// generic variable processor using STL functors for math operations
template<template<typename> class operation> 
void process( const Index& A, const Index& B ) {
    std::cout <<"\nprocess: standard version, no specialisation found."<< std::endl;
 } 
 
template<>
void process<std::plus>( const Index&, const Index& ) {
    std::cout <<"\nprocess: plus: standard version, no specialisation found."<< std::endl;
    double res = std::plus<double>()(3.,5.);
    std::cout <<" "<< res <<" ";
 } 
 
 
/**

Necessary combinations for all math operations, (numbperm(n,r)=72 - non-physical ones

node, cpoint, face, element, region,   scalar, vector, tensor

on the nodes

scalar,scalar,node  +,-,*,/, sqrt, log, exp etc.
vector,vector,node
tensor,tensor,node

*/ 
template<VARIABLE_TYPE t,PLACEMENT p> void process( const INDEX<t,p>& );

template<> 
void process( const INDEX<SCALAR,ELEMENT_INTEGRATION_POINT>& );
void process( const INDEX<SCALAR,FACE>& );
void process( const INDEX<SCALAR,INTER_FACE>& );
void process( const INDEX<SCALAR,REGION>& );

void process( const INDEX<VECTOR,NODE>& );
void process( const INDEX<VECTOR,ELEMENT_INTEGRATION_POINT>& );
void process( const INDEX<VECTOR,FACE>& );
void process( const INDEX<SCALAR,INTER_FACE>& );
void process( const INDEX<VECTOR,REGION>& );

void process( const INDEX<TENSOR,NODE>& );
void process( const INDEX<TENSOR,ELEMENT_INTEGRATION_POINT>& );
void process( const INDEX<TENSOR,FACE>& );
void process( const INDEX<SCALAR,INTER_FACE>& );
void process( const INDEX<TENSOR,REGION>& );


// illustration of Index - INDEX combinations
/*
void process( Region<3>& gref, const INDEX<SCALAR,NODE>&, const ScalarVariable& sc )
 {
    fill( gref.NodesBegin(), gref.NodesEnd(), sc );
    // minus, plus, etc.
    transform( gref.begin(), gref.end(), );
 }
 
 // as called
process( gref, permeability, sc );
*/


void  TemplatizedIndex_Example::Run()
 {
    Index idx;
    INDEX<SCALAR,NODE>  scalar(0);
    INDEX<VECTOR,NODE>  vector(1);
    
    add( scalar, scalar );
    add( vector, vector );
    add( scalar, vector );
    add( scalar, idx );
    
    process<plus>( idx, idx );
    
//    process( scalar );
 //   process( vector );
            
    set<uint32_t>  test_set;
    
    test_set.insert( 1U );
    test_set.insert( 2U );
    test_set.insert( 6U );
    test_set.insert( 3U );
    test_set.insert( 4U );
    
    for ( uint32_t i=20U; i>0; i-- )
      test_set.insert( i );
    
    for ( set<uint32_t>::const_reverse_iterator
          it=test_set.rbegin(); it!=test_set.rend(); it++ )
      cout << *it <<" ";

    cout <<"\nmain: last element in set: "<< (*test_set.rbegin()) << endl;
    
} // end Run

} // csmp

