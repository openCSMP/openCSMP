#include "FiniteElementManager1.h"
#include "Standard_IO_Handler.h"
#include "Exception.h"
#include "NameDemangler.h"

using namespace std;

namespace csmp {

// constructor providing options through choice of template arguments
template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::FiniteElementManager1()
  : line_elmt1(DIM), line_elmt2(DIM), line_elmt3(DIM),
    tria_elmt2(DIM), tria_elmt3(DIM),
    quad_elmt1(DIM), quad_elmt2(DIM), quad_elmt3(DIM)
    // other types do not require constructor arguments
  {
  }


template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
uint32_t  FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::NodesOfElementType( CSMP_FEM_TYPE etype )
 {
    return E( etype )->Nodes();
 }



template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
FiniteElement* const FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::LinearBarElement()
 {
//   static_assert( INTPOL_ORDER==1, "FiniteElementManager::LinearBarElement: Error: interpolation functions are nonlinear" );
   if constexpr ( !USING_LOCAL_COORDS ) return &line_elmt1;
   else return &line_elmt2;
   return nullptr;
 }
 
 
template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
FiniteElement* const FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::LinearTriangleElement()
 {
//   static_assert( INTPOL_ORDER==1, "FiniteElementManager::LinearTriangleElement: Error: interpolation functions are nonlinear" );
   if constexpr( !USING_LOCAL_COORDS ) {
        if constexpr( DIM != 3U ) return &tria_elmt1a;
        else return &tria_elmt1b;
     }
   else return &tria_elmt2;
   return nullptr;
 }


template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
FiniteElement* const FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::LinearTetrahedronElement()
 {
//   static_assert( INTPOL_ORDER==1, "FiniteElementManager::LinearTetrahedronElement: Error: interpolation functions are nonlinear" );
   if constexpr ( USING_LOCAL_COORDS == false )
     return &tet_elmt1;
   else
     return &tet_elmt2;
   return nullptr;
 }
 
 
 
 

/** supported types

linear isoparametric elements
ISOPARAMETRIC_LINEAR_TETRAHEDRON
ISOPARAMETRIC_LINEAR_HEXAHEDRON
ISOPARAMETRIC_LINEAR_PRISM
ISOPARAMETRIC_LINEAR_PYRAMID
ISOPARAMETRIC_LINEAR_TRIANGLE
ISOPARAMETRIC_LINEAR_QUADRILATERAL
ISOPARAMETRIC_LINEAR_BAR

quadratic isoparam elements
ISOPARAMETRIC_QUADRATIC_TETRAHEDRON
ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20
ISOPARAMETRIC_QUADRATIC_PYRAMID13
ISOPARAMETRIC_QUADRATIC_PRISM15
ISOPARAMETRIC_QUADRATIC_TRIANGLE
ISOPARAMETRIC_QUADRATIC_QUADRILATERAL
ISOPARAMETRIC_QUADRATIC_BAR

cubics (not implemented yet
ISOPARAMETRIC_CUBIC_TETRAHEDRON
ISOPARAMETRIC_CUBIC_HEXAHEDRON
ISOPARAMETRIC_CUBIC_PRISM
ISOPARAMETRIC_CUBIC_PYRAMID
ISOPARAMETRIC_CUBIC_TRIANGLE
ISOPARAMETRIC_CUBIC_QUADRILATERAL
ISOPARAMETRIC_CUBIC_BAR

barycentric quadratic forms
ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE
ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL
BARYCENTRIC_QUADRATIC_TETRAHEDRON
BARYCENTRIC_QUADRATIC_TRIANGLE
ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL

linear analytically integrated elements
LINEAR_TETRAHEDRON
LINEAR_CUBOID
LINEAR_TRIANGLE,
LINEAR_TRIANGLE3D
LINEAR_RECTANGLE
LINEAR_BAR

*/
template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
FiniteElement* const FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::E( int8_t etype )
 {
    // elements with linear interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_LINEAR_TETRAHEDRON: return &tet_elmt2;
               case ISOPARAMETRIC_LINEAR_HEXAHEDRON: return &hexa_elmt2;
               case ISOPARAMETRIC_LINEAR_PRISM: return &prism_elmt2;
               case ISOPARAMETRIC_LINEAR_PYRAMID: return &pyra_elmt2;
               case ISOPARAMETRIC_LINEAR_TRIANGLE: return &tria_elmt2;
               case ISOPARAMETRIC_LINEAR_QUADRILATERAL: return &quad_elmt2;
               case ISOPARAMETRIC_LINEAR_BAR: return &line_elmt2;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 2U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_LINEAR_TRIANGLE: return &tria_elmt2;
               case ISOPARAMETRIC_LINEAR_QUADRILATERAL: return &quad_elmt2;
               case ISOPARAMETRIC_LINEAR_BAR: return &line_elmt2;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 1U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_LINEAR_BAR: return &line_elmt2;
               default: return nullptr;
             }
          }
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 2 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON: return &tet_elmt3;
               case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27: return &hexa_elmt3;
               case ISOPARAMETRIC_QUADRATIC_PRISM18: return &prism_elmt3;
               case ISOPARAMETRIC_QUADRATIC_PYRAMID13: return &pyra_elmt3;
               case ISOPARAMETRIC_QUADRATIC_TRIANGLE: return &tria_elmt3;
               case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: return &quad_elmt3;
               case ISOPARAMETRIC_QUADRATIC_BAR: return &line_elmt3;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 2U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_QUADRATIC_TRIANGLE: return &tria_elmt3;
               case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: return &quad_elmt3;
               case ISOPARAMETRIC_QUADRATIC_BAR: return &line_elmt3;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 1U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_QUADRATIC_BAR: return &line_elmt3;
               default: return nullptr;
             }
          }
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS == false ) {
        if constexpr ( DIM == 3U ) {
           switch ( etype ) {
               case LINEAR_TETRAHEDRON: return &tet_elmt1;
               case LINEAR_CUBOID: return &hexa_elmt1;
               case LINEAR_TRIANGLE3D: return &tria_elmt1b;
               case LINEAR_QUADRILATERAL: return &quad_elmt1;
               case LINEAR_BAR: return &line_elmt1;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 2U ) {
           switch ( etype ) {
               case LINEAR_TRIANGLE: return &tria_elmt1a;
               case LINEAR_RECTANGLE: return &quad_elmt1;
               case LINEAR_BAR: return &line_elmt1;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 1U ) {
           switch ( etype ) {
               case LINEAR_BAR: return &line_elmt1;
               default: return nullptr;
             }
          }
      }

    cerr <<"\nFiniteElementManager1::Element(int8_t): Requested element is not available: ";
    cerr << etype <<" = "<< parseFiniteElementType(etype) <<" returning NULL pointer."<< std::endl;
      
    return nullptr;
    
 } // end E( int8_t )




template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
FiniteElement* const FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::E( CSMP_FEM_TYPE etype )
 {
    // elements with linear interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_LINEAR_TETRAHEDRON: return &tet_elmt2;
               case ISOPARAMETRIC_LINEAR_HEXAHEDRON: return &hexa_elmt2;
               case ISOPARAMETRIC_LINEAR_PRISM: return &prism_elmt2;
               case ISOPARAMETRIC_LINEAR_PYRAMID: return &pyra_elmt2;
               case ISOPARAMETRIC_LINEAR_TRIANGLE: return &tria_elmt2;
               case ISOPARAMETRIC_LINEAR_QUADRILATERAL: return &quad_elmt2;
               case ISOPARAMETRIC_LINEAR_BAR: return &line_elmt2;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 2U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_LINEAR_TRIANGLE: return &tria_elmt2;
               case ISOPARAMETRIC_LINEAR_QUADRILATERAL: return &quad_elmt2;
               case ISOPARAMETRIC_LINEAR_BAR: return &line_elmt2;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 1U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_LINEAR_BAR: return &line_elmt2;
               default: return nullptr;
             }
          }
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 2 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON: return &tet_elmt3;
               case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27: return &hexa_elmt3;
               case ISOPARAMETRIC_QUADRATIC_PRISM18: return &prism_elmt3;
               case ISOPARAMETRIC_QUADRATIC_PYRAMID13: return &pyra_elmt3;
               case ISOPARAMETRIC_QUADRATIC_TRIANGLE: return &tria_elmt3;
               case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: return &quad_elmt3;
               case ISOPARAMETRIC_QUADRATIC_BAR: return &line_elmt3;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 2U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_QUADRATIC_TRIANGLE: return &tria_elmt3;
               case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: return &quad_elmt3;
               case ISOPARAMETRIC_QUADRATIC_BAR: return &line_elmt3;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 1U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_QUADRATIC_BAR: return &line_elmt3;
               default: return nullptr;
             }
          }
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS == false ) {
        if constexpr ( DIM == 3U ) {
           switch ( etype ) {
               case LINEAR_TETRAHEDRON: return &tet_elmt1;
               case LINEAR_CUBOID: return &hexa_elmt1;
               case LINEAR_TRIANGLE3D: return &tria_elmt1b;
               case LINEAR_QUADRILATERAL: return &quad_elmt1;
               case LINEAR_BAR: return &line_elmt1;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 2U ) {
           switch ( etype ) {
               case LINEAR_TRIANGLE: return &tria_elmt1a;
               case LINEAR_RECTANGLE: return &quad_elmt1;
               case LINEAR_BAR: return &line_elmt1;
               default: return nullptr;
             }
          }
        if constexpr ( DIM == 1U ) {
           switch ( etype ) {
               case LINEAR_BAR: return &line_elmt1;
               default: return nullptr;
             }
          }
      }

    cerr <<"\nFiniteElementManager1::Element(CSMP_FEM_TYPE): Requested element is not available: ";
    cerr << etype <<" = "<< parseFiniteElementType(etype) <<" return NULL pointer."<< std::endl;
      
    return nullptr;
    
 } // end E( CSMP_FEM_TYPE )
 
 

/*
template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
uint32_t  FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::NodesOfElementType( CSMP_FEM_TYPE etype ) const
 {
    return E( etype )->Nodes();
 }
*/


template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
bool   FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::ContainsElementType( CSMP_FEM_TYPE etype ) const
 {
    // elements with linear interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_LINEAR_TETRAHEDRON: return true;
               case ISOPARAMETRIC_LINEAR_HEXAHEDRON: return true;
               case ISOPARAMETRIC_LINEAR_PRISM: return true;
               case ISOPARAMETRIC_LINEAR_PYRAMID: return true;
               case ISOPARAMETRIC_LINEAR_TRIANGLE: return true;
               case ISOPARAMETRIC_LINEAR_QUADRILATERAL: return true;
               case ISOPARAMETRIC_LINEAR_BAR: return true;
               default: return false;
             }
          }
        if constexpr ( DIM == 2U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_LINEAR_TRIANGLE: return true;
               case ISOPARAMETRIC_LINEAR_QUADRILATERAL: return true;
               case ISOPARAMETRIC_LINEAR_BAR: return true;
               default: return false;
             }
          }
        if constexpr ( DIM == 1U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_LINEAR_BAR: return true;
               default: return false;
             }
          }
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 2 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON: return true;
               case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27: return true;
               case ISOPARAMETRIC_QUADRATIC_PRISM18: return true;
               case ISOPARAMETRIC_QUADRATIC_PYRAMID13: return true;
               case ISOPARAMETRIC_QUADRATIC_TRIANGLE: return true;
               case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: return true;
               case ISOPARAMETRIC_QUADRATIC_BAR: return true;
               default: return false;
             }
          }
        if constexpr ( DIM == 2U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_QUADRATIC_TRIANGLE: return true;
               case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: return true;
               case ISOPARAMETRIC_QUADRATIC_BAR: return true;
               default: return false;
             }
          }
        if constexpr ( DIM == 1U ) {
           switch ( etype ) {
               case ISOPARAMETRIC_QUADRATIC_BAR: return true;
               default: return false;
             }
          }
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS == false ) {
        if constexpr ( DIM == 3U ) {
           switch ( etype ) {
               case LINEAR_TETRAHEDRON: return true;
               case LINEAR_CUBOID: return true;
               case LINEAR_TRIANGLE3D: return true;
               case LINEAR_QUADRILATERAL: return true;
               case LINEAR_BAR: return true;
               default: return false;
             }
          }
        if constexpr ( DIM == 2U ) {
           switch ( etype ) {
               case LINEAR_TRIANGLE: return true;
               case LINEAR_RECTANGLE: return true;
               case LINEAR_BAR: return true;
               default: return false;
             }
          }
        if constexpr ( DIM == 1U ) {
           switch ( etype ) {
               case LINEAR_BAR: return true;
               default: return false;
             }
          }
      }
      
    return false;
    
 } // end ContainsElementType(CSMP_FEM_TYPE)
 
 
 
 
template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
void  FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::CurrentElementTypes( list<CSMP_FEM_TYPE>& etypes ) const
 {
    etypes.erase( etypes.begin(), etypes.end() );
    
     // elements with linear interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
             etypes.push_back( ISOPARAMETRIC_LINEAR_TETRAHEDRON );
             etypes.push_back( ISOPARAMETRIC_LINEAR_HEXAHEDRON );
             etypes.push_back( ISOPARAMETRIC_LINEAR_PRISM );
             etypes.push_back( ISOPARAMETRIC_LINEAR_PYRAMID );
             etypes.push_back( ISOPARAMETRIC_LINEAR_TRIANGLE );
             etypes.push_back( ISOPARAMETRIC_LINEAR_QUADRILATERAL );
             etypes.push_back( ISOPARAMETRIC_LINEAR_BAR );
          }
        if constexpr ( DIM == 2U ) {
             etypes.push_back( ISOPARAMETRIC_LINEAR_TRIANGLE );
             etypes.push_back( ISOPARAMETRIC_LINEAR_QUADRILATERAL );
             etypes.push_back( ISOPARAMETRIC_LINEAR_BAR );
          }
        if constexpr ( DIM == 1U ) {
             etypes.push_back( ISOPARAMETRIC_LINEAR_BAR );
          }
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 2 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_TETRAHEDRON );
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 );
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_PRISM18 );
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_PYRAMID13 );
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_TRIANGLE );
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_QUADRILATERAL );
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_BAR );
          }
        if constexpr ( DIM == 2U ) {
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_TRIANGLE );
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_QUADRILATERAL );
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_BAR );
          }
        if constexpr ( DIM == 1U ) {
             etypes.push_back( ISOPARAMETRIC_QUADRATIC_BAR );
          }
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS == false ) {
        if constexpr ( DIM == 3U ) {
             etypes.push_back( LINEAR_TETRAHEDRON );
             etypes.push_back( LINEAR_CUBOID );
             etypes.push_back( LINEAR_TRIANGLE3D );
             etypes.push_back( LINEAR_QUADRILATERAL );
             etypes.push_back( LINEAR_BAR );
          }
        if constexpr ( DIM == 2U ) {
             etypes.push_back( LINEAR_TRIANGLE );
             etypes.push_back( LINEAR_RECTANGLE );
             etypes.push_back( LINEAR_BAR );
          }
        if constexpr ( DIM == 1U ) {
             etypes.push_back( LINEAR_BAR );
          }
      }
      
 } // end CurrentElementTypes
 
 
 
 
template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
void  FiniteElementManager1<DIM,INTPOL_ORDER,USING_LOCAL_COORDS>::Out() const
 {
     cout <<"\nFiniteElementManager1::Out: ";

     // elements with linear interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
             cout <<"\n\tcurrent volume elements: ";
             cout <<"\n\t\t"<< demangle(typeid(tet_elmt2).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(hexa_elmt2).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(prism_elmt2).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(pyra_elmt2).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(tria_elmt2).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(quad_elmt2).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(line_elmt2).name()) <<" ";
          }
        if constexpr ( DIM == 2U ) {
             cout <<"\n\tcurrent surface elements: ";
             cout <<"\n\t\t"<< demangle(typeid(tria_elmt2).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(quad_elmt2).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(line_elmt2).name()) <<" ";
          }
        if constexpr ( DIM == 1U ) {
            cout <<"\n\tcurrent line elements: ";
            cout <<"\n\t\t"<< demangle(typeid(line_elmt2).name()) <<" ";
          }
        cout << endl << endl;
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 2 && USING_LOCAL_COORDS ) {
        if constexpr ( DIM == 3U ) {
             cout <<"\n\tcurrent volume elements: ";
             cout <<"\n\t\t"<< demangle(typeid(tet_elmt3).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(hexa_elmt3).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(prism_elmt3).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(pyra_elmt3).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(tria_elmt3).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(quad_elmt3).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(line_elmt3).name()) <<" ";
          }
        if constexpr ( DIM == 2U ) {
             cout <<"\n\tcurrent surface elements: ";
             cout <<"\n\t\t"<< demangle(typeid(tria_elmt3).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(quad_elmt3).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(line_elmt3).name()) <<" ";
          }
        if constexpr ( DIM == 1U ) {
            cout <<"\n\tcurrent line elements: ";
            cout <<"\n\t\t"<< demangle(typeid(line_elmt3).name()) <<" ";
          }
        cout << endl << endl;
      }

    // elements with quadratic interpolation and shape functions
    if constexpr ( INTPOL_ORDER == 1 && USING_LOCAL_COORDS == false ) {
        if constexpr ( DIM == 3U ) {
             cout <<"\n\tcurrent volume elements: ";
             cout <<"\n\t\t"<< demangle(typeid(tet_elmt1).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(hexa_elmt1).name()) <<" ";
             // cout <<"\n\t\t"<< demangle(typeid(prism_elmt1).name()) <<" ";
             // cout <<"\n\t\t"<< demangle(typeid(pyra_elmt1).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(tria_elmt1b).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(quad_elmt1).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(line_elmt1).name()) <<" ";
          }
        if constexpr ( DIM == 2U ) {
             cout <<"\n\tcurrent surface elements: ";
             cout <<"\n\t\t"<< demangle(typeid(tria_elmt1a).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(quad_elmt1).name()) <<" ";
             cout <<"\n\t\t"<< demangle(typeid(line_elmt1).name()) <<" ";
          }
        if constexpr ( DIM == 1U ) {
            cout <<"\n\tcurrent line elements: ";
            cout <<"\n\t\t"<< demangle(typeid(line_elmt1).name()) <<" ";
          }
        cout << endl << endl;
      }
     
 } // end Out


//template<uint32_t DIM, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
template class FiniteElementManager1<1U,1,true>;
template class FiniteElementManager1<2U,1,true>;
template class FiniteElementManager1<3U,1,true>;

template class FiniteElementManager1<1U,2,true>;
template class FiniteElementManager1<2U,2,true>;
template class FiniteElementManager1<3U,2,true>;

template class FiniteElementManager1<1U,1,false>;
template class FiniteElementManager1<2U,1,false>;
template class FiniteElementManager1<3U,1,false>;

template class FiniteElementManager1<1U,2,false>;
template class FiniteElementManager1<2U,2,false>;
template class FiniteElementManager1<3U,2,false>;
 
} // end namespace csmp
 
 
 
 
 




