#include "FiniteElementManager.h"
#include "Standard_IO_Handler.h"
#include "Exception.h"

#include "LinearCuboid.h"
#include "LinearRectangle.h"

#include "LinearLineElement.h"
#include "IsoparametricQuadraticLineElement.h"
#include "LinearTriangle.h"
#include "LinearTriangle3D.h"
#include "LinearTetrahedron.h"
#include "IsoparametricQuadraticTriangle.h"
#include "IsoparametricQuadraticTetrahedron.h"

#include "IsoparametricLinearLineElement.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearPrism.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearQuadrilateral.h"

#include "IsoparametricQuadraticHexahedron.h"
#include "IsoparametricQuadraticPyramid.h"
#include "IsoparametricQuadraticPrism.h"
#include "IsoparametricQuadraticQuadrilateral.h"


using namespace std;

namespace csmp {

FiniteElementManager::FiniteElementManager()
  : dimensions(2), 
    interpolation(1),
    mixed_element_formulation(false),
    hexa_ptr(0), pyra_ptr(0), pris_ptr(0), tetr_ptr(0),
    quad_ptr(new IsoparametricLinearQuadrilateral()), 
    tria_ptr(new LinearTriangle()),
    line_ptr(new LinearLineElement(2))
 {
 }
 

FiniteElementManager::FiniteElementManager( size_t dim, 
                                            size_t interpolation_order,
                                            bool isoparametric )
  : dimensions(dim), 
    interpolation(interpolation_order),
    mixed_element_formulation(false),
    hexa_ptr(0), pyra_ptr(0), pris_ptr(0), tetr_ptr(0),
    quad_ptr(0), tria_ptr(0),
    line_ptr(0)
 {
    InitializeElements( dim, interpolation_order, isoparametric );
      
 } // end FiniteElementManager(custom constructor)




FiniteElementManager::FiniteElementManager( const FiniteElementManager& mgr )
  :  hexa_ptr(0), pyra_ptr(0), pris_ptr(0), tetr_ptr(0),
     quad_ptr(0), tria_ptr(0),
     line_ptr(0)
 {
    *this = mgr; 
 }
 


FiniteElementManager& FiniteElementManager::operator=( const FiniteElementManager& mgr )
 {
    if ( &mgr != this ) {
         dimensions                = mgr.dimensions;
         interpolation             = mgr.interpolation;
         mixed_element_formulation = mgr.mixed_element_formulation;
         
         // this assumes that the mixed element formulation is isoparametric
         InitializeElements( dimensions, interpolation, mixed_element_formulation );
      }
    return *this; 
 }

 
 

FiniteElementManager::~FiniteElementManager()
 {
    // volume elements
    delete hexa_ptr;
    delete pyra_ptr;
    delete pris_ptr;
    delete tetr_ptr;
    // surface elements
    delete quad_ptr;
    delete tria_ptr;
    // line elements
    delete line_ptr;
 }



FiniteElement*  FiniteElementManager::LinearBarElement() const
 {
    if ( line_ptr != NULL ) return line_ptr;
    std::cerr <<"\nFiniteElementManager::LinearBarElement: Not available."<< std::endl;
    return NULL;
 }
 
 
FiniteElement*  FiniteElementManager::LinearTriangleElement() const
 {
    if ( tria_ptr != NULL ) return tria_ptr;
    std::cerr <<"\nFiniteElementManager::LinearTriangleElement: Not available."<< std::endl;
    return NULL;
 }
 
 
FiniteElement*  FiniteElementManager::LinearTetrahedronElement() const
 {
    if ( tetr_ptr != NULL ) return tetr_ptr;
    std::cerr <<"\nFiniteElementManager::LinearTetrahedronElement: Not available."<< std::endl;
    return NULL;
 }


FiniteElement*  FiniteElementManager::E( int32 e_type ) const
 {
    if ( dimensions == 3 ) {
	     if ( hexa_ptr != NULL && e_type == hexa_ptr->ElementType() ) return hexa_ptr;
	     if ( pyra_ptr != NULL && e_type == pyra_ptr->ElementType() ) return pyra_ptr;
	     if ( pris_ptr != NULL && e_type == pris_ptr->ElementType() ) return pris_ptr;
	     if ( tetr_ptr != NULL && e_type == tetr_ptr->ElementType() ) return tetr_ptr;
      }
    if ( quad_ptr != NULL && e_type == quad_ptr->ElementType() ) return quad_ptr;
    if ( tria_ptr != NULL && e_type == tria_ptr->ElementType() ) return tria_ptr;
    if ( line_ptr != NULL && e_type == line_ptr->ElementType() ) return line_ptr;

    std::cerr <<"\nFiniteElementManager::Element: Requested element is not available: ";
    std::cerr << e_type <<" = "<< parseFiniteElementType(e_type) <<" returning NULL pointer."<< std::endl;
      
    return NULL;
 }


FiniteElement*  FiniteElementManager::E( CSMP_FEM_TYPE e_type ) const
 {
    if ( dimensions == 3 ) {
	     if ( hexa_ptr != NULL && e_type == hexa_ptr->ElementType() ) return hexa_ptr;
	     if ( pyra_ptr != NULL && e_type == pyra_ptr->ElementType() ) return pyra_ptr;
	     if ( pris_ptr != NULL && e_type == pris_ptr->ElementType() ) return pris_ptr;
	     if ( tetr_ptr != NULL && e_type == tetr_ptr->ElementType() ) return tetr_ptr;
      }
    if ( quad_ptr != NULL && e_type == quad_ptr->ElementType() ) return quad_ptr;
    if ( tria_ptr != NULL && e_type == tria_ptr->ElementType() ) return tria_ptr;
    if ( line_ptr != NULL && e_type == line_ptr->ElementType() ) return line_ptr;

    std::cerr <<"\nFiniteElementManager::Element: Requested element is not available: ";
    std::cerr << e_type <<" = "<< parseFiniteElementType(e_type) <<" return NULL pointer."<< std::endl;
      
    return NULL;
 }


size_t  FiniteElementManager::NodesOfElementType( CSMP_FEM_TYPE etype ) const
 {
    return E( etype )->Nodes();
 }



size_t  FiniteElementManager::Dimensions() const
 { return dimensions; }







void FiniteElementManager::InitializeElements( size_t dim, size_t interpolation_order, bool isoparametric )
 {
    dimensions = dim;
    mixed_element_formulation = false;
    
    // volume elements
    delete hexa_ptr;  hexa_ptr=0;
    delete pyra_ptr;  pyra_ptr=0;
    delete pris_ptr;  pris_ptr=0;
    delete tetr_ptr;  tetr_ptr=0;
    // surface elements
    delete quad_ptr;  quad_ptr=0;
    delete tria_ptr;  tria_ptr=0;
    // line elements
    delete line_ptr;  line_ptr=0;

    // 1D Models
    if ( dim == 1U ) {
         if ( isoparametric ) {
              if      ( interpolation_order == 1U ) line_ptr = new IsoparametricLinearLineElement(1);
              else if ( interpolation_order == 2U ) line_ptr = new IsoparametricQuadraticLineElement(1);
	          else
	          throw csmp::Exception( FATAL_ERROR, "FiniteElementManager(constructor):",
	                         "Desired order of interpolation functions is not available");
           }
         else {
              if      ( interpolation_order == 1 ) line_ptr = new LinearLineElement(1);
              else throw csmp::Exception( FATAL_ERROR, "FiniteElementManager(constructor):",
                                                "1D analytically integrated >=quadratic elements are not available");
           }
      }
                         
    // 2D Models                     
    else if ( dim == 2U ) {
         if ( interpolation_order == 1U ) {
              if ( isoparametric ) {
                   quad_ptr = new IsoparametricLinearQuadrilateral(2); 
                   tria_ptr = new IsoparametricLinearTriangle();
                   line_ptr = new IsoparametricLinearLineElement(2);
                }
              else {
                   // quad_ptr = new LinearQuadrilateral(2); 
                   tria_ptr = new LinearTriangle();
                   line_ptr = new LinearLineElement(2);
                }
           }
         else if ( interpolation_order == 2U ) {
              if ( isoparametric ) {
                   quad_ptr = new IsoparametricQuadraticQuadrilateral(2); 
                   tria_ptr = new IsoparametricQuadraticTriangle(2);
                   line_ptr = new IsoparametricQuadraticLineElement(2);
                }
              else {
                   throw csmp::Exception( WARNING, "FiniteElementManager::InitializeElements:",
                                  "In this version only isoparametric elements are available for quadratic interpolation");
                }
           }
         else
         throw csmp::Exception( FATAL_ERROR, "FiniteElementManager::InitializeElements:",
                         "Desired order of interpolation functions is not available in this version");
      }
      
    // 3D Models
    else if ( dim == 3U ) {
         if ( interpolation_order == 1U ) {
              if ( isoparametric ) {
                   hexa_ptr = new IsoparametricLinearHexahedron();
                   pyra_ptr = new IsoparametricLinearPyramid();
                   pris_ptr = new IsoparametricLinearPrism();
                   tetr_ptr = new IsoparametricLinearTetrahedron(); 
                   quad_ptr = new IsoparametricLinearQuadrilateral(3); 
                   tria_ptr = new IsoparametricLinearTriangle(3);
                   line_ptr = new IsoparametricLinearLineElement(3);
                }
              else {
				   hexa_ptr = new LinearCuboid();
				   quad_ptr = new LinearRectangle(3);
                   tetr_ptr = new LinearTetrahedron(); 
                   tria_ptr = new LinearTriangle3D();
                   line_ptr = new LinearLineElement(3);
                }
           }
         else if ( interpolation_order == 2U ) {
              if ( isoparametric ) {
                   hexa_ptr = new IsoparametricQuadraticHexahedron();
                   pyra_ptr = new IsoparametricQuadraticPyramid();
                   pris_ptr = new IsoparametricQuadraticPrism();
                   tetr_ptr = new IsoparametricQuadraticTetrahedron();
                   quad_ptr = new IsoparametricQuadraticQuadrilateral(3); 
                   tria_ptr = new IsoparametricQuadraticTriangle(3);
                   line_ptr = new IsoparametricQuadraticLineElement(3);
                }
              else {
                   throw csmp::Exception( WARNING, "FiniteElementManager::InitializeElements:",
                                  "In this version only isoparametric elements are available for quadratic interpolation");
                }
           }
         else
         throw csmp::Exception( FATAL_ERROR, "FiniteElementManager::InitializeElements:",
                         "Desired order of interpolation functions is not available");
      }
    else {
         throw csmp::Exception( FATAL_ERROR, "FiniteElementManager::InitializeElements:",
                         "Error in constructor argument 1");
      }

 } // end Initialize





    // can be changed, but intermediate nodes will become defunct
void      FiniteElementManager::InterpolationOrder( size_t interpolation_order )
 {
    interpolation = interpolation_order;
 }
 
 
size_t FiniteElementManager::InterpolationOrder() const
 {
    return interpolation;
 }
 
 



bool   FiniteElementManager::ContainsElementType( CSMP_FEM_TYPE e_type ) const
 {
    if ( dimensions == 3U ) {
	    if ( hexa_ptr != NULL && e_type == hexa_ptr->ElementType() ) return true;
	    if ( pyra_ptr != NULL && e_type == pyra_ptr->ElementType() ) return true;
	    if ( pris_ptr != NULL && e_type == pris_ptr->ElementType() ) return true;
	    if ( tetr_ptr != NULL && e_type == tetr_ptr->ElementType() ) return true;
      }
    if ( quad_ptr != NULL && e_type == quad_ptr->ElementType() ) return true;
	  if ( tria_ptr != NULL && e_type == tria_ptr->ElementType() ) return true;
	  if ( line_ptr != NULL && e_type == line_ptr->ElementType() ) return true;
      
    return false;   
 }
 
 

void  FiniteElementManager::CurrentElementTypes( std::list<CSMP_FEM_TYPE>& etypes ) const
 {
    etypes.erase( etypes.begin(), etypes.end() );
    if ( dimensions == 3U ) {
	    etypes.push_back( hexa_ptr->ElementType() );
	    etypes.push_back( pyra_ptr->ElementType() );
	    etypes.push_back( pris_ptr->ElementType() );
	    etypes.push_back( tetr_ptr->ElementType() );
      }
    if ( dimensions == 2U || dimensions == 3U ) {
	    etypes.push_back( quad_ptr->ElementType() );
	    etypes.push_back( tria_ptr->ElementType() );
      }
    etypes.push_back( line_ptr->ElementType() );
 }
 
 

void  FiniteElementManager::Out() const
 {
     cout <<"\nFiniteElementManager::Out: ";
     
     /// @todo (2-P) Remove typid by name fct
     if ( hexa_ptr or pyra_ptr or pris_ptr or tetr_ptr ) 
        cout <<"\n\tcurrent volume elements: ";
     if ( hexa_ptr  ) cout << typeid(*hexa_ptr).name() <<" ";
     if ( pyra_ptr ) cout << typeid(*pyra_ptr).name() <<" ";
     if ( pris_ptr ) cout << typeid(*pris_ptr).name() <<" ";
     if ( tetr_ptr ) cout << typeid(*tetr_ptr).name() <<" ";
     
     if ( quad_ptr or tria_ptr ) 
       cout <<"\n\tcurrent surface elements: ";
     if ( quad_ptr  ) cout << typeid(*quad_ptr).name() <<" ";
     if ( tria_ptr ) cout << typeid(*tria_ptr).name() <<" ";

     if ( line_ptr ) {
          cout <<"\n\tcurrent line elements: ";
          cout << typeid(*line_ptr).name() <<" ";
       }
     cout << endl;
 }



 
} // end namespace csmp
 
 
 
 
 




