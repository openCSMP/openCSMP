#include <sstream>
#include "AP_testUtilities.h"
//#include "FiniteVolumeTraits.h"
#include "FiniteVolumePolicy.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "NodeManifold.h"
#include "FiniteElementPolicy.h"

using namespace std;

namespace csmp {

typedef FiniteVolumePolicy<3U,Element> EFT3;

void rhinoOutput( const Model<3U>& sgroup )
 {
    system( "del iso_lin_*.txt" );
    
    const Region<3U>  sg(sgroup.Region("Model")); 
  
    // reading elements
    for ( auto eit=sg.CellsBegin(); eit!=sg.CellsEnd(); eit++ )
    {
      const CSMP_FEM_TYPE elType((*eit)->FE_Type());
      
      //print isoparametric linear triangle
      if(elType == ISOPARAMETRIC_LINEAR_TRIANGLE)
        printElementStencil( *(*eit), "iso_lin_tri.txt");
      
      //print isoparametric linear quadrilateral
      else if(elType == ISOPARAMETRIC_LINEAR_QUADRILATERAL)
        printElementStencil( *(*eit), "iso_lin_quad.txt");
      
      //print isoparametric linear tetrahedron
      else if(elType == ISOPARAMETRIC_LINEAR_TETRAHEDRON)
        printElementStencil( *(*eit), "iso_lin_tet.txt");
      
      //print isoparametric linear hexahedron
      else if(elType == ISOPARAMETRIC_LINEAR_HEXAHEDRON)
        printElementStencil( *(*eit), "iso_lin_hex.txt");
      
      //print isoparametric linear prism
      else if(elType == ISOPARAMETRIC_LINEAR_PRISM)
        printElementStencil( *(*eit), "iso_lin_pri.txt");
      
      //print isoparametric linear pyramid
      else if(elType == ISOPARAMETRIC_LINEAR_PYRAMID)
        printElementStencil( *(*eit), "iso_lin_pyr.txt");
          
    }
 } // end rhinoOutput
 


static void writeSurfaceFacet(std::stringstream & ss, const EFT3& efvt, const Element<3U>& e, uint32_t iFacet)
 {
   ss << "\nSrfPt ";
     
   for(auto iPoint = 0; iPoint < e.FV()->FacetPoints(iFacet); iPoint++ )
   {
    const Point<3U> pt( e.RstToXYZ(e.FV()->FacetPoint(iFacet,iPoint)) );
    ss << pt[0] << "," << pt[1] << "," << pt[2] << " ";
   }
   
   ss << "_Enter ";
 }
 

static void writeNormal (std::stringstream & ss, const EFT3& efvt, const Element<3U>& e, uint32_t iFacet, bool bInverted)
 {
   Point<3U> pt1( e.RstToXYZ(e.FV()->FacetIntegrationPoint(iFacet,0U)) );
   Point<3U> pt2;
   
   if(bInverted)
    pt2 = pt1 - efvt.FacetNormal( iFacet );
   else
    pt2 = pt1 + efvt.FacetNormal( iFacet );
   
   ss << "\nPolyline ";
   
   ss << pt1[0] << "," << pt1[1] << "," << pt1[2] << " ";
   ss << pt2[0] << "," << pt2[1] << "," << pt2[2] << " ";   
   
   ss << "_Enter ";
 }
  
// prints element as combination of surfaces and integration points 
void printElement( const Element<3U>& e, const string& sNameOfFile )
 {
  //draw element
  ofstream ofs(sNameOfFile.c_str(), ios::out|ios::app);
  
  //print faces
  {
  stringstream ss;
  createLayer (ss, "Faces", e.Idx());
  ofs << ss.str();
  }
  
  ofs << "\n;*****FACES ";
  vector<uint32_t> fnids;
  for ( auto iFace = 0U; iFace < e.Faces(); iFace++ )
  {
   e.FE()->NodesOfFace( iFace, fnids );
   ofs << "\nSrfPt ";
   for ( size_t iNode = 0U; iNode < fnids.size(); iNode++ )
   {
     ofs << e.N(fnids[iNode])->x() << ",";
     ofs << e.N(fnids[iNode])->y() << ",";
     ofs << e.N(fnids[iNode])->z() << " ";
   }
   ofs << "_Enter ";
  }
  
  //change layer
  ofs << "\n;*****INTEGRATION_POINTS ";
  
  {
  stringstream ss;
  createLayer (ss, "ip", e.Idx());
  ofs << ss.str();
  }
  
  //print facets
  for ( auto ip=0U; ip < e.FE()->IntegrationPoints(); ip++ )
  {
     ofs << "\nPoint ";
     Point<3U> pt((e).IntegrationPoint(ip));
     ofs << pt[0] << "," << pt[1] << "," << pt[2] << " ";
     ofs << "_Enter ";
   }
       
  ofs.close();
  
 } // end printElement

 void printElementStencil( const Element<3U>& e, const string& sNameOfFile )
 {
  //draw element
  ofstream ofs(sNameOfFile.c_str(), ios::out|ios::app);
  
  //print faces
  {
  stringstream ss;
  createLayer (ss, "Faces", e.Idx());
  ofs << ss.str();
  }
  
  ofs << "\n;*****FACES ";
  vector<uint32_t> fnids;
  for ( auto iFace = 0U; iFace < e.Faces(); iFace++ )
  {
   e.FE()->NodesOfFace( iFace, fnids );
   ofs << "\nSrfPt ";
   for ( size_t iNode = 0U; iNode < fnids.size(); iNode++ )
   {
     ofs << e.N(fnids[iNode])->x() << ",";
     ofs << e.N(fnids[iNode])->y() << ",";
     ofs << e.N(fnids[iNode])->z() << " ";
   }
   ofs << "_Enter ";
  }
  
  //change layer
  ofs << "\n;*****FACETS " << e.FV()->Facets();
  
  {
  stringstream ss;
  createLayer (ss, "Facets", e.Idx());
  ofs << ss.str();
  }
  
  //print facets
  for ( auto iFacet = 0U; iFacet < e.FV()->Facets(); iFacet++ )
  {
   if(e.FE()->IsSurface()) // only two facet points
   {
     ofs << "\nPolyline ";
     Point<3U> pt1( e.RstToXYZ(e.FV()->FacetPoint(iFacet,0U)) );
     Point<3U> pt2( e.RstToXYZ(e.FV()->FacetPoint(iFacet,1U)) );
     ofs << pt1[0] << "," << pt1[1] << "," << pt1[2] << " ";
     ofs << pt2[0] << "," << pt2[1] << "," << pt2[2] << " ";
     ofs << "_Enter ";
   }
   else if(e.FE()->IsVolume()) // there are four facet points
   {
     stringstream ss;
     writeSurfaceFacet(ss, e, e, iFacet); 
     ofs << ss.str();
   }
   
  }
  
  //change layer
  {
  stringstream ss;
  createLayer (ss, "Normals", e.Idx());
  ofs << ss.str();
  }
  
  ofs << "\n; *****NORMALS, out of the FACET Integration Points ";
  for ( auto iFacet = 0U; iFacet < e.FV()->Facets(); iFacet++ )
  {
   stringstream ss;
   writeNormal(ss, e, e, iFacet, false);
   ofs << ss.str();
  }
  
  ofs << "\n\n-Layer Current Default _Enter \n\n";
  ofs << "\n\n-Layer Off "  << e.Idx() << "Normals" << " _Enter \n\n";
  
  ofs.close();
 } // end printElementStencil
 


 void createLayer ( std::stringstream & ss, const std::string & layer, size_t id_element )
 {
  stringstream layer_color;
  layer_color << rand()%200 << "," << rand()%200 << "," << rand()%200;
  
  ss << "\n\n-Layer New "  << id_element << layer << " _Enter \n\n";
  ss << "\n-Layer Color "  << id_element << layer << " " << layer_color.str() << " _Enter \n";
  ss << "\n-Layer Current "<< id_element << layer << " _Enter \n";
  ss << "\n-Layer Material "<< id_element << layer << " Color " << layer_color.str() << " _Enter _Enter \n";
 }
 

 
 //only prints volumetric finite volumes
 void printFiniteVolumes( const Model<3U>& sg )
 {    
   uint32_t inside_node, outside_node;
    
   typedef std::vector< Point<3U> > FacetPoints;
   typedef std::pair<FacetPoints, Point<3U> > Facet;
   typedef std::vector< Facet > Node;
   
   vector< string > nodes (sg.Region("Model").Nodes());
   
   for ( auto i = 0; i < nodes.size(); i++ )
   {
    stringstream s;
    createLayer( s, "FV", i );
    nodes[i]=s.str();
   }
   
   const Region<3U>&  sgref(sg.Region("Model"));
   
   for ( auto eit=sgref.CellsBegin(); eit!=sgref.CellsEnd(); eit++ ) 
     for ( auto iFacet=0U; iFacet<(*eit)->FV()->Facets(); iFacet++ )
  	  if((*eit)->FE()->IsVolume()) // there are four facet points
        {    	    
    	    (*eit)->FV()->FacetEdgeNodes( iFacet, inside_node, outside_node );
  	      
    	    const size_t inode = (*eit)->N(inside_node)->Idx(); 
          const size_t onode = (*eit)->N(outside_node)->Idx(); 
          
          stringstream in_s,out_s;
          
    	    
    	    //get the facets
          writeSurfaceFacet(in_s, *(*eit), *(*eit), iFacet);
          writeSurfaceFacet(out_s, *(*eit), *(*eit), iFacet);
          
          //get the normals
          writeNormal(in_s, *(*eit), *(*eit), iFacet, false);
          writeNormal(out_s, *(*eit), *(*eit), iFacet, true);
          
          nodes[inode]+=in_s.str();
          nodes[onode]+=out_s.str();
          
         }
 
  ofstream ofs("finite_volumes.txt", ios::out|ios::trunc);
  
	for ( auto i = 0; i < nodes.size(); i++ )
	  ofs << nodes[i] << "\n";
  
  ofs.close();
    
  }
} // end csmp
