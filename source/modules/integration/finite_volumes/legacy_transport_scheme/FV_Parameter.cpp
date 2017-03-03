#include "FV_Parameter.h"

using namespace std;

namespace csmp {


FV_Parameter::FV_Parameter( size_t sectors, size_t facets, size_t dim, bool with_normals )
 : sector_volume_( sectors ),
   facet_v_and_A_( facets )
 {
    if ( !with_normals ) return;
    
    facet_unit_normal_.resize( facets );
    vector<vector<double64> >(facet_unit_normal_).swap(facet_unit_normal_);
    for ( auto it : facet_unit_normal_ ) {
         it.resize(dim); // n-dimensional
         vector<double64>( it ).swap( it );
      }
 }




void FV_Parameter::Resize( size_t sectors, size_t facets, size_t dim, bool with_normals )
 {
    sector_volume_.resize( sectors );  
    vector<double64>(sector_volume_).swap(sector_volume_);
    
    facet_v_and_A_.resize( facets );   
    vector<pair<double64,double64> >(facet_v_and_A_).swap(facet_v_and_A_);
    
    if ( with_normals ) {
          facet_unit_normal_.resize( facets );
          vector<vector<double64> >(facet_unit_normal_).swap(facet_unit_normal_);
          for ( auto it : facet_unit_normal_ ) {
               it.resize(dim); // n-dimensional
               vector<double64>( it ).swap( it );
            }
      }
 }




size_t  FV_Parameter::Bytes() const 
 {
    size_t bytes = (sector_volume_.size() + facet_v_and_A_.size()) * sizeof(double64) + sizeof(size_t);
    return (facet_unit_normal_.empty()) ? bytes : bytes + facet_unit_normal_.size() * facet_unit_normal_[0].size() * sizeof(double64); 
 }




/// dot product fn . vc
double64 FV_Parameter::FacetNormalProjection( size_t facet, 
                                               const std::vector<double64>& cxyz ) const
 {
    if ( facet_unit_normal_.empty() or 
         facet_unit_normal_[facet].empty() or 
         cxyz.empty() or 
         cxyz.size() != facet_unit_normal_[facet].size() ) {
         std::cout <<"\nFV_Parameter::FacetNormalProjection: incompatible vector sizes: "<< endl;
         cout <<"\nsizes of unit normal vector in stencil: "<< facet_unit_normal_.size();
         if ( !facet_unit_normal_.empty() ) cout <<", unit normal entry "<< facet <<": ";
         cout << facet_unit_normal_[facet].size();
         cout <<", size of input vector: "<< cxyz.size() << endl;
         if ( !facet_unit_normal_.empty() and !facet_unit_normal_[facet].empty() ) {
              cout <<"unit normal vector: "; 
              for ( auto fit : facet_unit_normal_[facet] )
                std::cout << fit <<" ";
           }
         std::cout << std::endl;
         return 0.; 
      }
 
    std::vector<double64>::const_iterator fit=facet_unit_normal_[facet].begin();
    std::vector<double64>::const_iterator it=cxyz.begin();
    double64  dotpr(0.);

    while ( it != cxyz.end() ) {
         dotpr += static_cast<double64>(*fit) * static_cast<double64>(*it);
         fit++;
         it++;
      }
    return dotpr;
 }




void FV_Parameter::Out(std::ostream& os) const
 {
    os <<"\nFV_Parameter::Out: Data of finite-volume stencil: ";

    os <<"area, unit normal, and velocity magnitude on the finite volume facets: "<< endl;
    for ( size_t i=0; i<facet_v_and_A_.size(); i++ ) {
         //      facet                        area
         os <<"f"<< i+1 <<": A "<< facet_v_and_A_[i].second <<", fn";
         // unit normal
         for ( size_t j=0; j<facet_unit_normal_[i].size(); j++ )
           os <<" "<< facet_unit_normal_[i][j];
         os <<", v ";
         // velocity
         os << facet_v_and_A_[i].first << endl;
      }

    os <<"sector volumes: "<< endl;
    for ( size_t i=0; i<sector_volume_.size(); i++ )
      os << i+1 <<": "<< sector_volume_[i] << endl;
    os << endl;
 } // end out


} // end namespace csmp















