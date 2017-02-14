#ifndef CSMP_FV_PARAMETER_H
#define CSMP_FV_PARAMETER_H

#include "VectorVariable.h"

namespace csmp {

/// storage for the sector volumes, facet area and projected velocities
/// so that these do not have to be recompute each timecrement
class FV_Parameter {
  public:
    FV_Parameter() {};
    FV_Parameter( size_t sectors, size_t facets, size_t dim, bool with_normals=false );
    FV_Parameter( const FV_Parameter& );
    FV_Parameter& operator=( const FV_Parameter& );
    FV_Parameter( FV_Parameter&& );
    FV_Parameter& operator=( FV_Parameter&& );
  
    // mutators
    void Resize( size_t sectors, size_t facets, size_t dim, bool with_normals=false );
    void SectorVolume( size_t sector, double64 vol );
    void FacetArea( size_t facet, double64 area );
    void FacetNormal( size_t facet, const std::vector<double64>& fnxyz );
    void FacetNormalVelocity( size_t facet, double64 flux );
    void Initialize( size_t facet, double64 flux, double64 area );
    
    // accessors
    size_t Sectors() const;
    size_t Facets() const;
    double64 SectorVolume( size_t sector ) const;
    double64 FacetArea( size_t facet ) const;
    double64 FacetNormalVelocity( size_t facet ) const;
    double64 FacetNormalComponent( size_t facet, size_t x_or_y_or_z ) const;
    double64 FacetNormalProjection( size_t facet, const std::vector<double64>& cxyz ) const;
    // universal versions
    double64 FacetNormalProjection( size_t facet, const VectorVariable<1U>& cxyz ) const;
    double64 FacetNormalProjection( size_t facet, const VectorVariable<2U>& cxyz ) const;
    double64 FacetNormalProjection( size_t facet, const VectorVariable<3U>& cxyz ) const;

    size_t Bytes() const;
    void Out() const { Out(std::cout); }
    void   Out(std::ostream& os) const;
    
  private:
    ///< velocites projected on facet normals and facet areas
    std::vector<double64>                       sector_volume_;
    ///< velocites projected on facet normals and facet areas
    std::vector<std::pair<double64,double64> >  facet_v_and_A_;
    std::vector<std::vector<double64> >         facet_unit_normal_;
};


inline FV_Parameter::FV_Parameter( const FV_Parameter& param )
 : sector_volume_(param.sector_volume_), 
   facet_v_and_A_(param.facet_v_and_A_),
   facet_unit_normal_(param.facet_unit_normal_) 
 {
 }


inline FV_Parameter::FV_Parameter( FV_Parameter&& param )
 : sector_volume_{param.sector_volume_},
   facet_v_and_A_{param.facet_v_and_A_},
   facet_unit_normal_{param.facet_unit_normal_}
 {
 }


inline FV_Parameter& FV_Parameter::operator=( const FV_Parameter& param )
 {
    if ( &param != this ) {
         sector_volume_     = param.sector_volume_;
         facet_v_and_A_     = param.facet_v_and_A_;
         facet_unit_normal_ = param.facet_unit_normal_;
      }
    return *this;
 }


inline FV_Parameter& FV_Parameter::operator=( FV_Parameter&& param )
 {
    if ( &param != this ) {
         sector_volume_     = {param.sector_volume_};
         facet_v_and_A_     = {param.facet_v_and_A_};
         facet_unit_normal_ = {param.facet_unit_normal_};
      }
    return *this;
 }


inline void FV_Parameter::SectorVolume( size_t sector, double64 vol )
 { sector_volume_[sector] = vol; }


inline void FV_Parameter::FacetNormalVelocity( size_t facet, double64 flux )
 { facet_v_and_A_[facet].first = flux; }


inline void FV_Parameter::FacetArea( size_t facet, double64 area )
 { facet_v_and_A_[facet].second = area; }


// assign the facet unit normal for the facet
inline void FV_Parameter::FacetNormal( size_t facet, const std::vector<double64>& fnxyz )
 {
    facet_unit_normal_[facet] = fnxyz;
 }

inline void FV_Parameter::Initialize( size_t facet, double64 flux, double64 area )
 { 
    facet_v_and_A_[facet].first  = flux;
    facet_v_and_A_[facet].second = area; 
 }


inline size_t FV_Parameter::Sectors() const
 { return sector_volume_.size(); }
 
 
inline size_t FV_Parameter::Facets() const
 { return facet_v_and_A_.size(); }


inline double64 FV_Parameter::SectorVolume( size_t sector ) const
 { return sector_volume_[sector]; }


inline double64 FV_Parameter::FacetNormalVelocity( size_t facet ) const
 { return facet_v_and_A_[facet].first; }


inline double64 FV_Parameter::FacetNormalComponent( size_t facet, size_t x_or_y_or_z ) const
 {
    return facet_unit_normal_[facet][x_or_y_or_z];
 }


/// universal versions of projection functions
inline double64 FV_Parameter::FacetNormalProjection( size_t facet, const VectorVariable<1U>& vc ) const
 {
    // dot product fn . vc
    return facet_unit_normal_[facet][0] * vc[0];
 }

inline double64 FV_Parameter::FacetNormalProjection( size_t facet, const VectorVariable<2U>& vc ) const
 {
    // dot product fn . vc
    return facet_unit_normal_[facet][0] * vc[0] + facet_unit_normal_[facet][1] * vc[1];
 }

inline double64 FV_Parameter::FacetNormalProjection( size_t facet, const VectorVariable<3U>& vc ) const
 {
    // dot product fn . vc
    return facet_unit_normal_[facet][0] * vc[0] + facet_unit_normal_[facet][1] * vc[1] + facet_unit_normal_[facet][2] * vc[2];
 }



inline double64 FV_Parameter::FacetArea( size_t facet ) const
 { return facet_v_and_A_[facet].second; }


} // end namespace csmp

#endif
