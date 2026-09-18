#include "FV_Parameter_Test.h"
#include "FV_Parameter.h"
#include "VectorVariable.h"

namespace csmp {



FV_Parameter_Test::FV_Parameter_Test()
{
  fTolerance = 1.e-20;
}


FV_Parameter_Test::~FV_Parameter_Test()
{
}


void FV_Parameter_Test::run()
{
    FV_ParameterCtor();
    FV_ParameterCopyCtor();
    FV_ParameterEqual();
    FV_ParameterInitialize();
    FV_ParameterResize();
    FV_ParameterSectorVolume();
    FV_ParameterFacetArea();
    FV_ParameterFacetNormal();
    FV_ParameterFacetNormalVelocity();
    FV_ParameterFacetNormalProjection();

}


void FV_Parameter_Test::FV_ParameterCtor()
{
    size_t sectors(2U);
    size_t facets(4U);
    uint32_t dim(2U);
    
    FV_Parameter fv_parameter1;
    FV_Parameter fv_parameter2( sectors, facets, dim, true );
    
    _test( fv_parameter1.Sectors() == 0U );
    _test( fv_parameter1.Facets() == 0U );
    _test( fv_parameter2.Sectors() == sectors );
    _test( fv_parameter2.Facets() == facets );
    
}



void FV_Parameter_Test::FV_ParameterCopyCtor()
{
    size_t sectors(2U);
    size_t facets(4U);
    uint32_t dim(2U);
    double sec_vol0(10.5);
    double sec_vol1(2.7);
    double fct_area0(0.8);
    double fct_area1(2.8);
    double fct_area2(5.8);
    double fct_area3(7.7);
    double fct_nvel0(-1.2);
    double fct_nvel1(3.2);
    double fct_nvel2(4.);
    double fct_nvel3(-7.);
    std::vector<double> fct_normal0(dim);
    std::vector<double> fct_normal1(dim);
    std::vector<double> fct_normal2(dim);
    std::vector<double> fct_normal3(dim);
    
    fct_normal0[0] = -1.;
    fct_normal0[1] = 2.;
    fct_normal1[0] = 3.;
    fct_normal1[1] = -1.;
    fct_normal2[0] = -1.5;
    fct_normal2[1] = -2.4;
    fct_normal3[0] = 1.5;
    fct_normal3[1] = 2.2;

    FV_Parameter fv_parameter1( sectors, facets, dim, true );
    
    fv_parameter1.SectorVolume(0U, sec_vol0);
    fv_parameter1.SectorVolume(1U, sec_vol1);
    fv_parameter1.FacetArea(0U, fct_area0);
    fv_parameter1.FacetArea(1U, fct_area1);
    fv_parameter1.FacetArea(2U, fct_area2);
    fv_parameter1.FacetArea(3U, fct_area3);
    fv_parameter1.FacetNormalVelocity(0U, fct_nvel0);
    fv_parameter1.FacetNormalVelocity(1U, fct_nvel1);
    fv_parameter1.FacetNormalVelocity(2U, fct_nvel2);
    fv_parameter1.FacetNormalVelocity(3U, fct_nvel3);
    fv_parameter1.FacetNormal(0U, fct_normal0);
    fv_parameter1.FacetNormal(1U, fct_normal1);
    fv_parameter1.FacetNormal(2U, fct_normal2);
    fv_parameter1.FacetNormal(3U, fct_normal3);
    
    FV_Parameter fv_parameter2(fv_parameter1);
    
    _test( fv_parameter2.Sectors() == sectors );
    _test( fv_parameter2.Facets() == facets );
    
    _equal( fv_parameter2.SectorVolume(0U), sec_vol0, fTolerance );
    _equal( fv_parameter2.SectorVolume(1U), sec_vol1, fTolerance );
    _equal( fv_parameter2.FacetArea(0U), fct_area0, fTolerance );
    _equal( fv_parameter2.FacetArea(1U), fct_area1, fTolerance );
    _equal( fv_parameter2.FacetArea(2U), fct_area2, fTolerance );
    _equal( fv_parameter2.FacetArea(3U), fct_area3, fTolerance );
    _equal( fv_parameter2.FacetNormalVelocity(0U), fct_nvel0, fTolerance );
    _equal( fv_parameter2.FacetNormalVelocity(1U), fct_nvel1, fTolerance );
    _equal( fv_parameter2.FacetNormalVelocity(2U), fct_nvel2, fTolerance );
    _equal( fv_parameter2.FacetNormalVelocity(3U), fct_nvel3, fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(0U, 0U), fct_normal0[0], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(0U, 1U), fct_normal0[1], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(1U, 0U), fct_normal1[0], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(1U, 1U), fct_normal1[1], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(2U, 0U), fct_normal2[0], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(2U, 1U), fct_normal2[1], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(3U, 0U), fct_normal3[0], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(3U, 1U), fct_normal3[1], fTolerance );

}


void FV_Parameter_Test::FV_ParameterEqual()
{
    size_t sectors(2U);
    size_t facets(4U);
    uint32_t dim(2U);
    double sec_vol0(10.5);
    double sec_vol1(2.7);
    double fct_area0(0.8);
    double fct_area1(2.8);
    double fct_area2(5.8);
    double fct_area3(7.7);
    double fct_nvel0(-1.2);
    double fct_nvel1(3.2);
    double fct_nvel2(4.);
    double fct_nvel3(-7.);
    std::vector<double> fct_normal0(dim);
    std::vector<double> fct_normal1(dim);
    std::vector<double> fct_normal2(dim);
    std::vector<double> fct_normal3(dim);
    
    fct_normal0[0] = -1.;
    fct_normal0[1] = 2.;
    fct_normal1[0] = 3.;
    fct_normal1[1] = -1.;
    fct_normal2[0] = -1.5;
    fct_normal2[1] = -2.4;
    fct_normal3[0] = 1.5;
    fct_normal3[1] = 2.2;

    FV_Parameter fv_parameter1( sectors, facets, dim, true );
    FV_Parameter fv_parameter2;
    
    fv_parameter1.SectorVolume(0U, sec_vol0);
    fv_parameter1.SectorVolume(1U, sec_vol1);
    fv_parameter1.FacetArea(0U, fct_area0);
    fv_parameter1.FacetArea(1U, fct_area1);
    fv_parameter1.FacetArea(2U, fct_area2);
    fv_parameter1.FacetArea(3U, fct_area3);
    fv_parameter1.FacetNormalVelocity(0U, fct_nvel0);
    fv_parameter1.FacetNormalVelocity(1U, fct_nvel1);
    fv_parameter1.FacetNormalVelocity(2U, fct_nvel2);
    fv_parameter1.FacetNormalVelocity(3U, fct_nvel3);
    fv_parameter1.FacetNormal(0U, fct_normal0);
    fv_parameter1.FacetNormal(1U, fct_normal1);
    fv_parameter1.FacetNormal(2U, fct_normal2);
    fv_parameter1.FacetNormal(3U, fct_normal3);
    
    fv_parameter2 = fv_parameter1;
    
    _test( fv_parameter2.Sectors() == sectors );
    _test( fv_parameter2.Facets() == facets );
    
    _equal( fv_parameter2.SectorVolume(0U), sec_vol0, fTolerance );
    _equal( fv_parameter2.SectorVolume(1U), sec_vol1, fTolerance );
    _equal( fv_parameter2.FacetArea(0U), fct_area0, fTolerance );
    _equal( fv_parameter2.FacetArea(1U), fct_area1, fTolerance );
    _equal( fv_parameter2.FacetArea(2U), fct_area2, fTolerance );
    _equal( fv_parameter2.FacetArea(3U), fct_area3, fTolerance );
    _equal( fv_parameter2.FacetNormalVelocity(0U), fct_nvel0, fTolerance );
    _equal( fv_parameter2.FacetNormalVelocity(1U), fct_nvel1, fTolerance );
    _equal( fv_parameter2.FacetNormalVelocity(2U), fct_nvel2, fTolerance );
    _equal( fv_parameter2.FacetNormalVelocity(3U), fct_nvel3, fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(0U, 0U), fct_normal0[0], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(0U, 1U), fct_normal0[1], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(1U, 0U), fct_normal1[0], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(1U, 1U), fct_normal1[1], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(2U, 0U), fct_normal2[0], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(2U, 1U), fct_normal2[1], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(3U, 0U), fct_normal3[0], fTolerance );
    _equal( fv_parameter2.FacetNormalComponent(3U, 1U), fct_normal3[1], fTolerance );

}


void FV_Parameter_Test::FV_ParameterInitialize()
{
    size_t sectors(2U);
    size_t facets(4U);
    uint32_t dim(2U);
    double fct_area0(0.8);
    double fct_area1(2.8);
    double fct_area2(5.8);
    double fct_area3(7.7);
    double fct_nvel0(-1.2);
    double fct_nvel1(3.2);
    double fct_nvel2(4.);
    double fct_nvel3(-7.);

    FV_Parameter fv_parameter1( sectors, facets, dim );
    
    fv_parameter1.Initialize(0U, fct_nvel0, fct_area0);
    fv_parameter1.Initialize(1U, fct_nvel1, fct_area1);
    fv_parameter1.Initialize(2U, fct_nvel2, fct_area2);
    fv_parameter1.Initialize(3U, fct_nvel3, fct_area3);

    _equal( fv_parameter1.FacetArea(0U), fct_area0, fTolerance );
    _equal( fv_parameter1.FacetArea(1U), fct_area1, fTolerance );
    _equal( fv_parameter1.FacetArea(2U), fct_area2, fTolerance );
    _equal( fv_parameter1.FacetArea(3U), fct_area3, fTolerance );
    _equal( fv_parameter1.FacetNormalVelocity(0U), fct_nvel0, fTolerance );
    _equal( fv_parameter1.FacetNormalVelocity(1U), fct_nvel1, fTolerance );
    _equal( fv_parameter1.FacetNormalVelocity(2U), fct_nvel2, fTolerance );
    _equal( fv_parameter1.FacetNormalVelocity(3U), fct_nvel3, fTolerance );

}


void FV_Parameter_Test::FV_ParameterResize()
{
    size_t sectors1(2U);
    size_t facets1(4U);
    uint32_t dim1(2U);
    size_t sectors2(4U);
    size_t facets2(8U);
    uint32_t dim2(3U);
    
    FV_Parameter fv_parameter1( sectors1, facets1, dim1 );
    fv_parameter1.Resize( sectors2, facets2, dim2, true );
    
    _test( fv_parameter1.Sectors() == sectors2 );
    _test( fv_parameter1.Facets() == facets2 );

    fv_parameter1.Resize( sectors1, facets1, dim1, true );
    
    _test( fv_parameter1.Sectors() == sectors1 );
    _test( fv_parameter1.Facets() == facets1 );

}


void FV_Parameter_Test::FV_ParameterSectorVolume()
{
    size_t sectors(2U);
    size_t facets(4U);
    uint32_t dim(2U);
    double sec_vol0(10.5);
    double sec_vol1(2.7);

    FV_Parameter fv_parameter1( sectors, facets, dim, true );
    
    fv_parameter1.SectorVolume(0U, sec_vol0);
    fv_parameter1.SectorVolume(1U, sec_vol1);
    
    _equal( fv_parameter1.SectorVolume(0U), sec_vol0, fTolerance );
    _equal( fv_parameter1.SectorVolume(1U), sec_vol1, fTolerance );

}


void FV_Parameter_Test::FV_ParameterFacetArea()
{
    size_t sectors(2U);
    size_t facets(4U);
    uint32_t dim(2U);
    double fct_area0(0.8);
    double fct_area1(2.8);
    double fct_area2(5.8);
    double fct_area3(7.7);

    FV_Parameter fv_parameter1( sectors, facets, dim, true );
    
    fv_parameter1.FacetArea(0U, fct_area0);
    fv_parameter1.FacetArea(1U, fct_area1);
    fv_parameter1.FacetArea(2U, fct_area2);
    fv_parameter1.FacetArea(3U, fct_area3);
    
    _equal( fv_parameter1.FacetArea(0U), fct_area0, fTolerance );
    _equal( fv_parameter1.FacetArea(1U), fct_area1, fTolerance );
    _equal( fv_parameter1.FacetArea(2U), fct_area2, fTolerance );
    _equal( fv_parameter1.FacetArea(3U), fct_area3, fTolerance );

}


void FV_Parameter_Test::FV_ParameterFacetNormal()
{
    size_t sectors(2U);
    size_t facets(4U);
    uint32_t dim(2U);

    std::vector<double> fct_normal0(dim);
    std::vector<double> fct_normal1(dim);
    std::vector<double> fct_normal2(dim);
    std::vector<double> fct_normal3(dim);
    
    fct_normal0[0] = -1.;
    fct_normal0[1] = 2.;
    fct_normal1[0] = 3.;
    fct_normal1[1] = -1.;
    fct_normal2[0] = -1.5;
    fct_normal2[1] = -2.4;
    fct_normal3[0] = 1.5;
    fct_normal3[1] = 2.2;

    FV_Parameter fv_parameter1( sectors, facets, dim, true );
    
    fv_parameter1.FacetNormal(0U, fct_normal0);
    fv_parameter1.FacetNormal(1U, fct_normal1);
    fv_parameter1.FacetNormal(2U, fct_normal2);
    fv_parameter1.FacetNormal(3U, fct_normal3);
    
    _equal( fv_parameter1.FacetNormalComponent(0U, 0U), fct_normal0[0], fTolerance );
    _equal( fv_parameter1.FacetNormalComponent(0U, 1U), fct_normal0[1], fTolerance );
    _equal( fv_parameter1.FacetNormalComponent(1U, 0U), fct_normal1[0], fTolerance );
    _equal( fv_parameter1.FacetNormalComponent(1U, 1U), fct_normal1[1], fTolerance );
    _equal( fv_parameter1.FacetNormalComponent(2U, 0U), fct_normal2[0], fTolerance );
    _equal( fv_parameter1.FacetNormalComponent(2U, 1U), fct_normal2[1], fTolerance );
    _equal( fv_parameter1.FacetNormalComponent(3U, 0U), fct_normal3[0], fTolerance );
    _equal( fv_parameter1.FacetNormalComponent(3U, 1U), fct_normal3[1], fTolerance );

}


void FV_Parameter_Test::FV_ParameterFacetNormalVelocity()
{
    size_t sectors(2U);
    size_t facets(4U);
    uint32_t dim(2U);
    double fct_nvel0(-1.2);
    double fct_nvel1(3.2);
    double fct_nvel2(4.);
    double fct_nvel3(-7.);

    FV_Parameter fv_parameter1( sectors, facets, dim, true );

    fv_parameter1.FacetNormalVelocity(0U, fct_nvel0);
    fv_parameter1.FacetNormalVelocity(1U, fct_nvel1);
    fv_parameter1.FacetNormalVelocity(2U, fct_nvel2);
    fv_parameter1.FacetNormalVelocity(3U, fct_nvel3);

    _equal( fv_parameter1.FacetNormalVelocity(0U), fct_nvel0, fTolerance );
    _equal( fv_parameter1.FacetNormalVelocity(1U), fct_nvel1, fTolerance );
    _equal( fv_parameter1.FacetNormalVelocity(2U), fct_nvel2, fTolerance );
    _equal( fv_parameter1.FacetNormalVelocity(3U), fct_nvel3, fTolerance );

}


void FV_Parameter_Test::FV_ParameterFacetNormalProjection()
{
    VectorVariable<1U> v1(PLAIN, 3);
    VectorVariable<2U> v2(PLAIN, PLAIN, 3., 2.);
    VectorVariable<3U> v3(PLAIN, PLAIN, PLAIN, 2., 5., -1.);
    
    std::vector<double> vec1(1U);
    std::vector<double> vec2(2U);
    std::vector<double> vec3(3U);

    vec1[0] = 3.;
    vec2[0] = 3.;
    vec2[1] = 2.;
    vec3[0] = 2.;
    vec3[1] = 5.;
    vec3[2] = -1.;

    size_t sectors(1U);
    size_t facets(1U);

    std::vector<double> fct_normal1(1U);
    std::vector<double> fct_normal2(2U);
    std::vector<double> fct_normal3(3U);

    fct_normal1[0] = 3.;
    fct_normal2[0] = -1.;
    fct_normal2[1] = -2.;
    fct_normal3[0] = 2.;
    fct_normal3[1] = -4.;
    fct_normal3[2] = 3.;

    FV_Parameter fv_parameter1( sectors, facets, 1U, true );
    FV_Parameter fv_parameter2( sectors, facets, 2U, true );
    FV_Parameter fv_parameter3( sectors, facets, 3U, true );
    
    fv_parameter1.FacetNormal(0U, fct_normal1);
    fv_parameter2.FacetNormal(0U, fct_normal2);
    fv_parameter3.FacetNormal(0U, fct_normal3);

    _equal( fv_parameter1.FacetNormalProjection(0U, v1), 9., fTolerance );
    _equal( fv_parameter2.FacetNormalProjection(0U, v2), -7., fTolerance );
    _equal( fv_parameter3.FacetNormalProjection(0U, v3), -19., fTolerance );

    _equal( fv_parameter1.FacetNormalProjection(0U, vec1), 9., fTolerance );
    _equal( fv_parameter2.FacetNormalProjection(0U, vec2), -7., fTolerance );
    _equal( fv_parameter3.FacetNormalProjection(0U, vec3), -19., fTolerance );

}



} // end namespace csmp
