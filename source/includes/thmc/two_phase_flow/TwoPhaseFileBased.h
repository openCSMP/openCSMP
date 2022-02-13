#ifndef TWOPHASEFILEBASED_H
#define TWOPHASEFILEBASED_H

#include "TwoPhaseModel.h"
#include "CubicSpline.h"

namespace csmp {

/// Relative Permeability model based on text file input of Kr and Pc vs. Seff
template<uint32_t dim>
class TwoPhaseFileBased : public TwoPhaseModel<dim> {
  public:
    explicit TwoPhaseFileBased( const char* fileName );
    TwoPhaseFileBased( const PropertyDatabase<dim>&, const char* fileName );
    TwoPhaseFileBased( const PropertyDatabase<dim>&, const std::vector<double>& seff,
                                                const std::vector<double>& krn_Phase,
                                                const std::vector<double>& krw_Phase,
                                                const std::vector<double>& pc_Phase  );

    virtual ~TwoPhaseFileBased();

    void ReadFile( const char* fileName );

    virtual void Initialize( const Element<dim>& );

    // relative permeabilities
    virtual double krn_Phase() const;
    virtual double krw_Phase() const;

    // derivatives of relative permeabilities
    virtual double dkrnds_Phase() const;
    virtual double dkrwds_Phase() const;

    // derivative of fractional flow (advection multipliers)
    virtual double dfds() const;

    // maximum absolute value returned by dfds
    virtual double MaxFractionalFlowDerivative() const;

    // derivatives of rel perm (advection multipliers)
    virtual double dGds( ) const;

    // capillary pressure
    virtual double pc_Phase( ) const;

    // capillary pressure derivatives
    virtual double dpcds_Phase( ) const;

    // inverse capillary pressure function
    virtual double Sw_Phase( double pc_Phase ) const;

    // inverse capillary pressure derivative
    virtual double dsdpc_Phase( double pc_Phase ) const;

    int32_t             writeData();

  private:
    TwoPhaseFileBased();
    std::streampos    findPosition( std::ifstream& ) const;
    int32_t             readData();

    std::string       catchPhrase_;
    std::ifstream     relpermFile_;

    std::vector<double>   seff_;
    std::vector<double>   krw_;
    std::vector<double>   krn_;
    std::vector<double>   pc_;

    const double    GRAVITY_,
                      LOWER_LIMIT_,
                      MAX_CAPILLARY_PRESSURE_,
                      MAX_CAPILLARY_PRESSURE_SLOPE_;

    CubicSpline       krwCurve_;
    CubicSpline       krnCurve_;
    CubicSpline       pcCurve_;
};

/**
@class TwoPhaseFileBased TwoPhaseFileBased "two_phase_flow/TwoPhaseFileBased.h"
@author P. Lang
@author Georg Seidl
@date May 2010

@section implementation Implementation
Uses csmp::CubicSpline to interpolate and compute derivatives
Maybe consider a local interpolation method and sacrifice derivatives for speed??

@section application Application
Use RelPerm Excel Workbook/GUI to create input file.

@code
TwoPhaseFileBased<DIM> relperm_model( model.Database(), "RelPerms.txt" );
@endcode

@section input Input File
This model requires a input file in the following format

@code
...SOME TEXT...

seff      krw     krn     pc
0 .       0       1       400000
...
.4        .5      .7      30000
...
1.        1.      0.      200
@endcode

*/

} // end namespace csmp

#endif // TWOPHASEFILEBASED_H
