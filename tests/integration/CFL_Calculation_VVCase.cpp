#include "CFL_Calculation_VVCase.h"

using namespace std;

namespace csmp {

template <uint32_t dim>
CFL_Calculation_VVCase<dim>::CFL_Calculation_VVCase()
{
}

template <uint32_t dim>
CFL_Calculation_VVCase<dim>::~CFL_Calculation_VVCase()
{
}


template <uint32_t dim>
CFL_Calculation_VVCase<dim>::CFL_Calculation_VVCase(const char* prefix)
{
    this->setName("SlightlyCompressibleSinglePhaseFlow2D_VVCase");
    prefix_=prefix;
}

template <uint32_t dim>
void CFL_Calculation_VVCase<dim>::run()
{
    if (dim!=0)
        cout<<"CFL_Calculation_VVCase::Test not fully implemented yet!"<<endl;

    exit(1);
}

template class CFL_Calculation_VVCase<1U>;
template class CFL_Calculation_VVCase<2U>;
template class CFL_Calculation_VVCase<3U>;
}
