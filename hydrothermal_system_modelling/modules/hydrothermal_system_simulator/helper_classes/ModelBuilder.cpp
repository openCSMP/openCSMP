// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ModelBuilder.h"

using namespace std;

namespace csmp
{

template<size_t dim>
ModelBuilder<dim>::ModelBuilder()
    : var_file("PhysicalVariables.txt")
{
}


template<size_t dim>
void ModelBuilder<dim>::CreateVariablesFile(bool with_gold,
                                            bool with_lithium,
                                            bool with_flux_visitor,
                                            bool with_tracer,
                                            bool with_degassing,
                                            bool with_split_boundary)
{

    std::ofstream variables_file(var_file.c_str());

    //! General variables
    variables_file << "name\tMParameter\tunit\tindex\tmin.\tmax.\tplace" << std::endl;
    variables_file << "element ID\teID\t-\t1\t0.00E0\t1.00E+10\telement" << std::endl;
    variables_file << "node ID\tnID\t-\t1\t0.00E0\t1.00E+10\tnode" << std::endl;
    variables_file << "nodal volatile mass\tnvm\tkg m-3\t1\t0\t1.0E+10\tnode" << std::endl;
    variables_file << "thickness\th\tm\t1\t1.00E-25\t1.00E+00\telement" << std::endl;
    variables_file << "permeability\tk\tm2\t1\t1.00E-25\t1.00E+00\telement" << std::endl;
    variables_file << "previous permeability\tpk\tm2\t1\t1.00E-25\t1.00E+00\telement" <<std::endl;
    variables_file << "previous previous permeability\tppk\tm2\t1\t1.00E-25\t1.00E+00\telement" << std::endl;
    variables_file << "nodal permeability\tk\tm2\t1\t1.00E-25\t1.00E+00\tnode" << std::endl;
    variables_file << "total heat capacity\tCPT\tJ kg-1 oC-1\t1\t0.00E+00\t1.00E+20\telement"<< std::endl;
    variables_file << "temperature\tT\toC\t1\t0.00E+00\t1.00E+03\tnode" << std::endl;
    variables_file << "initial temperature\tTinit\toC\t1\t0.00E+00\t1.00E+03\tnode" << std::endl;
    variables_file << "temperature pre P\tT\toC\t1\t0.00E+00\t1.00E+03\tnode" << std::endl;
    variables_file << "previous temperature\tpT\toC\t1\t0.00E+00\t1.00E+03\tnode" <<std::endl;
    variables_file << "temperature element\tT\toC\t1\t0.00E+00\t1.00E+03\telement" <<std::endl;
    variables_file << "nodal heat capacity\tnCP\tJ kg-1 oC-1\t1\t0.00E+00\t1.00E+20\tnode" <<std::endl;
    variables_file << "thermal conductivity\tKT\tJ kg-1 oC-1\t1\t-1.00E+01\t1.00E+01\telement"<< std::endl;
    variables_file << "nodal heat flux bottom\tQHB\tJ kg-1 s-1\t1\t-1.00E+10\t1.00E+10\tnode" << std::endl;
    variables_file << "nodal total compressibility\tnCT\tPa-1\t1\t1.00E-17\t1.00E-03\tnode" << std::endl;
    variables_file << "fluid pressure\tPF\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "test fluid pressure\ttPF\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "fluid pressure element\tPFe\tPa\t1\t0.00E+00\t1.00E+10\telement" << std::endl;
    variables_file << "reference pressure\trPF\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "failure pressure\tfPF\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "liquid mass mobility\trKl\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "vapor mass mobility\trKv\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "liquid mass mobility density\trKl\tkg s-1 m-1\t1\t1.00E-25\t1.00E+03\tnode" << std::endl;
    variables_file << "vapor mass mobility density\trKv\tkg s-1 m-1\t1\t1.00E-25\t1.00E+03\tnode" << std::endl;
    variables_file << "nodal fluid volume source\tnQT\tkg kg-1 s-1\t1\t-1.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "nodal fluid volume source2\tnQT2\tkg kg-1 s-1\t1\t-1.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "previous nodal fluid volume source\tpnQT\tkg kg-1 s-1\t1\t-1.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "fluid source rate\tfsr\tkg s-1\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "porosity\tPHI\tX\t1\t0.00E+00\t1.00E+00\telement" << std::endl;
    variables_file << "pore volume\tPV\tm3\t1\t0.00E+00\t1.00e+20\tnode" << std::endl;
    variables_file << "pore volume before\tPVb\tm3\t1\t0.00E+00\t1.00e+20\tnode" << std::endl;
    variables_file << "pore volume element\tPVe\tm3\t1\t0.00E+00\t1.00e+20\telement" << std::endl;
    variables_file << "fluid density\tRF\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "mass\tm\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "previous mass\tPm\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "nodal fluid mass accumulated\tmMf\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "mass KCl\tmKCl\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "previous mass KCl\tmKCl\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "boundary flow mass\tbfm\tkg m-3\t1\t-1.00E+07\t1.00E+07\tnode" << std::endl;
    variables_file << "total bfm mass out\tbfmot\tkg\t1\t-1.00E+30\t1.00E+30\tmodel" << std::endl;
    variables_file << "enthalpy content liquid\tCHl\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "enthalpy content vapor\tCHv\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "previous temperature\tTp\toC\t1\t0.00E+00\t1.00E+03\tnode" << std::endl;
    variables_file << "nodal porosity\tnPHI\tX\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "previous nodal porosity\tnpPHI\tX\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "nodal heat capacity rock\tnCPR\tJ kg-1 oC-1\t1\t0.00E+00\t2.00E+03\tnode" << std::endl;
    variables_file << "nodal density rock\tnRHOR\tkg m-3\t1\t1.00E-05\t1.00E+04\tnode" << std::endl;
    variables_file << "previous nodal density rock\tpnRHOR\tkg m-3\t1\t1.00E-05\t1.00E+04\tnode" << std::endl;
    variables_file << "previous enthalpy content liquid\tPCHl\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "previous enthalpy content vapor\tPCHv\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "fluid source h\tfsh\tJ s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "fluid mass liquid\tMl\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "fluid mass vapor\tMv\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "previous fluid mass liquid\tPMl\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "previous fluid mass vapor\tPMv\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "previous fluid density\tPRF\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "total enthalpy\tHt\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "previous total enthalpy\tHtp\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "previous previous total enthalpy\tHtpp\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "enthalpy\tH\tJ m-3\t1\t-2.00E+20\t2.00E+20\tnode" << std::endl;
    variables_file << "previous enthalpy\tPH\tJ m-3\t1\t-2.00E+20\t2.00E+20\tnode" << std::endl;
    variables_file << "boundary flow enthalpy\tbfH\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "fluid enthalpy\tHF\tJ kg-1\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "enthalpy\tH\tJ m-3\t1\t-2.00E+20\t2.00E+20\tnode" << std::endl;
    variables_file << "saturation liquid\tSl\tX\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "saturation vapor\tSv\tX\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "density liquid\tRHOl\tkg m-3\t1\t1.00E+00\t2.00E+03\tnode" << std::endl;
    variables_file << "density vapor\tRHOv\tkg m-3\t1\t1.00E+00\t2.00E+03\tnode" << std::endl;
    variables_file << "viscosity liquid\tMUl\tPa s-1\t1\t1.00E-07\t1.00E-03\tnode" << std::endl;
    variables_file << "viscosity vapor\tMUv\tPa s-1\t1\t1.00E-07\t1.00E-03\tnode" << std::endl;
    variables_file << "fluid viscosity\tMUf\tPa s-1\t1\t1.00E-07\t1.00E-03\tnode" << std::endl;
    variables_file << "enthalpy liquid\tHl\tJ kg-1\t1\t1.00E+02\t1.00E+20\tnode" << std::endl;
    variables_file << "enthalpy vapor\tHv\tJ kg-1\t1\t1.00E+02\t1.00E+20\tnode" << std::endl;
    variables_file << "fluid heat capacity\tCPF\tJ kg-3 C-1\t1\t0.00E+00\t1.00E+04\tnode" << std::endl;
    variables_file << "compressibility\tB\tPa-1\t1\t0.00E+00\t1.00E-03\tnode" << std::endl;
    variables_file << "volumetric enthalpy liquid\tVHl\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "volumetric enthalpy vapor\tVHv\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "density liquid transport\tRHOlT\tkg m-3\t1\t1.00E+00\t2.00E+08\tnode" << std::endl;
    variables_file << "density vapor transport\tRHOvT\tkg m-3\t1\t1.00E+00\t2.00E+08\tnode" << std::endl;
    variables_file << "fluid mass\tMF\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "previous fluid mass\tpMF\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "volume factor\tOPF\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "fluid state\tFS\t-\t1\t-1.00E+00\t1.00E+01\tnode" << std::endl;
    variables_file << "bulk fluid density\tBRF\tkg m-3\t1\t0.00E+00\t2.00E+03\tnode" << std::endl;
    variables_file << "nodal compressibility rock\tBR\tPa-1\t1\t0.00E+00\t1.00E-03\tnode" << std::endl;
    variables_file << "liquid enthalpy mobility\trKl\tm2 s-1\t1\t1.00E-25\t1.00E+10\tnode" << std::endl;
    variables_file << "vapor enthalpy mobility\trKv\tm2 s-1\t1\t1.00E-25\t1.00E+10\tnode" << std::endl;
    variables_file << "lithostatic pressure\tPL\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "liquid enthalpy mobility density\trKld\tm2 s-1\t1\t1.00E-25\t1.00E+10\tnode" << std::endl;
    variables_file << "vapor enthalpy mobility density\trKvd\tm2 s-1\t1\t1.00E-25\t1.00E+10\tnode" << std::endl;
    variables_file << "density rock\tRHOR\tkg m-3\t1\t1.00E-05\t1.00E+04\telement" << std::endl;
    variables_file << "velocity\tV\tm s-1\t2\t-1.00E+10\t1.00E+10\telement" << std::endl;
    variables_file << "pore velocity\tPV\tm s-1\t2\t-1.00E+02\t1.00E+02\telement" << std::endl;
    variables_file << "nodal source liquid\tnQTl\tkg kg-1 s-1\t1\t-1.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "nodal source vapor\tnQTv\tkg kg-1 s-1\t1\t-1.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "previous fluid pressure\tPPF\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "previous previous fluid pressure\tPPPF\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "previous fluid state\tpFS\t-\t1\t-0.1\t1.00E+01\tnode" << std::endl;
    variables_file << "diffusivity\tDI\tX\t1\t0.00E+00\t1.00E+00\telement" << std::endl;
    variables_file << "velocity liquid\tVl\tm s-1\t2\t-1.00E+02\t1.00E+02\telement" << std::endl;
    variables_file << "velocity vapor\tVv\tm s-1\t2\t-1.00E+02\t1.00E+02\telement" << std::endl;
    variables_file << "pore velocity liquid\tPVl\tm s-1\t2\t-1.00E+02\t1.00E+02\telement" << std::endl;
    variables_file << "pore velocity vapor\tPVv\tm s-1\t2\t-1.00E+02\t1.00E+02\telement" << std::endl;
    variables_file << "nodal pore velocity\tpvn\tm s-1\t2\t-1.00E+02\t1.00E+02\tnode" << std::endl;
    variables_file << "relperm viscosity liquid\tMUl-1\ts Pa-1\t1\t1.00E-07\t1.00E+07\tnode" << std::endl;
    variables_file << "relperm viscosity vapor\tMUv-1\ts Pa-1\t1\t1.00E-07\t1.00E+07\tnode" << std::endl;
    variables_file << "KgradP\tKgradP\tm2 Pa m-1\t2\t1.00E-10\t1.00E+10\telement" << std::endl;
    variables_file << "gradP scaling\tgradP_sc\t-\t1\t1.00E-10\t1.00E+10\telement" << std::endl;
    variables_file << "courant liquid\tcfl_l\ts\t1\t1.00E-10\t1.00E+20\telement" << std::endl;
    variables_file << "courant vapor\tcfl_v\ts\t1\t1.00E-10\t1.00E+20\telement" << std::endl;
    variables_file << "upwind control liquid\tuc_l\t-\t1\t0\t30\telement" << std::endl;
    variables_file << "upwind control vapor\tuc_v\t-\t1\t0\t30\telement\t" << std::endl;
    variables_file << "salt fraction fluid\tXF\tWt. NaCl\t1\t0.00E+00\t1.00E+02\tnode" << std::endl;
    variables_file << "salt fraction liquid\tXl\tWt. NaCl\t1\t0.00E+00\t1.00E+02\tnode" << std::endl;
    variables_file << "salt fraction vapor\tXv\tWt. NaCl\t1\t0.00E+00\t1.00E+02\tnode" << std::endl;
    variables_file << "salt fraction halite\tXh\tWt. NaCl\t1\t0.00E+00\t1.00E+02\tnode" << std::endl;
    variables_file << "salt content liquid\tCXl\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "salt content vapor\tCXv\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "salt content halite\tCXv\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "salt content liquid KCl\tKCl_CXl\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "salt content vapor KCl\tKCl_CXv\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "salt content halite KCl\tKCl_CXh\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "salt content fluid\tCXf\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous salt content fluid\tpCXf\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous salt content liquid\tpCXl\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous salt content vapor\tpCXv\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous mass salt\tpms\tkg\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous salt content liquid KCl\tKCl_CXlp\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous salt content vapor KCl\tKCl_CXv\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous salt content halite KCl\tKCl_CXhp\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "boundary flow salt\tbfss\tkg\t1\t-1.00E+09\t1.00E+09\tnode" << std::endl;
    variables_file << "salinity\tS\tWt. NaCl\t1\t0.00E+00\t1.00E+02\tnode" << std::endl;
    variables_file << "salinity element\tSe\tWt. NaCl\t1\t0.00E+00\t1.00E+02\telement" << std::endl;
    variables_file << "previous salinity\tSp\tWt. NaCl\t1\t0.00E+00\t1.00E+02\tnode" << std::endl;
    variables_file << "liquid salt mobility\trKl\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "vapor salt mobility\trKv\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "liquid salt mobility KCl\trKl_KCl\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "vapor salt mobility KCl\trKv_KCl\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "fluid source wt\tfsw\tkg s-1\t1\t0.0\t1.0e10\tnode" << std::endl;
    variables_file << "saturation halite\tSh\tX\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "saturation halite element\tShe\tX\t1\t0.00E+00\t1.00E+00\telement" << std::endl;
    variables_file << "density halite\tRHOh\tkg m-3\t1\t1.00E+00\t2.00E+03\tnode" << std::endl;
    variables_file << "solid mass halite\tMh\tkg m-3\t1\t0.00E+00\t1.00E+04\tnode" << std::endl;
    variables_file << "enthalpy halite\tHh\tJ kg-1\t1\t1.00E+02\t1.00E+20\tnode" << std::endl;
    variables_file << "volumetric enthalpy halite\tVHh\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "enthalpy content halite\tCHv\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "volumetric salinity vapor\tVXv\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "volumetric salinity liquid\tVXl\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "volumetric salinity halite\tVXh\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "fracturing reference\tcf\t-\t1\t0.00E+00\t3.00E+01\telement" << std::endl;
    // variables_file << "fracturation reference\tcf\t-\t1\t0.00E+00\t3.00E+01\telement" << std::endl;
    // variables_file << "time factor\ttf\t-\t1\t0.00E+00\t1.01E+00\tnode" << std::endl;
    variables_file << "fracturing events\tfe\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
    variables_file << "permeability increase\tki\tm2\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
    variables_file << "permeability log increase\tkli\tm2\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
    variables_file << "permeability decrease\tkd\tm2\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
    variables_file << "permeability log decrease\tkld\tm2\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
    variables_file << "energy flux per second\tEfps\tJ m-3 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "energy flux integral\tEfI\tJ m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "energy flux\tEf\tJ m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "fluid flux per second\tffps\tkg m-3 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "fluid flux integral\tffiI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "fluid flux with direction integral\tffwdi\tm s-1\t2\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "fluid flux with direction\tffwd\tm s-1\t2\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "fluid flux\tffi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "salt flux per second\tSfps\tkg m-3 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "salt flux integral\tSfI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "salt flux\tSfI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "previous fluid flux integral\tpffI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "previous salt flux integral\tpSfI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "previous energy flux integral\tpEfI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "previous compressibility\tBP\tPa-1\t1\t0.00E+00\t1.00E-03\tnode" << std::endl;
    variables_file << "reference compressibility\tBP\tPa-1\t1\t0.00E+00\t1.00E-03\tnode" << std::endl;

    variables_file << "dangerous phase change\tDPC\t-\t1\t0.00E+00\t1.00E+03\tnode" << std::endl;
    variables_file << "after phasechange counter\tAPC\t-\t1\t0.00E+00\t2.00E+01\tnode" << std::endl;
    variables_file << "reference specific enthalpy top\trset\tJ kg-1\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "dt\tdt\ts\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "bulk volume\tVl\tm3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "bulk volume element\tVle\tm3\t1\t0.00E+00\t1.00E+10\telement" << std::endl;
    variables_file << "previous pore volume\tVl\tm3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "previous bulk volume\tVl\tm3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "initial bulk volume\tVl0\tm3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "initial pore volume\tpVl\tm3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    variables_file << "region ID\trID\t-\t1\t0.00E+00\t1.00E+01\telement" << std::endl;
    variables_file << "with salinity\twsal\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with passive tracer\twpt\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with magma chamber\twmc\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with gold\twAu\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with quartz\twQtz\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with flux visitor\twqv\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with lithium\twLi\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with magma model\twmm\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with point source\twpt\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;

    variables_file << "with temperature halo\tThalo\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with lithostatic pressure\tlithPflag\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "open top\toTop\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "barycenter depth\tm\tm\t1\t0.00E+00\t1.00E+06\telement" << std::endl;

    variables_file << "with anisotropic permeability\tanik\t-\t1\t0\t1.0E+00\tmodel" << std::endl;
    variables_file << "permeability anisotropy factor\tanikfact\t-\t1\t1.0E-10\t1.0E+10\telement" << std::endl;
    variables_file << "vertical permeability\tkv\tm2\t1\t1.00E-25\t1.00E+00\telement" << std::endl;
    variables_file << "horizontal permeability\tkh\tm2\t1\t1.00E-25\t1.00E+00\telement" << std::endl;
    variables_file << "log vertical permeability increase\tlogkinc\tlogm2\t1\t-1.00E+01\t2.00E+01\telement" << std::endl;
    variables_file << "log horizontal permeability increase\tlogkinc\tlogm2\t1\t-1.00E+01\t2.00E+01\telement" << std::endl;

    variables_file << "depth dependent permeability\tdk\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "with reference depth\twrefd\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "reference depth\trefd\t-\t1\t-1.00E+10\t1.00E+10\tmodel" << std::endl;
    variables_file <<"temperature dependent permeability\tTk\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "mineral dependent permeability\tmk\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "pore fluid factor dependent permeability\tPffk\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "hydrofracturing\thf\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "immediate fracture closure\tfclosure\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "change BDT temperature\tcBDTT\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "near critically pressured\tcritP\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "BDT start temperature\tBDTTs\toC\t1\t0.00E+00\t1.00E+04\tmodel" << std::endl;
    variables_file << "BDT ductile temperature\tBDTTd\toC\t1\t0.00E+00\t1.00E+04\tmodel" << std::endl;
    variables_file << "BDT end temperature\tBDTTe\toC\t1\t0.00E+00\t1.00E+04\tmodel" << std::endl;
    variables_file << "BDT start log permeability\tBDTks\tm2\t1\t-25.00E+00\t-1.00E+00\tmodel" << std::endl;
    variables_file << "BDT start permeability element\tBDTks\tm2\t1\t-25.00E+00\t-1.00E+00\telement" << std::endl;
    variables_file << "BDT ductile log permeability\tBDTks\tm2\t1\t-25.00E+00\t-1.00E+00\tmodel" << std::endl;
    variables_file << "BDT end log permeability\tBDTkd\tm2\t1\t-25.00E+00\t-1.00E+00\tmodel" << std::endl;
    variables_file <<"maximum log permeability model\tkmax\tm2\t1\t-25.00E+00\t-1.00E+00\tmodel" << std::endl;
    variables_file <<"initial wtpercent volatile in magma\tiwtpVolMag\tperc\t1\t0.00E+00\t1.00E+02\tmodel" <<std::endl;
    // variables_file << "magma volatile saturation\tvolsatmag\t-\t1\t0.00E+00\t1.00E+00\tmodel"<< std::endl;

    variables_file << "crystallization curve exponent\tccexp\t-\t1\t1.00E+00\t1.00E+02\tmodel"<< std::endl;

    variables_file << "liquidus temperature\tliqT\toC\t1\t0.00E+00\t1.00E+03\tnode" <<std::endl;
    // variables_file << "previous liquidus temperature\tpliqT\toC\t1\t0.00E+00\t1.00E+03\tnode"<< std::endl; // JK used??
    variables_file << "solidus temperature\tsoldT\toC\t1\t0.00E+00\t1.00E+03\tnode" <<std::endl;
    // variables_file << "previous solidus temperature\tpsoldT\toC\t1\t0.00E+00\t1.00E+03\tnode"<< std::endl; // JK used??

    variables_file << "model time\ttime\ts \t1\t0.00E+00\t1.00E+20\tmodel" << std::endl;

    variables_file << "pore volume change factor\tdelta_pV\t-\t1\t0.\t1.E+9\tnode" << std::endl;

    variables_file << "threshold pressure\tPthresh\tPa\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "failure overpressure\tfOP\tPa\t1\t0.\t1.E+9\tnode" << std::endl;
    variables_file << "threshold overpressure\tOPthresh\tPa\t1\t0.\t1.E+9\tnode" << std::endl;
    variables_file << "fracturable flag\tfrac_flag\t-\t1\t0.\t1.E+9\telement" << std::endl;

    variables_file << "nodal depth\tnd\tm\t1\t0.\t1.E+9\tnode" << std::endl;
    variables_file << "volume\tV\tm3\t1\t0.00E+00\t1.00e+20\tnode" << std::endl;
    variables_file << "gravity flag\tg\tx\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "salinity top\tStop\tWt. NaCl\t1\t0.00E+00\t1.00E+02\tmodel" <<std::endl;
    variables_file << "cfl scaling\tcfl\t-\t1\t0.00E+00\t1.00E+04\tmodel" << std::endl;
    variables_file << "heat capacity melt\tCPrm\tJ kg-1 oC-1\t1\t0.00E+00\t1.00E+20\tmodel" << std::endl;
    variables_file << "with custom latent heat\twclh\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "latent heat of fusion\tlhf\tJ kg-1\t1\t0.00E+00\t1.00E+20\tmodel" << std::endl;
    variables_file << "permeability tensor\tkT\tm2\t3\t1.00E-25\t1.00E+00\telement" << std::endl;
    variables_file << "differential stress\tdiff_stress\tPa\t1\t0.00E+00\t1.00E+10\tmodel" << std::endl;
    variables_file << "cohesion\tC\t-\t1\t0.00E+00\t1.00E+10\tmodel" << std::endl;
    variables_file << "initial T diffusion time\tTdiff_t0\ts\t1\t0.00E+00\t1.00E+20\tmodel" << std::endl;
    variables_file << "initial intrusion temperature\tTmagma\toC\t1\t0.00E+00\t1.00E+04\tmodel" << std::endl;
    variables_file << "no standard output\tstdcout\t-\t1\t0.\t1\tmodel" << std::endl;
    variables_file << "coordinate x\tx\tm\t1\t-1.00E+05\t1.00E+05\tnode" << std::endl;
    variables_file << "coordinate y\ty\tm\t1\t-1.00E+05\t1.00E+05\tnode" << std::endl;
    variables_file << "coordinate z\tz\tm\t1\t-1.00E+05\t1.00E+05\tnode" << std::endl;
    variables_file << "save model increment\tsave_t\t-\t1\t0.00E+00\t1.00E+20\tmodel" << std::endl;
    variables_file << "fluid flux in\tffi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;                           // JK used??
    variables_file << "fluid flux out\tffo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;                          // JK used??
    variables_file << "total fluid flux in\ttffi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;                    // JK used??
    variables_file << "total fluid flux out\ttffo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;                   // JK used??

    variables_file << "previous total fluid flux in\ttffi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;           // JK used??
    variables_file << "previous total fluid flux out\ttffo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;          // JK used??
    variables_file << "conductivity\tK\tm2 Pa-1 s-1\t1\t0.00E+00\t1.00E+00\telement integration point" << std::endl;
    variables_file << "mass conductivity\trK\ts\t1\t0.00E+00\t1.00E+00\telement integration point" << std::endl;        // JK used??
    variables_file << "temperature dependent differential stress flag\tTdep_diffP\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;

    // variables_file << "unit normal vector\tun\t-\t2\t-1.00E+00\t1.00E+00\tnode" << std::endl;

    variables_file << "permeability ID\tkID\t-\t1\t0.00E+00\t1.00E+01\telement" << std::endl;
    variables_file << "log permeability increase\tlogkinc\tlogm2\t1\t-1.00E+01\t2.00E+01\telement" << std::endl;

    variables_file << "crystallization curve\tXcurvel\t-\t1\t0.00E+00\t2.00E+00\tmodel" << std::endl;

    // Porosity-permeability coupling
    variables_file << "with depth dependent porosity permeability coupling\twPhiKCouplingt\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "lower porosity limit\tphiMint\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "upper porosity limit\tphiMaxt\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;

    //! Quartz
    variables_file << "quartz solubility\tQS\tmol kg-1\t1\t-1.00E+20\t1.00E+20\tnode" << std::endl;
    variables_file << "quartz solid\tQtzm\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "initial quartz solid\tinQtzm\tkg m-3\t1\t0.00E+00\t1.00E+21\tnode" << std::endl;
    variables_file << "vein quartz solid\tvQtzm\tkg m-3\t1\t0.00E+00\t1.00E+22\tnode" << std::endl;
    variables_file << "minimum quartz solid\tminQtzm\tkg m-3\t1\t0.00E+00\t1.00E+23\tnode" << std::endl;
    variables_file << "delta quartz precipitated\tdQtzp\tkg m-3\t1\t-1.00E+20\t1.00E+24\tnode" << std::endl;



    //! Alteration Visitor
    variables_file << "with alteration\twalt\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "min temperature alteration window\tTminAlt\tdegC\t1\t0.00E+00\t4.00E+02\tmodel" << std::endl;
    variables_file << "max temperature alteration window\tTmaxAlt\tdegC\t1\t0.00E+00\t4.00E+02\tmodel" << std::endl;
    variables_file << "characteristic alteration time\taltTime\tyears\t1\t0.00E+00\t1.00E+04\tmodel" << std::endl;
    variables_file << "characteristic fluid flux\tcq\tkg m-2 s-1\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;

    variables_file << "initial permeability\tk_init\tm2\t1\t1.00E-25\t1.00E+00\telement" << std::endl;
    variables_file << "initial vertical permeability\tkv_init\tm2\t1\t1.00E-25\t1.00E+00\telement" << std::endl;
    variables_file << "initial horizontal permeability\tkh_init\tm2\t1\t1.00E-25\t1.00E+00\telement" << std::endl;
    variables_file << "initial porosity\tphi_init\t-\t1\t0.00E+00\t1.00E+00\telement" << std::endl;
    variables_file << "alteration index\tAidx\t-\t1\t0.00E+00\t1.00E+00\telement" << std::endl;
    variables_file << "density liquid element\tRHOle\tkg m-3\t1\t1.00E+00\t2.00E+03\telement" << std::endl;
    variables_file << "density vapor element\tRHOve\tkg m-3\t1\t1.00E+00\t2.00E+03\telement" << std::endl;
    variables_file << "saturation liquid element\tSle\tX\t1\t0.00E+00\t1.00E+00\telement" << std::endl;
    variables_file << "saturation vapor element\tSve\tX\t1\t0.00E+00\t1.00E+00\telement" << std::endl;

    variables_file << "condensation tracker\tcondt\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;


    //! Subcycling
    variables_file << "stash pressure\tsPF\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;

    //! Flux visitor
    if (with_flux_visitor)
    {
        variables_file << "LDE ID\teID\t-\t1\t0.00E+00\t1.00E+06\telement" << std::endl;
        variables_file << "unit normal vector LDE\tun\t-\t2\t-1.00E+00\t1.00E+00\telement" << std::endl;

        variables_file << "average velocity liquid LDE\taVl\t-\t2\t-1.00E+20\t1.00E+20\telement" << std::endl;
        variables_file << "average velocity vapor LDE\taVv\t-\t2\t-1.00E+20\t1.00E+20\telement" << std::endl;
        variables_file << "average velocity fluid LDE\taVf\t-\t2\t-1.00E+20\t1.00E+20\telement" << std::endl;
        variables_file << "fluid flux liquid element\tffl\t-\t2\t-1.00E+20\t1.00E+20\telement" << std::endl;
        variables_file << "fluid flux vapor element\tffv\t-\t2\t-1.00E+20\t1.00E+20\telement" << std::endl;
        variables_file << "fluid flux fluid element\tfff\t-\t2\t-1.00E+20\t1.00E+20\telement" << std::endl;
        variables_file << "energy flux liquid element\tefl\t-\t2\t-1.00E+20\t1.00E+20\telement" << std::endl;
        variables_file << "energy flux vapor element\tefv\t-\t2\t-1.00E+20\t1.00E+20\telement" << std::endl;
        variables_file << "energy flux fluid element\teff\t-\t2\t-1.00E+20\t1.00E+20\telement" << std::endl;

        variables_file << "integrated fluid flux up liquid element\tifflud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated fluid flux up vapor element\tiffvud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated fluid flux up fluid element\tifffud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated fluid flux down liquid element\tiffld\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated fluid flux down vapor element\tiffvd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated fluid flux down fluid element\tifffd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;

        variables_file << "integrated energy flux up liquid element\tieflud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated energy flux up vapor element\tiefvud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated energy flux up fluid element\tieffud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated energy flux down liquid element\tiefld\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated energy flux down vapor element\tiefvd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated energy flux down fluid element\tieffd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;

        variables_file << "integrated fluid flux out liquid element\tifflud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated fluid flux out vapor element\tiffvud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated fluid flux out fluid element\tifffud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated fluid flux in liquid element\tiffld\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated fluid flux in vapor element\tiffvd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated fluid flux in fluid element\tifffd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;

        variables_file << "integrated energy flux out liquid element\tieflud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated energy flux out vapor element\tiefvud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated energy flux out fluid element\tieffud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated energy flux in liquid element\tiefld\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated energy flux in vapor element\tiefvd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated energy flux in fluid element\tieffd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;

        variables_file << "integrated li flux up liquid element\tiliflud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated li flux up vapor element\tilifvud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated li flux up fluid element\tiliffud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated li flux down liquid element\tilifld\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated li flux down vapor element\tilifvd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated li flux down fluid element\tiliffd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;

        variables_file << "li concentration liquid element\tliconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\telement" << std::endl;
        variables_file << "li concentration vapor element\tliconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\telement" << std::endl;
        variables_file << "li concentration fluid element\tliconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\telement" << std::endl;

        variables_file << "integrated magmatic flux up fluid element\tisnffud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated magmatic flux down fluid element\tisnffd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated magmatic flux in fluid element\tisnffud\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated magmatic flux out fluid element\tisnffd\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;

        variables_file << "integrated salt flux up fluid element\tisffud\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "integrated salt flux down fluid element\tisffd\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated salt flux in fluid element\tisnffud\t-\t1\t-1.00E+20\t0.00E+00\telement" << std::endl;
        variables_file << "integrated salt flux out fluid element\tisnffd\t-\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
    }


    //! Gold Visitor
    if (with_gold)
    {
        variables_file << "split region SO2 source\tqSO2\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
        variables_file << "split region AuHS0 source\tqAuHS0\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;

        variables_file << "average crystallinity\tavx\t-\t1\t0.00E+00\t1.00E+02\tmodel" << std::endl;
        variables_file << "delta average crystallinity\tdavx\t-\t1\t-1.00E+02\t1.00E+02\tmodel" << std::endl;

        variables_file << "relative gold saturation\trgs\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
        variables_file << "total input gold mass\tmauinp\tkg\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        variables_file << "total input sulfur mass\tmsinp\tkg\t1\t0.00E+00\t1.00E+24\tmodel\t" << std::endl;
        variables_file << "redox model\trxm\t-\t1\t0.00E+00\t1.00E+02\tmodel" << std::endl;
        variables_file << "ph model\tphm\t-\t1\t0.00E+00\t1.00E+02\tmodel" << std::endl;
        variables_file << "gold source model\tgsm\t-\t1\t0.00E+00\t1.00E+02\tmodel" << std::endl;
        variables_file << "with low permeability\twlk\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
        variables_file << "with initial solid gold\twisg\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;

        variables_file << "gold source region ID\tgsrID\t-\t1\t0.00E+00\t1.00E+01\tnode" << std::endl;
        variables_file << "ph buffered difference to neutral\tdph\t-\t1\t-1.00E+20\t1.00E+24\tmodel" << std::endl;
        variables_file << "magmatic so2 concentration\tso2mag\tmol kg-1\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        // variables_file << "magmatic co2 concentration\tco2mag\tmol kg-1\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        variables_file << "magmatic gold concentration\taumag\tmol kg-1\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        variables_file << "h2 activity\tah2\tmol kg-1\t1\t0.00E+00\t1.00E+24\tmodel\t" << std::endl;
        variables_file << "so2 partition coefficient fluid melt\tDso2fm\t-\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        // variables_file << "co2 partition coefficient fluid melt\tDco2fm\t-\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        variables_file << "gold partition coefficient fluid melt\tDaufm\t-\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        variables_file << "so2 partition coefficient crystals melt\tDso2cm\t-\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        // variables_file << "co2 partition coefficient crystals melt\tDco2cm\t-\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        variables_file << "gold partition coefficient crystals melt\tDaucm\t-\t1\t0.00E+00\t1.00E+24\tmodel" << std::endl;
        variables_file << "gold slope partition coefficient vapor liquid\tDgradauvl\t-\t1\t-1.00E+24\t1.00E+24\tmodel" << std::endl;
        variables_file << "with inhibited gold redissolution\twaudissinh\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
        variables_file << "with gold vapor liquid partitioning\twauvlpart\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
        variables_file << "with ph buffering\twphbuff\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;

        // variables_file << "h2 content liquid\th2Cl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "h2 content vapor\th2Cv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "h2 content fluid\th2Cf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 content fluid\tco2Cf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 content liquid\tco2Cl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 content vapor\tco2Cv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 content fluid\tso2Cf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 content liquid\tso2Cl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 content vapor\tso2Cv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 content fluid\tauhs0Cf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 content liquid\tauhs0Cl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 content vapor\tauhs0Cv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 content fluid\tsio2Cf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 content liquid\tsio2Cl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 content vapor\tsio2Cv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid h2 mobility\th2lm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor h2 mobility\th2vm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid so2 mobility\tso2lm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor so2 mobility\tso2vm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "liquid co2 mobility\tco2lm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "vapor co2 mobility\tco2vm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid auhs0 mobility\tauhs0lm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor auhs0 mobility\tauhs0vm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid sio2 mobility\tsio2lm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor sio2 mobility\tsio2vm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;

        variables_file << "delta quartz precipitated element\tdQtzpe\tkg m-3\t1\t-1.00E+20\t1.00E+25\telement" << std::endl;
        variables_file << "delta gold precipitated element\tdAupe\tkg m-3\t1\t-1.00E+20\t1.00E+36\telement\t" << std::endl;
        variables_file << "chalcedony solid\tChalcm\tkg m-3\t1\t0.00E+00\t1.00E+26\tnode" << std::endl;
        variables_file << "initial chalcedony solid\tinChalcm\tkg m-3\t1\t0.00E+00\t1.00E+27\tnode" << std::endl;
        variables_file << "vein chalcedony solid\tvChalcm\tkg m-3\t1\t0.00E+00\t1.00E+28\tnode" << std::endl;
        variables_file << "minimum chalcedony solid\tminChalcm\tkg m-3\t1\t0.00E+00\t1.00E+29\tnode" << std::endl;
        variables_file << "delta chalcedony precipitated\tdeltaChalcprec\tkg m-3\t1\t-1.00E+20\t1.00E+31\tnode\t" << std::endl;
        variables_file << "amorphous silica solid\tsiamm\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "initial amorphous silica solid\tinsiamm\tkg m-3\t1\t0.00E+00\t1.00E+21\tnode" << std::endl;
        variables_file << "delta amorphous silica precipitated\tdsiamp\tkg m-3\t1\t-1.00E+20\t1.00E+24\tnode" << std::endl;
        variables_file << "delta amorphous silica precipitated element\tdsiampe\tkg m-3\t1\t-1.00E+20\t1.00E+25\telement" << std::endl;
        variables_file << "gold solid\tAum\tkg m-3\t1\t0.00E+00\t1.00E+31\tnode" << std::endl;
        variables_file << "gold solid by boiling\tAumb\tkg m-3\t1\t-1.00E+20\t1.00E+32\tnode" << std::endl;
        variables_file << "bulk gold grade\taugradeb\tmg Au/kg rock\t1\t0.00E+00\t1.00E+33\tnode" << std::endl;
        variables_file << "vein gold grade\taugradev\tmg Au/kg rock\t1\t0.00E+00\t1.00E+34\tnode" << std::endl;
        variables_file << "delta gold precipitated\tdAup\tkg m-3\t1\t-1.00E+20\t1.00E+35\tnode" << std::endl;
        variables_file << "gold saturation index\tauSI\t-\t1\t0.00E+00\t1.00E+37\tnode" << std::endl;
        variables_file << "total solid gold in system\ttotAusys\tkg\t1\t0.00E+00\t1.00E+38\tmodel" << std::endl;
        variables_file << "gold content fluid\tAuCf\tkg m-3\t1\t0.00E+00\t1.00E+39\tnode" << std::endl;
        variables_file << "gold content liquid\tAuCl\tkg m-3\t1\t0.00E+00\t1.00E+40\tnode" << std::endl;
        variables_file << "gold content vapor\tAuCv\tkg m-3\t1\t0.00E+00\t1.00E+41\tnode" << std::endl;
        variables_file << "gold content melt\tAuCm\tkg m-3\t1\t0.00E+00\t1.00E+42\tnode" << std::endl;
        variables_file << "liquid gold mobility\tAulm\tm2 s-1\t1\t0.00E+00\t1.00E+44\tnode" << std::endl;
        variables_file << "vapor gold mobility\tAuvm\tm2 s-1\t1\t0.00E+00\t1.00E+45\tnode" << std::endl;
        variables_file << "gold fraction fluid\tAufracf\t-\t1\t0.00E+00\t1.00E+46\tnode" << std::endl;
        variables_file << "gold fraction liquid\tAufracl\t-\t1\t0.00E+00\t1.00E+47\tnode" << std::endl;
        variables_file << "gold fraction vapor\tAufracv\t-\t1\t0.00E+00\t1.00E+48\tnode" << std::endl;
        variables_file << "gold concentration liquid\tAuconcl\tkg kg-1\t1\t0.00E+00\t1.00E+49\tnode" << std::endl;
        variables_file << "gold concentration vapor\tAuconcv\tkg kg-1\t1\t0.00E+00\t1.00E+49\tnode" << std::endl;
        variables_file << "gold concentration fluid\tAuconcf\tkg kg-1\t1\t0.00E+00\t11.00E+30\tnode" << std::endl;
        variables_file << "h2s content fluid\th2sCf\tkg m-3\t1\t0.00E+00\t12.00E+30\tnode" << std::endl;
        variables_file << "h2s content liquid\th2sCl\tkg m-3\t1\t0.00E+00\t13.00E+30\tnode" << std::endl;
        variables_file << "h2s content vapor\th2sCv\tkg m-3\t1\t0.00E+00\t14.00E+30\tnode" << std::endl;
        variables_file << "h2s content melt\th2sCm\tkg m-3\t1\t0.00E+00\t15.00E+30\tnode" << std::endl;
        variables_file << "previous h2s content fluid\th2sCfp\tkg m-3\t1\t0.00E+00\t16.00E+30\tnode" << std::endl;
        variables_file << "previous h2s content liquid\th2sClp\tkg m-3\t1\t0.00E+00\t17.00E+30\tnode" << std::endl;
        variables_file << "previous h2s content vapor\th2sCvp\tkg m-3\t1\t0.00E+00\t18.00E+30\tnode" << std::endl;
        variables_file << "liquid h2s mobility\th2slm\tm2 s-1\t1\t0.00E+00\t19.00E+30\tnode" << std::endl;
        variables_file << "vapor h2s mobility\th2svm\tm2 s-1\t1\t0.00E+00\t20.00E+30\tnode" << std::endl;
        variables_file << "h2s fraction fluid\th2sfracf\t-\t1\t0.00E+00\t21.00E+00\tnode" << std::endl;
        variables_file << "h2s fraction liquid\th2sfracl\t-\t1\t0.00E+00\t22.00E+00\tnode" << std::endl;
        variables_file << "h2s fraction vapor\th2sfracv\t-\t1\t0.00E+00\t23.00E+00\tnode" << std::endl;
        variables_file << "h2s concentration liquid\th2sconcl\tmol kg-1\t1\t0.00E+00\t24.00E+30\tnode" << std::endl;
        variables_file << "h2s concentration vapor\th2sconcv\tmol kg-1\t1\t0.00E+00\t25.00E+30\tnode" << std::endl;
        variables_file << "h2s concentration fluid\th2sconcf\tmol kg-1\t1\t0.00E+00\t26.00E+30\tnode" << std::endl;

        // variables_file << "h2 content melt\th2Cm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "h2 content crystals\th2Cc\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous h2 content fluid\th2Cfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous h2 content liquid\th2Clp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous h2 content vapor\th2Cvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "h2 fraction fluid\th2fracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "h2 fraction liquid\th2fracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "h2 fraction vapor\th2fracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "h2 concentration liquid\th2concl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "h2 concentration vapor\th2concv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "h2 concentration fluid\th2concf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "h2 partial pressure\th2pP\tbar\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;

        // variables_file << "co2 content melt\tco2Cm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 content crystals\tco2Cc\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous co2 content fluid\tco2Cfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous co2 content liquid\tco2Clp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous co2 content vapor\tco2Cvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 fraction fluid\tco2fracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "co2 fraction liquid\tco2fracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "co2 fraction vapor\tco2fracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "co2 concentration liquid\tco2concl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 concentration vapor\tco2concv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 concentration fluid\tco2concf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 concentration melt\tco2concm\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 concentration crystals\tco2concc\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 mass liquid\tco2ml\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 mass vapor\tco2mv\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 mass fluid\tco2mf\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 mass melt\tco2mm\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 mass crystals\tco2mc\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co2 mass total\tco2mtot\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "hco3m content fluid\thco3mCf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "hco3m content liquid\thco3mCl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "hco3m content vapor\thco3mCv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "hco3m content melt\thco3mCm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous hco3m content fluid\thco3mCfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous hco3m content liquid\thco3mClp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous hco3m content vapor\thco3mCvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "liquid hco3m mobility\thco3mlm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "vapor hco3m mobility\thco3mvm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "hco3m fraction fluid\thco3mfracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "hco3m fraction liquid\thco3mfracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "hco3m fraction vapor\thco3mfracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "hco3m concentration liquid\thco3mconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "hco3m concentration vapor\thco3mconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "hco3m concentration fluid\thco3mconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co3mm content fluid\tco3mmCf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co3mm content liquid\tco3mmCl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co3mm content vapor\tco3mmCv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co3mm content melt\tco3mmCm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous co3mm content fluid\tco3mmCfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous co3mm content liquid\tco3mmClp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "previous co3mm content vapor\tco3mmCvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "liquid co3mm mobility\tco3mmlm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "vapor co3mm mobility\tco3mmvm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co3mm fraction fluid\tco3mmfracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "co3mm fraction liquid\tco3mmfracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "co3mm fraction vapor\tco3mmfracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        // variables_file << "co3mm concentration liquid\tco3mmconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co3mm concentration vapor\tco3mmconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        // variables_file << "co3mm concentration fluid\tco3mmconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs2m content fluid\tauhs2mCf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs2m content liquid\tauhs2mCl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs2m content vapor\tauhs2mCv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs2m content melt\tauhs2mCm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous auhs2m content fluid\tauhs2mCfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous auhs2m content liquid\tauhs2mClp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous auhs2m content vapor\tauhs2mCvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid auhs2m mobility\tauhs2mlm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor auhs2m mobility\tauhs2mvm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs2m fraction fluid\tauhs2mfracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "auhs2m fraction liquid\tauhs2mfracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "auhs2m fraction vapor\tauhs2mfracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "auhs2m concentration liquid\tauhs2mconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs2m concentration vapor\tauhs2mconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs2m concentration fluid\tauhs2mconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 content melt\tsio2Cm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 content crystals\tsio2Cc\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous sio2 content fluid\tsio2Cfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous sio2 content liquid\tsio2Clp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous sio2 content vapor\tsio2Cvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 fraction fluid\tsio2fracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "sio2 fraction liquid\tsio2fracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "sio2 fraction vapor\tsio2fracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "sio2 concentration liquid\tsio2concl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 concentration vapor\tsio2concv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 concentration fluid\tsio2concf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 concentration melt\tsio2concm\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 concentration crystals\tsio2concc\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 mass liquid\tsio2ml\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 mass vapor\tsio2mv\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 mass fluid\tsio2mf\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 mass melt\tsio2mm\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 mass crystals\tsio2mc\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "sio2 mass total\tsio2mtot\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;

        variables_file << "hsm concentration liquid\thsmconcl\tmol kg-1\t1\t0.00E+00\t7.05E+30\tnode" << std::endl;
        variables_file << "hsm concentration vapor\thsmconcv\tmol kg-1\t1\t0.00E+00\t7.06E+30\tnode" << std::endl;
        variables_file << "hsm concentration fluid\thsmconcf\tmol kg-1\t1\t0.00E+00\t7.07E+30\tnode" << std::endl;

        variables_file << "auhs0 content melt\tauhs0Cm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 content crystals\tauhs0Cc\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous auhs0 content fluid\tauhs0Cfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous auhs0 content liquid\tauhs0Clp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous auhs0 content vapor\tauhs0Cvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 fraction fluid\tauhs0fracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "auhs0 fraction liquid\tauhs0fracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "auhs0 fraction vapor\tauhs0fracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "auhs0 concentration liquid\tauhs0concl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 concentration vapor\tauhs0concv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 concentration fluid\tauhs0concf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 concentration melt\tauhs0concm\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 concentration crystals\tauhs0concc\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 mass liquid\tauhs0ml\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 mass vapor\tauhs0mv\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 mass fluid\tauhs0mf\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 mass melt\tauhs0mm\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 mass crystals\tauhs0mc\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auhs0 mass total\tauhs0mtot\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auoh0 content fluid\tauoh0Cf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auoh0 content liquid\tauoh0Cl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auoh0 content vapor\tauoh0Cv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auoh0 content melt\tauoh0Cm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous auoh0 content fluid\tauoh0Cfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous auoh0 content liquid\tauoh0Clp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous auoh0 content vapor\tauoh0Cvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid auoh0 mobility\tauoh0lm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor auoh0 mobility\tauoh0vm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auoh0 fraction fluid\tauoh0fracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "auoh0 fraction liquid\tauoh0fracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "auoh0 fraction vapor\tauoh0fracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "auoh0 concentration liquid\tauoh0concl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auoh0 concentration vapor\tauoh0concv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "auoh0 concentration fluid\tauoh0concf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "aucl2m content fluid\taucl2mCf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "aucl2m content liquid\taucl2mCl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "aucl2m content vapor\taucl2mCv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "aucl2m content melt\taucl2mCm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous aucl2m content fluid\taucl2mCfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous aucl2m content liquid\taucl2mClp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous aucl2m content vapor\taucl2mCvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid aucl2m mobility\taucl2mlm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor aucl2m mobility\taucl2mvm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "aucl2m fraction fluid\taucl2mfracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "aucl2m fraction liquid\taucl2mfracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "aucl2m fraction vapor\taucl2mfracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "aucl2m concentration liquid\taucl2mconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "aucl2m concentration vapor\taucl2mconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "aucl2m concentration fluid\taucl2mconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hsm content fluid\thsmCf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hsm content liquid\thsmCl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hsm content vapor\thsmCv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hsm content melt\thsmCm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous hsm content fluid\thsmCfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous hsm content liquid\thsmClp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous hsm content vapor\thsmCvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid hsm mobility\thsmlm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor hsm mobility\thsmvm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hsm fraction fluid\thsmfracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "hsm fraction liquid\thsmfracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "hsm fraction vapor\thsmfracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "hsm concentration liquid\thsmconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hsm concentration vapor\thsmconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hsm concentration fluid\thsmconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "smm content fluid\tsmmCf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "smm content liquid\tsmmCl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "smm content vapor\tsmmCv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "smm content melt\tsmmCm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous smm content fluid\tsmmCfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous smm content liquid\tsmmClp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous smm content vapor\tsmmCvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid smm mobility\tsmmlm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor smm mobility\tsmmvm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "smm fraction fluid\tsmmfracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "smm fraction liquid\tsmmfracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "smm fraction vapor\tsmmfracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "smm concentration liquid\tsmmconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "smm concentration vapor\tsmmconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "smm concentration fluid\tsmmconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;

        variables_file << "so2 content melt\tso2Cm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 content crystals\tso2Cc\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous so2 content fluid\tso2Cfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous so2 content liquid\tso2Clp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous so2 content vapor\tso2Cvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 fraction fluid\tso2fracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "so2 fraction liquid\tso2fracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "so2 fraction vapor\tso2fracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "so2 concentration liquid\tso2concl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 concentration vapor\tso2concv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 concentration fluid\tso2concf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 concentration melt\tso2concm\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 concentration crystals\tso2concc\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 mass liquid\tso2ml\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 mass vapor\tso2mv\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 mass fluid\tso2mf\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 mass melt\tso2mm\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 mass crystals\tso2mc\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so2 mass total\tso2mtot\tkg\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hso4m content fluid\thso4mCf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hso4m content liquid\thso4mCl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hso4m content vapor\thso4mCv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hso4m content melt\thso4mCm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous hso4m content fluid\thso4mCfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous hso4m content liquid\thso4mClp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous hso4m content vapor\thso4mCvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid hso4m mobility\thso4mlm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor hso4m mobility\thso4mvm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hso4m fraction fluid\thso4mfracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "hso4m fraction liquid\thso4mfracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "hso4m fraction vapor\thso4mfracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "hso4m concentration liquid\thso4mconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hso4m concentration vapor\thso4mconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "hso4m concentration fluid\thso4mconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so4mm content fluid\tso4mmCf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so4mm content liquid\tso4mmCl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so4mm content vapor\tso4mmCv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so4mm content melt\tso4mmCm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous so4mm content fluid\tso4mmCfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous so4mm content liquid\tso4mmClp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous so4mm content vapor\tso4mmCvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid so4mm mobility\tso4mmlm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor so4mm mobility\tso4mmvm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so4mm fraction fluid\tso4mmfracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "so4mm fraction liquid\tso4mmfracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "so4mm fraction vapor\tso4mmfracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "so4mm concentration liquid\tso4mmconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so4mm concentration vapor\tso4mmconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "so4mm concentration fluid\tso4mmconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "h2o activity liquid\tah2ol\tmol kg-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "h2o activity vapor\tah2ov\tmol kg-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "ph\tph\t-\t1\t-1.00E+02\t1.00E+02\tnode" << std::endl;
        variables_file << "sulfidation state\tsulfstate\t-\t1\t0.00E+00\t1.00E+02\tnode" << std::endl;
        variables_file << "h2s hsm ratio liquid\th2shsmratio\tmol/mol\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "log h2 fugacity\tloghh2\tbar\t1\t-1.00E+5\t1.00E+5\tnode\t" << std::endl;
    }


    //! Lithium Model
    if (with_lithium)
    {
        variables_file << "partition coefficient lithium fluid melt\tDfmli\t-\t1\t0.00E+00\t1.00E+30\tmodel" << std::endl;
        variables_file << "partition coefficient lithium crystal melt\tDcmli\t-\t1\t0.00E+00\t1.00E+30\tmodel" << std::endl;
        variables_file << "partition coefficient lithium vapor liquid slope\tdDvlli\t-\t1\t-1.00E+30\t1.00E+30\tmodel" << std::endl;
        variables_file << "li concentration melt initial\tlicminit\t-\t1\t0.00E+00\t1.00E+30\tmodel\t" << std::endl;
        variables_file << "li content fluid\tliCf\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "li content liquid\tliCl\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "li content vapor\tliCv\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "liquid li mobility\tlilm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "vapor li mobility\tlivm\tm2 s-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;

        variables_file << "li content melt\tliCm\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "li content crystals\tliCc\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous li content fluid\tliCfp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous li content liquid\tliClp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "previous li content vapor\tliCvp\tkg m-3\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "li fraction fluid\tlifracf\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "li fraction liquid\tlifracl\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "li fraction vapor\tlifracv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "li concentration liquid\tliconcl\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "li concentration vapor\tliconcv\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "li concentration fluid\tliconcf\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "li concentration melt\tliconcm\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "li concentration crystals\tliconcc\tmol kg-1\t1\t0.00E+00\t1.00E+30\tnode" << std::endl;
        variables_file << "li fraction fluid rescaled\tlifracfsc\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "li mass out total\tlimot\tkg \t1\t-1.00E+30\t1.00E+30\tmodel" << std::endl;
        variables_file << "li mass liquid\tliml\tkg \t1\t0.00E+30\t1.00E+30\tnode" << std::endl;
        variables_file << "li mass vapor\tlimv\tkg \t1\t0.00E+30\t1.00E+30\tnode" << std::endl;
        variables_file << "li mass melt\tlimm\tkg \t1\t0.00E+30\t1.00E+30\tnode" << std::endl;
        variables_file << "li mass fluid\tlimf\tkg \t1\t0.00E+30\t1.00E+30\tnode" << std::endl;
        variables_file << "li mass crystals\tlimc\tkg \t1\t0.00E+30\t1.00E+30\tnode" << std::endl;
        variables_file << "li mass total output\tlimto\tkg \t1\t0.00E+30\t1.00E+30\tnode" << std::endl;
        variables_file << "partition coefficient lithium vapor liquid node\tDvllinode\t-\t1\t0.00E+30\t1.00E+30\tnode" << std::endl;
    }


    //! Passive Tracer
    //! also used with well model
    // if (with_tracer)
    {
        variables_file << "tracer content fluid\ttcf\tkg m-3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "tracer content liquid\ttcl\tkg m-3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "tracer content vapor\ttcv\tkg m-3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "liquid tracer mobility\tltm\t-\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "vapor tracer mobility\tltm\t-\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;

        variables_file << "previous tracer content fluid\tptcf\tkg m-3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "previous tracer content liquid\tptcl\tkg m-3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "previous tracer content vapor\tptcv\tkg m-3\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "tracer fraction fluid\ttfracf\t- \t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "tracer fraction fluid rescaled\ttfracfscaled\t- \t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "tracer fraction liquid\ttfracl\t- \t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "tracer fraction vapor\ttfracv\t- \t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "tracer mass out total\ttrmot\tkg \t1\t-1.00E+30\t1.00E+30\tmodel" << std::endl;
    }


    //! Magmatic degassing (i.e., MagmaModel or PointSource)
    variables_file << "with Cl partitioning\twClpart\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "magmatic fluid mass liquid\tmMl\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "magmatic fluid mass vapor\tmMv\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "magmatic fluid mass\tmMf\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "magmatic salt content liquid\tmCXl\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "magmatic salt content vapor\tmCXv\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "magmatic mass salt\tmms\tkg\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "magmatic liquid mass mobility\tmrKl\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "magmatic vapor mass mobility\tmrKv\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "magmatic liquid salt mobility\trKlm\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "magmatic vapor salt mobility\trKvm\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "magmatic fluid flux integral\tmffI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "shell magmatic fluid flux integral\tsmffI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "magmatic fluid flux\tmff\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "shell magmatic fluid flux\tsmff\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "magmatic salt content halite\tmCXh\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;

    variables_file << "sigma1 MagmaModel\tsig1\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl; // also used for CVFEM rock class crystallization

    variables_file << "previous magmatic fluid mass\tmMfp\tkg m-3\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "previous magmatic fluid mass vapor\tmMvp\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "previous magmatic fluid mass liquid\tmMlp\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "previous magmatic salt content liquid\tmCXlp\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous magmatic salt content vapor\tmCXvp\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous magmatic salt content halite\tmCXhp\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous magmatic mass salt\tpmms\tkg\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "previous magmatic fluid flux integral\tpmffI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "previous shell magmatic fluid flux integral\tpsffI\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "previous salinity magmatic fluids\tpSmf\tWt. NaCl\t1\t0.00E+00\t1.00E+02\tnode" << std::endl;

    if (with_degassing)
    {
        variables_file << "crystal density\tnRHOc\tkg m-3\t1\t1.00E-05\t1.00E+04\tmodel" <<std::endl;
        variables_file << "compressibility magma\tBM\tPa-1\t1\t0.00E+00\t1.00E-03\tmodel" <<std::endl;
        variables_file << "power law\tpowL\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;

        variables_file << "with channel transport\twct\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl; // point source
        variables_file << "with half radial\twhr\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl; // point source
        variables_file << "nodal crystallization\tncryst\t- \t1\t0.00E+00\t1.00E+00\tmodel" << std::endl; // point source
        variables_file << "radius cupola\tr\t- \t1\t1.00E+00\t1.00E+04\tmodel" << std::endl; // point source
        variables_file << "chamber axial ratio z x\tw_z\t- \t1\t0.00E+00\t1.00E+00\tmodel" << std::endl; // point source

        variables_file << "porous flag node\tphi_flag\t-\t1\t0.\t1.E+9\tnode" << std::endl;
        variables_file << "porous flag element\tphi_flag\t-\t1\t0.\t1.E+9\telement" << std::endl;

        variables_file << "magmatic ratio\t-\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "magmatic water ratio\t-\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "magmatic salt ratio\t-\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;

        variables_file << "magmatic liquid flux in\tmlfi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "magmatic liquid flux out\tmlfo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "magmatic vapor flux in\tmvfi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "magmatic vapor flux out\tmvfo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "shell liquid flux in\tslfi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "shell liquid flux out\tslfo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "shell vapor flux in\tsvfi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "shell vapor flux out\tsvfo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "NaCl KCl ratio\t-\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;

        variables_file << "magmatic production rate\tmpr\tkg s-1\t1\t0.00E+00\t1.00E+20\tmodel" << std::endl;
        variables_file << "magmatic injection rate\tmir\tkg s-1\t1\t0.00E+00\t1.00E+20\tmodel" << std::endl;
        variables_file << "total magmatic fluid\ttmf\tkg\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "total fluid injected\ttfi\tkg\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "magma pressure\tmPL\tPa\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "injection location\til\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;

        variables_file << "volume deformation\tVdef\tm3\t1\t-1.00E+20\t1.00E+20\tnode" << std::endl;
        variables_file << "total volume deformation\tVdef_tot\tm3\t1\t-1.00E+20\t1.00E+20\tnode" << std::endl;
        variables_file << "rock density scaling\trhoR_scale\tm3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "mini nodal density rock\tmin_rhoR\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "mini pore volume\tminVp\tm3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "mini mass fluid\tmin_mF\tkg\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "melt volume\tVm\tm3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "crystal volume\tVc\tm3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "melt volume change\tdelta_Vm\tm3\t1\t-1.00E+20\t1.00E+20\tnode" << std::endl;
        variables_file << "crystal volume change\tdelta_Vc\tm3\t1\t-1.00E+20\t1.00E+20\tnode" << std::endl;
        variables_file << "crystallinity\tcryst\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "crystallinity element\tcryst_e\t-\t1\t0.00E+00\t1.00E+00\telement" << std::endl;
        variables_file << "delta crystallinity\tdelta_cryst\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "crystallized chamber flag\tcc_flag\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
        variables_file << "crystallized chamber flag element\tcc_flag_e\t-\t1\t0.00E+00\t1.00E+00\telement" << std::endl;
        variables_file << "mass volatile produced\tmV\tkg\t1\t-1.00E+20\t1.00E+20\tnode" << std::endl;
        variables_file << "total mass volatile produced\tmV_tot\tkg\t1\t-1.00E+20\t1.00E+20\tnode" << std::endl;
        variables_file << "total mass volatile produced model\tmV_tot_model\tkg \t1\t-1.00E+20\t1.00E+20\tmodel" << std::endl;
        variables_file << "total mass volatile injected model\tmV_tot_model\tkg \t1\t-1.00E+20\t1.00E+20\tmodel" << std::endl;
        variables_file << "integrated mass volatile produced model\tmV_model_int\tkg \t1\t-1.00E+20\t1.00E+20\tmodel" << std::endl;


        variables_file << "magma density\trhoMag\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "melt density\trhoM\tkg m-3\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "temperature difference\tdelta_T\toC\t1\t-1.00E04\t1.00E+04\tnode" << std::endl;
        variables_file << "volatile dissolved fraction\tVdiss_frac\t-\t1\t0.\t1.E+9\tnode" << std::endl;
        variables_file << "volatile dissolved mass\tVdiss_m\tkg\t1\t0.\t1.E+9\tnode" << std::endl;
        variables_file << "volatile dissolved eq\tVdiss_eq\t-\t1\t0.\t1.\tnode" << std::endl;
        variables_file << "volatile saturation\tVsat\t-\t1\t0.\t1.E+9\tnode" << std::endl;
        variables_file << "relative saturation\tVsatr\t-\t1\t0.\t1.E+0\tnode" << std::endl;

        variables_file << "crystal fraction\tcryst_frac\t-\t1\t0.\t1.\tnode" << std::endl;
        variables_file << "crystal fraction element\tcryst_frac\t-\t1\t0.\t1.\telement" << std::endl;
        variables_file << "melt fraction\tfracM\t-\t1\t0.\t1.\tnode" << std::endl;
        variables_file << "melt fraction element\tfracMe\t-\t1\t0.\t1.\telement" << std::endl;
        variables_file << "volatile fraction\tfracV\t-\t1\t0.\t1.\tnode" << std::endl;
        variables_file << "volatile fraction element\tfracV\t-\t1\t0.\t1.\telement" << std::endl;
        variables_file << "volatile fraction critical element\tV_frac_crit\t-\t1\t0.\t1.E+9\telement" << std::endl;
        variables_file << "relative permeability\tkcti\tm2\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "channels permeability\tkcti\tm2\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "intrinsic permeability\tkcti\tm2\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "crystallization reference\tcf\t-\t1\t0.00E+00\t3.00E+01\telement" << std::endl;

        variables_file << "previous magmatic liquid flux in\tmlfi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "previous magmatic liquid flux out\tmlfo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "previous magmatic vapor flux in\tmvfi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "previous magmatic vapor flux out\tmvfo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "previous shell liquid flux in\tslfi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "previous shell liquid flux out\tslfo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "previous shell vapor flux in\tsvfi\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "previous shell vapor flux out\tsvfo\tkg m-3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
        variables_file << "new magmatic fluid\tnmf\tkg\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;

        variables_file << "half radial volume\thrV\tm3\t1\t0.00E+00\t1.00E+20\telement" << std::endl;
        variables_file << "nodal crystallization reference\tcf\t-\t1\t0.00E+00\t30\tnode" << std::endl;
        variables_file << "nodal half radial volume\thrV\tm3\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;

        variables_file << "chlorine mass melt\tclmm\t-\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "chlorine mass fluid\tclmf\t-\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "chlorine mass crystals\tclmc\t-\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "chlorine concentration fluid\tclcf\t-\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "chlorine concentration crystals\tclcc\t-\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "chlorine concentration melt\tclcm\t-\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
        variables_file << "partition coefficient chlorine fluid melt\tDClfm\t-\t1\t0.00E+00\t1.00E+10\tnode" << std::endl;
    }

    if (with_degassing || with_gold)
    {
        variables_file << "melt mass\tmM\tkg\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "crystal mass\tmC\tkg\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
        variables_file << "magmatic fluid salinity\tSmf\tsal_magmaf\t1\t0.00E+00\t1.00E+02\tmodel" << std::endl;
    }

    //! Split boundaries
    variables_file << "with split boundary\twsb\t-\t1\t0.00+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "split region gravity mass source\tgrav_mass_source\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
    if (with_split_boundary){
        variables_file << "open space\tos\tm\t1\t0.00E00\t1.0E+00\tmodel" << std::endl;
        variables_file << "fault permeability anisotropy factor\tkfaf\tm\t1\t1.00E-10\t1.00E+05\tmodel" << std::endl;

        variables_file << "thickness split boundary\thsb\tm\t1\t1.00E-25\t1.00E+03\tmodel" << std::endl;
        variables_file << "split region nodal area\tnA\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
        variables_file << "split region mass source\tqm\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
        variables_file << "split region energy source\tqe\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
        variables_file << "split region salt source\tqs\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
        variables_file << "split region potential gradient l\tsplit_grad_pot_l\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
        variables_file << "split region potential gradient v\tsplit_grad_pot_v\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
        variables_file << "split region relperm l\tsplit_relperm_l\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
        variables_file << "split region relperm v\tsplit_relperm_v\t-\t1\t0.00+00\t1.00E+20\tnode" << std::endl;
        variables_file << "leftover nfvs\tlonfvs\t-\t1\t-1.00E+20\t1.00E+20\tnode" << std::endl;
        variables_file << "previous leftover nfvs\tplonfvs\t-\t1\t-1.00E+20\t1.00E+20\tnode" << std::endl;
    }

    variables_file << "specific enthalpy liquid top\trhlT\tJ kg-1\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "specific enthalpy vapor top\trhvT\tJ kg-1\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;

    variables_file << "with depth dependent fault permeability\twkfD\t-\t1\t0.00E+00\t1.00E+00\tmodel" << std::endl;
    variables_file << "minimum permeability fault\tkf_min\tm2\t1\t1.00E-30\t1.00E-01\telement" << std::endl;


    //! With air phase
    variables_file << "nodal compressibility air\tnBa\tPa-1\t1\t1.00E-17\t1.00E-03\tnode" << std::endl;
    variables_file << "gas mix compressibility\tBgm\tPa-1\t1\t0.00E+00\t1.00E-03\tnode" << std::endl;
    variables_file << "partial pressure air\tpa\tPa\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "fluid mass air\tMa\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "previous fluid mass air\tpMa\tkg m-3\t1\t0.00E+00\t3.00E+09\tnode" << std::endl;
    variables_file << "saturation air\tSa\tX\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "density air\tRHOa\tkg m-3\t1\t0.00E+00\t2.00E+03\tnode" << std::endl;
    variables_file << "gas mix density\tRHOgm\tkg m-3\t1\t0.00E+00\t2.00E+03\tnode" << std::endl;
    variables_file << "density air transport\tRHOaT\tkg m-3\t1\t1.00E+00\t2.00E+08\tnode" << std::endl;
    variables_file << "viscosity air\tMUa\tPa s\t1\t1.00E-07\t1.00E-03\tnode" << std::endl;
    variables_file << "gas mix viscosity\tMUgm\tPa s\t1\t1.00E-07\t1.00E-03\tnode" << std::endl;
    variables_file << "enthalpy air\tHa\tJ kg-1\t1\t1.00E+02\t1.00E+20\tnode" << std::endl;
    variables_file << "gas mix enthalpy\tHgm\tJ kg-1\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "volumetric enthalpy air\tVHa\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "enthalpy content air\tCHa\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "previous enthalpy content air\tPCHa\tJ m-3\t1\t0.00E+00\t2.00E+20\tnode" << std::endl;
    variables_file << "air mass mobility\tmoba\tm2 s-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "air mass mobility density\tmobda\tkg s-1 m-1\t1\t0.00E+00\t1.00E+20\tnode" << std::endl;
    variables_file << "air enthalpy mobility\thKa\tm2 s-1\t1\t1.00E-25\t1.00E+10\tnode" << std::endl;
    variables_file << "air enthalpy mobility density\thKad\tm2 s-1\t1\t1.00E-25\t1.00E+10\tnode" << std::endl;
    variables_file << "relperm viscosity air\tMUa_inv\ts Pa-1\t1\t1.00E-07\t1.00E+07\tnode" << std::endl;
    variables_file << "velocity air\tVa\tm s-1\t2\t-1.00E+03\t1.00E+03\telement" << std::endl;
    variables_file << "pore velocity air\tPVa\tm s-1\t2\t-1.00E+03\t1.00E+03\telement" << std::endl;
    variables_file << "nodal source air\tnSa\tkg kg-1 s-1\t1\t-1.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "boundary flow air\tbfa\tkg m-3\t1\t-1.00E+09\t1.00E+09\tnode" << std::endl;
    variables_file << "courant air\tcfl_a\ts\t1\t1.00E-10\t1.00E+20\telement" << std::endl;
    variables_file << "upwind control air\tuc_a\t-\t1\t0\t30\telement" << std::endl;
    variables_file << "gas mix molar fraction air\tmolfa\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "gas mix molar fraction vapor\tmolfv\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "specific enthalpy air top\trhaT\tJ kg-1\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "mass fraction inflow vapor in air top\tmfvT\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "mass fraction inflow air top\tmfaT\t-\t1\t0.00E+00\t1.00E+00\tnode" << std::endl;
    variables_file << "specific enthalpy liquid top\trhlT\tJ kg-1\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "specific enthalpy vapor top\trhvT\tJ kg-1\t1\t0.00E+00\t1.00E+07\tnode" << std::endl;
    variables_file << "rain recharge top\train\tmm yr-1\t1\t0.00E+00\t1.00E+05\tnode" << std::endl;
    variables_file << "partial pressure air\tpa\tPa\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;
    variables_file << "partial pressure vapor\tpv\tPa\t1\t0.00E+00\t1.00E+09\tnode" << std::endl;


    variables_file.close();
}

template<size_t dim>
void ModelBuilder<dim>::CreateLimiterLogFile()
{
    std::ofstream LimiterLog_file("PressureLimiter.log");
    LimiterLog_file << "PressureLimiter" <<  std::endl << std::endl;
    LimiterLog_file << "Limiter type\tmodel time\tx-coord\ty-coord\tz-coord\tfluid pressure\ttemperature"
                    << std::endl << std::endl;
    LimiterLog_file.close();
}

template<size_t dim>
void ModelBuilder<dim>::DeleteVariablesFile()
{
#ifdef WIN32
    system("del PhysicalVariables.txt");
#else
    system("rm PhysicalVariables.txt");
#endif
}

template<size_t dim>
ModelBuilder<dim>::~ModelBuilder()
{
}

template class ModelBuilder<1U>;
template class ModelBuilder<2U>;
template class ModelBuilder<3U>;

} // end namespace csmp
