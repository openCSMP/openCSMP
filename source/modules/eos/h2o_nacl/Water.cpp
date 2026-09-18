// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Water.h"

using namespace std;

namespace csmp
{
  /** constructor
      
      Water needs to be constructed with a reference to the temperature [C] and pressure [Pa] variables as constructor argument, otherwise it will not work properly. Default construction has been disbaled by making the constructor private.
      
      @todo read up on smart pointers in order to be able to do safe installation of PROST newProp calls into initializer list. Di, the initialization of the ak and bij arrays is still in the constructor body.
  */
  Water::Water(const double& externaltemperature, 
               const double& externalpressure)
    : temperature(externaltemperature),
      pressure(externalpressure),
      kelvin(273.15e0),
      mol_h2o(18015.0e0),
      tcurrent(1.0), 
      pcurrent(1.0),
      d(500.0e0),
      dp(1.0e-08),
      dvdppascal(1.0), 
      comprpascal(1.0), 
      delp(1.0e-2),
      delt(1.0e-2),
      dvdt(1.0), 
      viscosity(1.0)
  {
    ak[0] = 0.0181583, ak[1] = 0.0177624, ak[2] =  0.0105287, ak[3] = -0.0036744;
    bij[0][0] = 0.501938, bij[1][0] =  0.162888, bij[2][0] = -0.130356, bij[3][0] =  0.907919,
	    bij[4][0] = -0.551119, bij[5][0] =  0.146543, bij[0][1] =  0.235622,  bij[1][1] =  0.789393,
	    bij[2][1] =  0.673665, bij[3][1] =  1.207552, bij[4][1] =  0.0670665, bij[5][1] = -0.084337,
	    bij[0][2] = -0.274637, bij[1][2] = -0.743539, bij[2][2] = -0.959456,  bij[3][2] = -0.687343,
	    bij[4][2] = -0.497089, bij[5][2] =  0.195286, bij[0][3] =  0.145831,  bij[1][3] =  0.263129,
	    bij[2][3] =  0.347247, bij[3][3] =  0.213486, bij[4][3] =  0.100754,  bij[5][3] = -0.032932,
	    bij[0][4] = -0.0270448,bij[1][4] = -0.0253093,bij[2][4] = -0.0267758, bij[3][4] = -0.0822904,
	    bij[4][4] =  0.0602253,bij[5][4] = -0.0202595;
    properties         = newProp('T', 'p', 1);
    delPproperties     = newProp('T', 'p', 1);
    delTproperties     = newProp('T', 'p', 1);
    minTproperties     = newProp('T', 'p', 1);
    liqpropsp          = newProp('T', 'p', 1);
    vappropsp          = newProp('T', 'p', 1); 
    liqpropst          = newProp('T', 'p', 1);
    vappropst          = newProp('T', 'p', 1);
    liqpropstplusdelt  = newProp('T', 'p', 0);
    vappropstplusdelt  = newProp('T', 'p', 0);
    liqpropstminusdelt = newProp('T', 'p', 0);
    vappropstminusdelt = newProp('T', 'p', 0);
  }
  

  /** destructor
      Default destruction would probably not work. Main taks is calling PROST's freeProp function
  */
  Water::~Water()
  {
    freeProp(properties);
    freeProp(delPproperties);
    freeProp(delTproperties);
    freeProp(minTproperties);
    freeProp(liqpropsp);
    freeProp(vappropsp);
    freeProp(liqpropst);
    freeProp(vappropst);
    freeProp(liqpropstplusdelt);
    freeProp(vappropstplusdelt);
    freeProp(liqpropstminusdelt);
    freeProp(vappropstminusdelt);
  }

  /** update temperature and pressure to current values and fill PROST's property structures with the correct values; in addition compute viscosity for new conditions

      @todo the viscosity call seems error-prone close to the saturation curve as no phase state determination has been done when it is called
  */
  void Water::CheckStatus()
  {
    tcurrent = temperature;
    pcurrent = pressure;
    water_tp( tcurrent+273.15, pcurrent, d, dp, properties );
    water_tp( tcurrent+273.15, pcurrent+delp*1.0e5, d, dp, delPproperties );
    water_tp( tcurrent+273.15+delt, pcurrent, d, dp, delTproperties );
    water_tp( tcurrent+273.15-delt, pcurrent, d, dp, minTproperties );
    viscosity = Viscosity(tcurrent+273.15,properties->d);
      
    if(tcurrent >= cp_h2o.Temperature())
      {
        dvdppascal  = ( 18015.0e0/delPproperties->d - (18015.0e0/properties->d) )/delp/1.0e5;
        comprpascal = -1.0e0/(18015.0e0/properties->d)*dvdppascal;
        dvdt        = ( 18015.0e0/delTproperties->d - 18015.0e0/properties->d )/delt;
        tcurrent    = temperature;
        pcurrent    = pressure;
        if(tcurrent-delt < cp_h2o.Temperature())
          {
            sat_t( tcurrent+273.15-delt, liqpropstminusdelt, vappropstminusdelt );
          }
        else
          {
            water_td(tcurrent+273.15-delt,cp_h2o.Density(),liqpropstminusdelt); 
            water_td(tcurrent+273.15-delt,cp_h2o.Density(),vappropstminusdelt);
          }
        return;
      }
    else
      { // subcritical
        sat_p( pcurrent, liqpropsp, vappropsp );
        sat_t( tcurrent+273.15, liqpropst, vappropst );
        
        if(pcurrent == vappropst->p || pcurrent == liqpropst->p )
          {
            dvdppascal = (18015.0e0/delPproperties->d - 18015.0e0/liqpropst->d)/
              (delPproperties->p - pcurrent);
            comprpascal = -liqpropst->d/18015.0e0*dvdppascal;
            water_tp( tcurrent+273.15-delt, pcurrent, d, dp, delTproperties );
            dvdt        =  ( -18015.0e0/delTproperties->d + 18015.0e0/liqpropst->d )/delt;
            tcurrent    = temperature;
            pcurrent    = pressure;
            water_tp( tcurrent+273.15+delt, pcurrent, d, dp, delTproperties );
            sat_t( tcurrent+273.15-delt, liqpropstminusdelt, vappropstminusdelt );
            if(temperature+delt < cp_h2o.Temperature())
              sat_t( tcurrent+273.15+delt, liqpropstplusdelt, vappropstplusdelt );
            else
              { // metastable extension of saturation curve = critical isochore
                water_td(tcurrent+273.15+delt,cp_h2o.Density(),liqpropstplusdelt); 
                water_td(tcurrent+273.15+delt,cp_h2o.Density(),vappropstplusdelt);
              } 
            return;
          }
        else if(pcurrent < vappropst->p && pcurrent+delp*1.0e5 >= vappropst->p)
          {
            // vapour near saturation line !!!
            // Careful, as T approaches tcrit, the delp has to be adjusted
            // because the saturation curve is not precisely calculated by prost
            dvdppascal = (18015.0e0/vappropst->d - 18015.0e0/properties->d)/
              (vappropst->p - pcurrent);
            comprpascal = -1.0e0/(18015.0e0/properties->d)*dvdppascal;
            dvdt        = ( 18015.0e0/delTproperties->d - 18015.0e0/properties->d )/delt;
            tcurrent    = temperature;
            pcurrent    = pressure;
            sat_t( tcurrent+273.15-delt, liqpropstminusdelt, vappropstminusdelt );
            if(temperature+delt < cp_h2o.Temperature())
              sat_t( tcurrent+273.15+delt, liqpropstplusdelt, vappropstplusdelt );
            else
              { // metastable extension of saturation curve = critical isochore
                water_td(tcurrent+273.15+delt,cp_h2o.Density(),liqpropstplusdelt); 
                water_td(tcurrent+273.15+delt,cp_h2o.Density(),vappropstplusdelt);
              } 
            return;
          }
        else if(pcurrent > liqpropst->p && pcurrent-delp*1.0e5 <= liqpropst->p)
          {
            // liquid near saturation line !!!
            // Careful, as T approaches tcrit, the delp has to be adjusted
            // because the saturation curve is not precisely calculated by prost
            dvdppascal = ((18015.0e0/properties->d) - 18015.0e0/liqpropst->d)/
              (pcurrent-liqpropst->p);
            comprpascal = -1.0e0/(18015.0e0/properties->d)*dvdppascal;
            water_tp( tcurrent+273.15-delt, pcurrent, d, dp, delTproperties );
            dvdt        =  ( -18015.0e0/delTproperties->d + 18015.0e0/properties->d )/delt;
            tcurrent    = temperature;
            pcurrent    = pressure;
            water_tp( tcurrent+273.15+delt, pcurrent, d, dp, delTproperties );
            sat_t( tcurrent+273.15-delt, liqpropstminusdelt, vappropstminusdelt );
            if(temperature+delt < cp_h2o.Temperature())
              sat_t( tcurrent+273.15+delt, liqpropstplusdelt, vappropstplusdelt );
            else
              { // metastable extension of saturation curve = critical isochore
                water_td(tcurrent+273.15+delt,cp_h2o.Density(),liqpropstplusdelt); 
                water_td(tcurrent+273.15+delt,cp_h2o.Density(),vappropstplusdelt);
              } 
            return;
          }
        else
          {
            // liquid or vapour !!!
            dvdppascal = (18015.0e0/delPproperties->d - 18015.0e0/properties->d)/delp/1.0e5;
            comprpascal = -1.0e0/(18015.0e0/properties->d)*dvdppascal;
            dvdt        = ( 18015.0e0/delTproperties->d - 18015.0e0/properties->d )/delt;
            tcurrent    = temperature;
            pcurrent    = pressure;
            sat_t( tcurrent+273.15-delt, liqpropstminusdelt, vappropstminusdelt );
            if(temperature+delt < cp_h2o.Temperature())
              sat_t( tcurrent+273.15+delt, liqpropstplusdelt, vappropstplusdelt );
            else
              { // metastable extension of saturation curve = critical isochore
                water_td(tcurrent+273.15+delt,cp_h2o.Density(),liqpropstplusdelt); 
                water_td(tcurrent+273.15+delt,cp_h2o.Density(),vappropstplusdelt);
              } 
            return;
          }
      }
  }
    

  /** dynamic visocsity of water

      @todo This is an outdated formulation (probably 1985 IAPS), accurate enough for all simulations but update to most recent IAPWS release would be good
  */
  double Water::Viscosity(const double& T, const double& rho)
  {
    const double tstar(647.27e0), rhostar(317.763e0);
    double trat,trat1,rhorat,rhorat1,n0,n;
    int i,j,k;
    
    trat    = T/tstar;
    trat1   = 1.0e0/trat-1.0e0;
    rhorat  = rho/rhostar;
    rhorat1 = rhorat-1.0e0;
    
    n0=0;
    for(k=0;k<4;++k)
      {
        n0 += (ak[k]*pow((1.0e0/trat),k));
      }  
    n0 = 1.0e0/n0;
    n0 *= (sqrt(trat));
    
    n=0;
    for(i=0;i<6;++i)
      {
        for(j=0;j<5;++j)
          {
            n += (bij[i][j]*pow(trat1,i)*pow(rhorat1,j));
          }
      }
    n = n0*exp(rhorat*n);
    
    return n*1.0e-06; // converts to kg/s/m
  }
  

}//csmp
