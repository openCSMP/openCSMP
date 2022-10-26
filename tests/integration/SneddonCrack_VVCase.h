#ifndef SNEDDON_CRACK_VVCASE_H
#define SNEDDON_CRACK_VVCASE_H

#include "Test.h"
#include "GlobalVerbose.h"

#include "CSMP_highLevelUtilities.h"
#include "ANSYS_Model2D.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_Integrator.h"
//#include "NumIntegral_BT_C_B_dV.h"
#include "NumIntegral_BT_D_B_dV.h"
//#include "StressAndStrainOutput.h"
#include "ExtractTensorVariableComponent.h"
#include "NumIntegral_PT_op_dS.h"
#include "ModelSubDomain.h"
#include "Fracture.h"
#include <stdio.h>
#include <functional>
#include <cmath>

//Interrellations for conductivity operator
#include "InterFaceFractureVisitor.h"

//Numerical Integrals Lubricaiton Equation
#include "NumIntegral_dNT_mixed_op_dN_dV.h"           // to interpolate the aperture
#include "NumIntegral_dNT_op_dN_dV.h"                 // to interpolate the aperture cubed  (try with quadratic base functions...)
#include "NumIntegral_NT_lhsop_N_dV.h"                // Mass Matrix LHS
#include "NumIntegral_SetRHS_to_Zero.h"               // Zero right hand side
#include "NumIntegral_NT_op_N_dV.h"                   // Lumped Mass Matrix RHS or Source term with AccumulateLater()
#include "NumIntegral_PT_op_dV.h"
#include "PT_op.h"                                    //b force num int

//#include "RecoveryBasedOnDisplacement.h"

namespace csmp
  {

  class SneddonCrack_VVCase: public Test
    {
    public:
      SneddonCrack_VVCase(const char* prefix);
      virtual void run();

      void set_min_avg_error(double epsilon){ min_avg_error = epsilon; };


    private:

      double min_avg_error = 8.0;

      template<typename T>
      void WriteSolutionToFile( std::vector<std::vector<T>> data, std::string name);

      template<size_t dim, typename T>
      void Verify(std::map<Point<dim>,double>& PressureProfile, bool at_angle, T F , double tol){
        double avg_percent_error = 0.0;
        double length = 0.0;
        if (!at_angle){                                                                  //If not anglend use this as verification
          for (typename std::map<Point<dim>,double>::iterator it = PressureProfile.begin(); it != PressureProfile.end(); it++){
              double x = it->first[0];
              length = 2.0;
              if (F(x,length) == 0.0 || it->second == 0.0)
                continue;

              std::cout << "\nPos: "  << it->first[0] << " Numerical  pressure: " << it->second << std::endl;
              std::cout << "Pos: "    << it->first[0] << " Analytical pressure: " << F(x, length) << std::endl;
              std::cout << "Average Normalised Error: " << it->second/F(x,length) - 1.0 << std::endl;
              _equal(it->second/F(x,length), 1.0, tol);
              avg_percent_error += std::abs( it->second - F(x,length))/F(x,length);
              }
          }
        else {                                                                          //More accurate angled version of error calculation
            const Point<dim> midpoint (0.0,0.0);
            typename std::map<Point<dim>,double>::iterator mid_i = PressureProfile.find(midpoint);
            assert(mid_i != PressureProfile.end());

            for (typename std::map<Point<dim>,double>::iterator it = PressureProfile.begin(); it != --(PressureProfile.end()); 0){
                auto it_before = it++;
                length += (it_before->first).DistanceTo((it)->first);                 //calculate length of domain
              }


            for (typename std::map<Point<dim>,double>::iterator it = PressureProfile.begin(); it != PressureProfile.end(); it++){
                double x = 0.0;
                for (typename std::map<Point<dim>,double>::iterator point_it = it; point_it != mid_i; 0 ){
                  if (it->first[0] < 0) {
                        auto point_it_before = point_it++;
                        x -= point_it_before->first.DistanceTo(point_it->first);
                      }
                  else if (it->first[0] > 0) {
                      auto point_it_before = point_it--;
                      x += point_it_before->first.DistanceTo(point_it->first);
                    }
                  else {
                      x = 0.0;
                    }
                  }
                if (F(x,length) == 0.0 || it->second == 0.0)
                  continue;
                std::cout << "\nPos: "  << it->first[0] << " Numerical  pressure: " << it->second << std::endl;
                std::cout << "Pos: "    << it->first[0] << " Analytical pressure: " << F(x,length) << std::endl;
                std::cout << "Average Normalised Error: " << it->second/F(x,length) - 1.0 << std::endl;
                _equal(it->second/F(x,length), 1.0, tol);
                avg_percent_error += std::abs( it->second - F(x,length))/F(x,length);
                }
          }

        std::cout << "\n---------------------------------------\n"
                   "Verify::Average Normalised Error: " << avg_percent_error/(PressureProfile.size()) << std::endl;
     }



    };




  } // csmp

#endif
