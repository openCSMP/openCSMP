//
// Created by Michael Liem on 08.06.22.
//

#ifndef CSMP_API_MEASUREDATAFORENKF_H
#define CSMP_API_MEASUREDATAFORENKF_H

#include "CSMP_definitions.h"

// the CSMP model
#include "Model.h"
#include "Region.h"
#include "Boundary.h"

// output interfaces
#include "TextFileIO.h"
#include "TextInterface.h"
#include "VTK_Interface.h"

// Transport scheme
#include "NodeCenteredFiniteVolumeTransport.h"

#include <filesystem>

namespace csmp {
  class MeasureDataForEnKF {

  public:

    static void MeasurePressure(                   Model<2U> &model, const std::string &region_str, const std::string &output_path);
    static void MeasureVolumeFlowRate(             Model<2U> &model, const std::string &region_str, const std::string &output_path);
    static void MeasureVolumeFlowRateAlongFracture(Model<2U> &model, const std::string &output_path);
    static void MeasureConcentration(              Model<2U> &model, const std::string &region_str, const std::string &output_path, size_t timestep);
    static void MeasurePressureAndConcentration(   Model<2U> &model, const std::string &region_str, const std::string &output_path, size_t timestep);

    static void WriteSteadyStateVTK(   Model<2U>& model, const std::string& output_path );
    static void WriteConcentrationVTK( Model<2U>& model, const std::string& output_path, int time_step );
    static void WriteFractureNumbers(  Model<2U>& model, const std::string& output_path );

    static void CheckMatlabConsistency( Model<2U>& model, const std::string& output_path );

    static void PrintOutOfRangeConcentration( Model<2U>& model, const std::string & region_str, const std::string &output_path, size_t timestep, const double c_min=0., const double c_max=1. );

    static void CreateRectangularHollowRegion( Model<2U> &model, const std::string & region_str, Point<2U> x1, Point<2U> x2, double_t thickness);
    static void CreateMeasRegion( Model<2U> &model, const int meas_region_case );

  };
}


#endif //CSMP_API_MEASUREDATAFORENKF_H
