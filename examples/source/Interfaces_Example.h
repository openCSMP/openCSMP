// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_INTERFACES_EXAMPLE_H
#define CSMP_INTERFACES_EXAMPLE_H

#include "Example.h"
#include "Model.h"


namespace csmp {

class  Interfaces_Example : public Example
{

public:
    virtual void Run();
    virtual void Specifications();

private:
    void BuildFromANSYS2DModel(const std::string& model_name);
    void BuildFromANSYS3DModel(const std::string& model_name);
    void BuildFromTRIANGLEModel(const std::string& model_name);
    void BuildFromTriangulatorModel(const std::string& model_name, double extent1, double extent2);
    void BuildFromQuadrilateratorModel(const std::string& model_name, double x, double y);
    template<uint32_t dim> void BuildFromSKUAModel(const std::string& model_name);
    void BuildFromEclipseModel(const std::string& model_name);
    void BuildFromGeoModellerModel(const std::string& model_name);
    void BuildFromRhinoModel(const std::string& model_name);
    template<uint32_t dim> void OutputRegionIDToVTU (Model<dim>& model);

    bool output_vtu_ = false;
};

} // csmp

#endif // CSMP_INTERFACES_EXAMPLE_H
