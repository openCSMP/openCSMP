// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_STATISTICAL_DISTRIBUTION_GENERATOR_H
#define CSMP_STATISTICAL_DISTRIBUTION_GENERATOR_H

#include "Matrix.h"

namespace csmp {

/**
       Generation of various statistical distributions and output of corresponding histograms to file.
       
@author Siroos Azizmohammadi

@note old name MathLibrary.

*/
class StatisticalDistributionGenerator
{
public:
    static double RandomNormalNonTrigonometric( double mu = 0.0, double sigma = 1.0 );
    static double RandomNormalTrigonometric( double mu = 0.0, double sigma = 1.0 );
    static double RandomNormal( double mu = 0.0, double sigma = 1.0 );
    static double RandomNormalLimit( double mu, double sigma, double lower, double upper );
    static double RandomLogNormal( double mu = 0.0, double sigma = 1.0 );
    static double RandomLogNormalLimit( double mu, double sigma, double lower, double upper );

    void OutputHistogramToFile( std::string fileName );
    
    static void         EigenDecomposition( Matrix A, Matrix& V, std::vector<double>& d );
    void                Histogram( std::vector<double>& prop, const size_t no_bins );
    std::vector<double> Bins() const;
    std::vector<size_t> Frequencies() const;
    

private:
    size_t abs(size_t n) { return n > 0 ? n : -n; }

    std::vector<double> hist_bins_;
    std::vector<size_t> hist_freqs_ ;
};

} // csmp

#endif // CSMP_STATISTICAL_DISTRIBUTION_GENERATOR_H
