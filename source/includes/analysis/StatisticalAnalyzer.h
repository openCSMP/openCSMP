#ifndef CSMP_STATISTICAL_ANALYZER_H
#define CSMP_STATISTICAL_ANALYZER_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Model;

typedef std::vector<std::pair<double64,double64> > HistogramBins;

/** 
     @brief Class for the generation of volume-weighted and other histograms from CSMP models and their unique regions.
     
     @todo refactor such that Maple specific output is done by MapleInterface class.
     
     @uthor Stephan Matthai
     @date 2001
*/
template<size_t dim>
class StatisticalAnalyzer {
  public:
    StatisticalAnalyzer( const Model<dim>& );
    ~StatisticalAnalyzer();
    
    /// specify via ascii input file, the x-axis range of the histogram columns, i.e. the bin size
    void DefineBins( const char* bin_ascii_file, HistogramBins& ) const;

    /// specify the binsize using bin ranges specified by a sequence of double values of user-defined length
    void DefineBins( HistogramBins& bins, double64 first_val, ... ) const;

    /// creats histogram bins by providing the data range and number of devisions
    void DefineBins( const double64 minimum, 
                     const double64 maximum, 
                     const size_t number_of_bins,
                     HistogramBins& data ) const;
                     
    /// normalizes number of counts by element volume
    void RegionPropertyHistograms( const char* prop,          ///< analysed element variable
                                   const HistogramBins& bins, ///< bin ranges from < to <=
                                   //       region-name           bin-end value, counts     total counts 
                                   std::map<std::string,std::pair<HistogramBins, size_t> >& results ) const;

    /// normalize by counts; @attention this makes sense only for regular grids
    void RegionPropertyHistograms( const char* prop, // desired output property
                                  // bin ranges from < to <=
                                  const HistogramBins& bins,
                                  //   groupname          bin-end value, counts    total counts 
                                  std::map<std::string,std::pair<HistogramBins,size_t> >& results,
                                  // normalize to total counts = group size
                                  bool normalize ) const;
    
    /// analyzes properties element by element interpolating variables to barycenter if necessary
    void RegionPropertyHistogramsElement( const char* prop, 
                                          const HistogramBins& bins,
                                          std::map<std::string,std::pair<HistogramBins,size_t> >& results) const;
    
    /// analyzes properties element by element, integral over interpolation points of each element; written for velocity
    void RegionPropertyHistogramsIntegrationPoint( const char* prop,
                                                  const HistogramBins& bins,
                                                  std::map<std::string,std::pair<HistogramBins,size_t> >& results, const std::string& flow_domain) const;
                                  
    /// as previous method, but binning second property based on first property occurance                         
    void RegionPropertyHistogramsElementProperty2BinningBasedOnProperty1( const char* prop1,
                                                                          const char* prop2,
                                                                          const HistogramBins& bins,
                                          std::map<std::string,std::pair<HistogramBins,size_t> >& results1,
                                          std::map<std::string,std::pair<HistogramBins,size_t> >& results2,
                                          bool weighted_by_porosity = false ) const;


// OUTPUT TO FILES ---------------------------------------------------------------------------------------------------

    /// for each model region, a separate histogram is created using the user-specified bins
    void OutputRegionPropertyHistograms( const char* prop, 
                                         // bin ranges from < to <=
                                         const std::map<std::string,std::pair<HistogramBins,size_t> >&  results ) const;
                     
    /// writes Maple statlist format with weights; results must be copied and pasted into Maple workbook                                    
    void OutputRegionPropertyHistogramsMaple( const char* prop, 
                                              const HistogramBins& bins,
                                              // bin ranges from < to <=
                                              const std::map<std::string,std::pair<HistogramBins,size_t> >&  results,
                                              bool log10_of_bin_values ) const;
  
    /// for copy/paste to Maple: list of [bin-center,histogram column-height] pairs to plot curves 
    void OutputRegionPropertyAbundancePolygonsMaple( const char* prop, 
                                                     const HistogramBins& bins,
                                                     // bin ranges from < to <=
                                                     const std::map<std::string,std::pair<HistogramBins,size_t> >&  results,
                                                     bool log10_of_bin_values ) const;
 
    /// like previous method, but name can include filename etc. 
    void OutputRegionPropertyAbundancePolygonsMaple( const char* prop, 
                                                     const char* file_name_prefix,
                                                     const std::vector<std::pair<double64,double64> >& bins,
                                                     // bin ranges from < to <=
                                                     const std::map<std::string,std::pair<HistogramBins,size_t> >&  results,
                                                     bool log10_of_bin_values ) const;
private:
    const Model<dim>&             sref;
    const PropertyDatabase<dim>&  pref;  /// @todo (3-D) Why explicitly in there if access through model available
    
    void WriteHistogramToTextfile( const char* fname, const HistogramBins& hist, size_t points ) const;

    void WriteHistogramToMapleTextfile( const char* fname, const char* datasetname, 
                                        const HistogramBins& bins,
                                        const HistogramBins& hist, 
                                        size_t points,
                                        bool log10_of_bin_data=false ) const; 

    void WriteCurvePointsToMapleTextfile( const char* fname, const char* datasetname, 
                                          const HistogramBins& bins,
                                          const HistogramBins& hist, 
                                          size_t points,
                                          bool log10_of_bin_data=false ) const; 
};

} // csmp

#endif 



