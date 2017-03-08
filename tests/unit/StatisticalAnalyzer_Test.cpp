#include "StatisticalAnalyzer_Test.h"


namespace csmp {


StatisticalAnalyzer_Test::StatisticalAnalyzer_Test()
: model_ (new Model1D<1U>( "StatisticalAnalyzer_Test", "CSMP-2phase-variables.txt", 1., 100U ))

{
  fTolerance = 1.e-8;

  std::cout << "StatisticalAnalyzer_Test is creating test_histogram.bins\n";

  std::ofstream fout("test_histogram.bins");

  fout << "'test_histogram.bins' Testfile for CSMP's StatisticalAnalyzer capability\n";
  fout << "11\n0.\n1.\n2.\n3.\n4.\n5.\n6.\n7.\n8.\n9.\n10.\n";
  fout.close();


}


StatisticalAnalyzer_Test::~StatisticalAnalyzer_Test()
{
    std::remove("test_histogram.bins");
    delete model_;
}


void StatisticalAnalyzer_Test::run()
{
    StatisticalAnalyzerDefineBins();
    StatisticalAnalyzerRegionPropertyHistogramsElement();
    StatisticalAnalyzerRegionPropertyHistograms();
    StatisticalAnalyzer_TestRegionPropertyHistogramsElementProperty2BinningBasedOnProperty1();
}



void StatisticalAnalyzer_Test::StatisticalAnalyzerDefineBins()
{
    StatisticalAnalyzer<1U> statistical_analyzer(*model_);
    HistogramBins bins;

    statistical_analyzer.DefineBins("test_histogram.bins", bins);

    _test( bins.size() == 10U );
    _equal( bins[0].first, 0., fTolerance );
    _equal( bins[0].second, 1., fTolerance );
    _equal( bins[1].first, 1., fTolerance );
    _equal( bins[1].second, 2., fTolerance );
    _equal( bins[2].first, 2., fTolerance );
    _equal( bins[2].second, 3., fTolerance );
    _equal( bins[3].first, 3., fTolerance );
    _equal( bins[3].second, 4., fTolerance );
    _equal( bins[4].first, 4., fTolerance );
    _equal( bins[4].second, 5., fTolerance );
    _equal( bins[5].first, 5., fTolerance );
    _equal( bins[5].second, 6., fTolerance );
    _equal( bins[6].first, 6., fTolerance );
    _equal( bins[6].second, 7., fTolerance );
    _equal( bins[7].first, 7., fTolerance );
    _equal( bins[7].second, 8., fTolerance );
    _equal( bins[8].first, 8., fTolerance );
    _equal( bins[8].second, 9., fTolerance );
    _equal( bins[9].first, 9., fTolerance );
    _equal( bins[9].second, 10., fTolerance );

    bins.clear();
    /* Removing this particular test
      This method uses old C style for defining arguments (va_lists)
      and it may fail on some compilers
    statistical_analyzer.DefineBins(bins,  2., 4., 6., 8., 10., 12.);
    std::cout << "bins.size() is: " << bins.size() << "\n";
    _test( bins.size() == 5U );
    _equal( bins[0].first, 2., fTolerance );
    _equal( bins[0].second, 4., fTolerance );
    _equal( bins[1].first, 4., fTolerance );
    _equal( bins[1].second, 6., fTolerance );
    _equal( bins[2].first, 6., fTolerance );
    _equal( bins[2].second, 8., fTolerance );
    _equal( bins[3].first, 8., fTolerance );
    _equal( bins[3].second, 10., fTolerance );
    _equal( bins[4].first, 10., fTolerance );
    _equal( bins[4].second, 12., fTolerance );
    */
    statistical_analyzer.DefineBins(5., 10., 5, bins);

    _test( bins.size() == 5U );
    _equal( bins[0].first, 5., fTolerance );
    _equal( bins[0].second, 6., fTolerance );
    _equal( bins[1].first, 6., fTolerance );
    _equal( bins[1].second, 7., fTolerance );
    _equal( bins[2].first, 7., fTolerance );
    _equal( bins[2].second, 8., fTolerance );
    _equal( bins[3].first, 8., fTolerance );
    _equal( bins[3].second, 9., fTolerance );
    _equal( bins[4].first, 9., fTolerance );
    _equal( bins[4].second, 10., fTolerance );


}



void StatisticalAnalyzer_Test::StatisticalAnalyzerRegionPropertyHistogramsElement()
{
   Index por_key = model_->Database().StorageKey("porosity");
   Index sat_key = model_->Database().StorageKey("saturation water");
   size_t counter(1U);

   for (std::vector<Element<1U>*>::const_iterator
        eit = model_->Region("Model").ElementsBegin();
        eit != model_->Region("Model").ElementsEnd(); eit++, counter++)
   {
      (*eit)->Store(por_key, makeScalar(PLAIN, 0.01 * counter));
   }

    StatisticalAnalyzer<1U> statistical_analyzer(*model_);
    HistogramBins bins;
    std::map<std::string,std::pair<HistogramBins,size_t> > results;

    statistical_analyzer.DefineBins( 0., 1., 5, bins);
    statistical_analyzer.RegionPropertyHistogramsElement( "porosity", bins, results );

    _test( (*results.begin()).first == "Model" );
    _equal( (*results.begin()).second.first[0U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[1U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[2U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[3U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[4U].second, 0.2, fTolerance );

    for (std::vector<Node<1U>*>::const_iterator nit = model_->Region("Model").NodesBegin();
        nit != model_->Region("Model").NodesEnd(); nit++)
    {
      (*nit)->Store(sat_key, makeScalar(PLAIN, (*nit)->Coordinate()[0U]));
    }

    statistical_analyzer.RegionPropertyHistogramsElement("saturation water", bins, results);

    _test( (*results.begin()).first == "Model" );
    _equal( (*results.begin()).second.first[0U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[1U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[2U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[3U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[4U].second, 0.2, fTolerance );


}



void StatisticalAnalyzer_Test::StatisticalAnalyzerRegionPropertyHistograms()
{
   Index por_key = model_->Database().StorageKey("porosity");
   Index sat_key = model_->Database().StorageKey("saturation water");
   size_t counter(1U);

   for (std::vector<Element<1U>*>::const_iterator eit = model_->Region("Model").ElementsBegin();
        eit != model_->Region("Model").ElementsEnd(); eit++, counter++)
   {
      (*eit)->Store(por_key, makeScalar(PLAIN, 0.01 * counter));
   }

    StatisticalAnalyzer<1U> statistical_analyzer(*model_);
    HistogramBins bins;
    std::map<std::string,std::pair<HistogramBins,size_t> > results;

    statistical_analyzer.DefineBins( 0., 1., 5, bins);
    statistical_analyzer.RegionPropertyHistograms("porosity", bins, results);

    _test( (*results.begin()).first == "Model" );
    _equal( (*results.begin()).second.first[0U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[1U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[2U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[3U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[4U].second, 0.2, fTolerance );

    statistical_analyzer.RegionPropertyHistograms("porosity", bins, results, false);
    _test( (*results.begin()).first == "Model" );
    _equal( (*results.begin()).second.first[0U].second, 20, fTolerance );
    _equal( (*results.begin()).second.first[1U].second, 20, fTolerance );
    _equal( (*results.begin()).second.first[2U].second, 20, fTolerance );
    _equal( (*results.begin()).second.first[3U].second, 20, fTolerance );
    _equal( (*results.begin()).second.first[4U].second, 20, fTolerance );

    statistical_analyzer.RegionPropertyHistograms("porosity", bins, results, true);
    _test( (*results.begin()).first == "Model" );
    _equal( (*results.begin()).second.first[0U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[1U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[2U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[3U].second, 0.2, fTolerance );
    _equal( (*results.begin()).second.first[4U].second, 0.2, fTolerance );


    for (std::vector<Node<1U>*>::const_iterator nit = model_->Region("Model").NodesBegin();
        nit != model_->Region("Model").NodesEnd(); nit++)
    {
      (*nit)->Store(sat_key, makeScalar(PLAIN, (*nit)->Coordinate()[0U]));
    }

    statistical_analyzer.RegionPropertyHistograms("saturation water", bins, results, false);
    _test( (*results.begin()).first == "Model" );
    _equal( (*results.begin()).second.first[0U].second, 21, fTolerance );
    _equal( (*results.begin()).second.first[1U].second, 20, fTolerance );
    _equal( (*results.begin()).second.first[2U].second, 20, fTolerance );
    _equal( (*results.begin()).second.first[3U].second, 20, fTolerance );
    _equal( (*results.begin()).second.first[4U].second, 20, fTolerance );


    statistical_analyzer.RegionPropertyHistograms("saturation water", bins, results, true);
    _test( (*results.begin()).first == "Model" );
    _equal( (*results.begin()).second.first[0U].second, 21./101., fTolerance );
    _equal( (*results.begin()).second.first[1U].second, 20./101., fTolerance );
    _equal( (*results.begin()).second.first[2U].second, 20./101., fTolerance );
    _equal( (*results.begin()).second.first[3U].second, 20./101., fTolerance );
    _equal( (*results.begin()).second.first[4U].second, 20./101., fTolerance );

    for (std::map<std::string,std::pair<HistogramBins,size_t> >::iterator it = results.begin();
         it != results.end(); it++)
    {
        for (HistogramBins::iterator hit = (*it).second.first.begin();
             hit != (*it).second.first.end(); hit++)
        {
            std::cout << (*it).first << "\t" << (*hit).first << "\t" << (*hit).second << "\n";
        }
    }

}




void StatisticalAnalyzer_Test::StatisticalAnalyzer_TestRegionPropertyHistogramsElementProperty2BinningBasedOnProperty1()
{
    Index por_key = model_->Database().StorageKey("porosity");
    Index sat_key = model_->Database().StorageKey("saturation water");
    size_t counter(1U);

    for (std::vector<Element<1U>*>::const_iterator eit = model_->Region("Model").ElementsBegin();
        eit != model_->Region("Model").ElementsEnd(); eit++, counter++)
    {
      (*eit)->Store(por_key, makeScalar(PLAIN, 1. - (*eit)->BaryCenter()[0U]));
    }

    for (std::vector<Node<1U>*>::const_iterator nit = model_->Region("Model").NodesBegin();
        nit != model_->Region("Model").NodesEnd(); nit++)
    {
      (*nit)->Store(sat_key, makeScalar(PLAIN, (*nit)->Coordinate()[0U]));
    }


    StatisticalAnalyzer<1U> statistical_analyzer(*model_);
    HistogramBins bins;
    std::map<std::string,std::pair<HistogramBins,size_t> > result1, result2;

    statistical_analyzer.DefineBins( 0., 1., 5, bins);
    statistical_analyzer.RegionPropertyHistogramsElementProperty2BinningBasedOnProperty1(
                         "porosity", "saturation water", bins, result1, result2, false);

    _test( (*result1.begin()).first == "Model" );
    _equal( (*result1.begin()).second.first[0U].second, 0.2, fTolerance );
    _equal( (*result1.begin()).second.first[1U].second, 0.2, fTolerance );
    _equal( (*result1.begin()).second.first[2U].second, 0.2, fTolerance );
    _equal( (*result1.begin()).second.first[3U].second, 0.2, fTolerance );
    _equal( (*result1.begin()).second.first[4U].second, 0.2, fTolerance );

    _test( (*result2.begin()).first == "Model" );
    _equal( (*result2.begin()).second.first[0U].second, 0.18, fTolerance );
    _equal( (*result2.begin()).second.first[1U].second, 0.14, fTolerance );
    _equal( (*result2.begin()).second.first[2U].second, 0.1, fTolerance );
    _equal( (*result2.begin()).second.first[3U].second, 0.06, fTolerance );
    _equal( (*result2.begin()).second.first[4U].second, 0.02, fTolerance );

}



} // end namespace csmp
