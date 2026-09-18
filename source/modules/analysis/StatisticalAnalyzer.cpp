// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "StatisticalAnalyzer.h"
#include "Node.h"
#include "InterFace.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
StatisticalAnalyzer<dim>::StatisticalAnalyzer( const Model<dim>& sg )
   : sref(sg),
     pref(sg.Database())
  {
  } // end


template<uint32_t dim>
StatisticalAnalyzer<dim>::~StatisticalAnalyzer()
 {
 }




/**

Goes through all the unique regions of the Model accumulating element volume
for elements with values within the range of the corresponding 'bin'. 
Returns a map of bins and flow volume for each region, using the
region name as identifier.

The column height is normalized by the total volume of the target region.

@todo make this method more general so that it can deal with multi-dimensional elements 
in a single region.

*/
template<uint32_t dim>
void StatisticalAnalyzer<dim>::RegionPropertyHistograms( const char* prop,
                                                         const HistogramBins& bins,
                                                         map<string,pair<HistogramBins,size_t> >& results )
const
 {
      ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( sref.UniqueRegions() == 0 ) {
           csmp_error.Note( ERROR, "StatisticalAnalyzer::RegionPropertyHistograms: ",
                             "No regions are defined. Nothing was done.");
           return;
       }
     csmp::Index prop_key = pref.StorageKey(prop);

     if ( prop_key.place != ELEMENT ) {
           csmp_error.Note( ERROR, "StatisticalAnalyzer::RegionPropertyHistograms: Area/volume normalization ",
                                     "only works for element properties." );
           return;
       }

     // vector<pair<double,double> >
     HistogramBins        result( bins.size(), make_pair(0.,0.) );
     double               volume;
     VectorVariable<dim>  vc;
     TensorVariable<dim>  ts;
     double               val;
     
     // defining upper bin limits in result vector (lowest limit is 0.0)
     typename HistogramBins::iterator rit=result.begin();
     
     for ( typename HistogramBins::const_iterator
           vit=bins.begin(); vit!=bins.end(); vit++, rit++ )
       {
          (*rit).first  = (*vit).second; 
          (*rit).second = 0.;
       }
     
     // 1. go through all groups and do the binning
     // -------------------------------------------
     for ( typename map<std::string,Region<dim> >::const_iterator
           grit=sref.UniqueRegionsBegin(); grit!=sref.UniqueRegionsEnd(); grit++ )
       {  
          double total_volume = 0.; 
          size_t  n(0);
           for ( auto it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++ )
             {
                total_volume += volume = (*it)->Volume();
                if ( prop_key.type == SCALAR ) {
                      val = (*it)->Read ( prop_key );
                   }
                else if ( prop_key.type == VECTOR ) {
                      (*it)->Read ( prop_key, vc ); 
                      val = vc.Length();
                   }
                else if ( prop_key.type == TENSOR ) {
                      (*it)->Read ( prop_key, ts );
                      val = ts.Determinant();
                   }
                else // ARRAY VARIABLE
                throw csmp::Exception( ERROR, "StatisticalAnalyzer<dim>::RegionPropertyHistograms:",
                                      "Array variables are not handled yet.");
               
                // adding the value to the corresponding column of the histogram
                size_t  i(0U);
                for ( typename HistogramBins::const_iterator
                      vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                  {
                      if (vit == bins.begin()) {
                            if ( val >= (*vit).first && val <= (*vit).second )
                              result[i].second += volume;
                        }
                      else {
                            if ( val > (*vit).first && val <= (*vit).second )
                              result[i].second += volume;
                        }
                  }
                n++;
             }

           // normalizing by area, i.e. how much of total area has this characteristic
           // ------------------------------------------------------------------------
           for ( size_t i{0U}; i<result.size(); i++ ) result[i].second /= total_volume;
            
           // 3. storing result map and zeroing vector for next region
           //---------------------------------------------------------
           results[ (*grit).first ] = pair<HistogramBins,size_t>(result,n);
         
           // zeroing out the column values
           for ( rit=result.begin(); rit!=result.end(); rit++ ) (*rit).second = 0.;
          
       } // end for all regions
     
 } // end RegionPropertyHistograms (by area/volume)    
  
  


/**

Goes through all regions of the current Model and counts property values
within the ranges defined in 'bins'. Returns a map of bins and counts for each region
using the region name as key / identifier.

@param normalize normalizes the total number of counts can be normalized by the number of counts per
region if this is set to true.

*/
template<uint32_t dim>
void StatisticalAnalyzer<dim>::RegionPropertyHistograms( const char* prop,
                                                         const HistogramBins& bins,
                                                         map<std::string,pair<HistogramBins,size_t> >& results,
                                                         bool normalize )
const
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( sref.UniqueRegions() == 0 ) {
           csmp_error.Note( ERROR, "StatisticalAnalyzer::RegionPropertyHistograms: ",
                             "No regions are defined. Nothing was done.");
           return;
       }

     csmp::Index  prop_key = pref.StorageKey(prop);
     HistogramBins  result( bins.size() );
     ScalarVariable       sc;
     VectorVariable<dim>  vc;
     TensorVariable<dim>  ts;
     double             val;
     size_t               n(1U);

     // defining upper bin limits in result vector (lowest limit is always 0.)
     HistogramBins::iterator rit=result.begin();
     for ( typename HistogramBins::const_iterator
           vit=bins.begin(); vit!=bins.end(); vit++, rit++ ) {
          (*rit).first  = (*vit).second;
          (*rit).second = 0.;
       }

     // 1. go through all groups and do the binning
     // -------------------------------------------
     for ( auto grit=sref.UniqueRegionsBegin(); grit!=sref.UniqueRegionsEnd(); grit++ )
       {
          switch ( prop_key.place )
            {
                case NODE:
                     for ( auto it=(*grit).second.NodesBegin(); it!=(*grit).second.NodesEnd(); it++ )
                       {
                          // 1.1 reading the property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          if ( prop_key.type == SCALAR ) {
                                (*it)->Read ( prop_key, sc );
                                val = sc();
                             }
                          else if ( prop_key.type == VECTOR ) {
                                (*it)->Read ( prop_key, vc );
                                val = vc.Length();
                             }
                          else { // TENSOR
                                (*it)->Read ( prop_key, ts );
                                val = ts.Determinant();
                             }
                          // 1.2 binning the value
                          // ---------------------
                          size_t i{0U};
                          for ( typename HistogramBins::const_iterator
                                vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                          {
                              if (vit == bins.begin())
                              {
                                  if ( val >= (*vit).first && val <= (*vit).second )
                                      result[i].second++;
                              }
                              else
                              {
                                  if ( val > (*vit).first && val <= (*vit).second )
                                      result[i].second++;
                              }
                          }

                       }
                     n = (*grit).second.Nodes();
                   break;
                case ELEMENT_INTEGRATION_POINT:
                     for ( auto it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++ )
                       for ( auto j{0U}; j<(*it)->IntegrationPoints(); j++ )
                         {
                            if ( prop_key.type == SCALAR ) {
                                  val = (*it)->Read( j, prop_key );
                               }
                            else if ( prop_key.type == VECTOR ) {
                                  (*it)->Read( j, prop_key, vc );
                                  val = vc.Length();
                               }
                            else { // TENSOR
                                  (*it)->Read( j, prop_key, ts );
                                  val = ts.Determinant();
                               }

                            size_t i{0U};
                            for ( typename HistogramBins::const_iterator
                                  vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                            {
                                if (vit == bins.begin())
                                {
                                    if ( val >= (*vit).first && val <= (*vit).second )
                                        result[i].second++;
                                }
                                else
                                {
                                    if ( val > (*vit).first && val <= (*vit).second )
                                        result[i].second++;
                                }
                            }

                         }
                     n = (*grit).second.IntegrationPoints();
                   break;
                case ELEMENT: 

                     for ( auto it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++ )
                       {

                          if ( prop_key.type == SCALAR ) {
                                (*it)->Read ( prop_key, sc );
                                val = sc();
                             }
                          else if ( prop_key.type == VECTOR ) {
                                (*it)->Read ( prop_key, vc );
                                val = vc.Length();
                             }
                          else { // TENSOR
                                (*it)->Read ( prop_key, ts );
                                val = ts.Determinant();
                             }
                             
                          size_t i{0U};
                          for ( typename HistogramBins::const_iterator
                                vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                          {
                              if (vit == bins.begin())
                              {
                                  if ( val >= (*vit).first && val <= (*vit).second )
                                      result[i].second++;
                              }
                              else
                              {
                                  if ( val > (*vit).first && val <= (*vit).second )
                                      result[i].second++;
                              }
                          }
                       }

                     // number of samples for total normalization
                     // -----------------------------------------
                     n = (*grit).second.Cells();
                  break;
                default:
                  csmp_error.Note( ERROR, "StatisticalAnalyzer", "property placement not handled yet");
            }
         // 2. normalization of results
         // ---------------------------
         if ( normalize )
           for ( auto it=result.begin(); it!=result.end(); it++ )
             (*it).second /= n;

         // 3. storing result map and zeroing vector for next group
         //--------------------------------------------------------
         results[ (*grit).first ] = pair<HistogramBins,size_t>(result,n);
         for ( auto it=result.begin(); it!=result.end(); it++ ) (*it).second = 0.;

     } // end for all groups

 } // end RegionPropertyHistograms





/**

Goes through all Regions of current Model and calculates the area/volume of the region with property values
within the ranges defined in 'bins'. Normalizes the computed values by the total area/volume of the Region.
Returns a map of bins and area/volume normalized values for each region
using the region name as key / identifier.

*/

template<uint32_t dim>
void StatisticalAnalyzer<dim>::RegionPropertyHistogramsElement( const char* prop,
                                                                const HistogramBins& bins,
                                                                //  region_name  up.bin.lt, value, n-samples
                                                                map<std::string,pair<HistogramBins,size_t> >& results )
const
 {
     csmp::Index  prop_key = pref.StorageKey(prop);
     HistogramBins        result( bins.size() );
     double             total_volume, volume;
     ScalarVariable       sc;
     VectorVariable<dim>  vc;
     TensorVariable<dim>  ts;
     double             val;

     // defining upper bin limits in result vector (lowest limit is always 0.)
     HistogramBins::iterator rit=result.begin();
     for ( HistogramBins::const_iterator
           vit=bins.begin(); vit!=bins.end(); vit++, rit++ ) {
          (*rit).first  = (*vit).second;
          (*rit).second = 0.;
       }

     if ( !results.empty() ) results.clear();

     // 1. go through all groups and do the binning
     // -------------------------------------------
     for ( typename map<std::string,Region<dim> >::const_iterator
           grit=sref.UniqueRegionsBegin(); grit!=sref.UniqueRegionsEnd(); grit++ )
       {
          switch ( prop_key.place )
            {
                case NODE:
                     total_volume = 0.;
                     for ( auto it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++ )
                       {
                          total_volume += fabs( volume = (*it)->Volume() );
                          // 1.1 reading the property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          if ( prop_key.type == SCALAR ) {
                                (*it)->PropertyValueAtBaryCenter( prop_key, sc );
                                val = sc();
                             }
                          else if ( prop_key.type == VECTOR ) {
                                (*it)->PropertyValueAtBaryCenter( prop_key, vc );
                                val = vc.Length();
                             }
                          else { // TENSOR
                                (*it)->PropertyValueAtBaryCenter( prop_key, ts );
                                val = ts.Determinant();
                             }
                          // 1.2 binning the value
                          // ---------------------
                          size_t i{0U};
                          for ( typename HistogramBins::const_iterator
                                vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                          {
                              if (vit == bins.begin())
                              {
                                  if ( val >= (*vit).first && val <= (*vit).second )
                                      result[i].second += fabs( volume );
                              }
                              else
                              {
                                  if ( val > (*vit).first && val <= (*vit).second )
                                      result[i].second += fabs( volume );
                              }
                          }
                       }
                     // normalizing by volume, i.e. how much of total volume has this characteristic
                     // ------------------------------------------------------------------------
                     for ( rit=result.begin(); rit !=result.end(); rit++ ) (*rit).second /= total_volume;
                   break;

                case ELEMENT_INTEGRATION_POINT:
                     total_volume = 0.;
                     for ( auto it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++ )
                       {
                          total_volume += fabs( volume = (*it)->Volume() );
                          // 1.1 reading the property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          if ( prop_key.type == SCALAR ) {
                                val = 0.;
                                for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                     val += (*it)->PropertyValueAtIntegrationPoint( prop_key, n );
                                  }
                                // averaging the ip values
                                val /= static_cast<double>((*it)->FE()->IntegrationPoints());
                             }
                          else if ( prop_key.type == VECTOR ) {
                                val = 0.;
                                for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                     (*it)->PropertyValueAtIntegrationPoint( prop_key, n, vc );
                                     val += vc.Length();
                                  }
                                val /= static_cast<double>((*it)->FE()->IntegrationPoints());
                             }
                          else { // TENSOR
                                val = 0.;
                                for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                     (*it)->PropertyValueAtIntegrationPoint( prop_key, n, ts );
                                     val += ts.Determinant();
                                  }
                                val /= static_cast<double>((*it)->FE()->IntegrationPoints());
                             }
                          // 1.2 binning the value
                          // ---------------------
                          size_t i{0U};
                          for ( typename HistogramBins::const_iterator
                                vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                          {
                              if (vit == bins.begin())
                              {
                                  if ( val >= (*vit).first && val <= (*vit).second )
                                      result[i].second += fabs( volume );
                              }
                              else
                              {
                                  if ( val > (*vit).first && val <= (*vit).second )
                                      result[i].second += fabs( volume );
                              }
                          }
                       }
                     // normalizing by volume, i.e. how much of total volume has this characteristic
                     // ------------------------------------------------------------------------
                     for ( rit=result.begin(); rit !=result.end(); rit++ ) (*rit).second /= total_volume;
                   break;

                case ELEMENT:
                     total_volume = 0.;
                     for ( auto it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++ )
                       {
                          total_volume += fabs( volume = (*it)->Volume() );
                          // 1.1 reading the property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          if ( prop_key.type == SCALAR ) {
                                val = (*it)->Read( prop_key );
                             }
                          else if ( prop_key.type == VECTOR ) {
                                (*it)->Read( prop_key, vc );
                                val = vc.Length();
                             }
                          else { // TENSOR
                                (*it)->Read( prop_key, ts );
                                val = ts.Determinant();
                             }
                          // 1.2 binning the value
                          // ---------------------
                          size_t i{0U};
                          for ( typename HistogramBins::const_iterator
                                vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                          {
                              if (vit == bins.begin())
                              {
                                  if ( val >= (*vit).first && val <= (*vit).second )
                                      result[i].second += fabs( volume );
                              }
                              else
                              {
                                  if ( val > (*vit).first && val <= (*vit).second )
                                      result[i].second += fabs( volume );
                              }
                          }
                       }
                     // normalizing by volume, i.e. how much of total volume has this characteristic
                     // ------------------------------------------------------------------------
                     for ( rit=result.begin(); rit !=result.end(); rit++ ) (*rit).second /= total_volume;
                   break;

                 default:
                   cout <<"\nStatisticalAnalyzer<dim>::RegionPropertyHistogramsElement: ";
                   cout <<" placement of property analyzed not recognized: "<< prop << endl;
            }

          // 2. storing result map and zeroing vector for next group
          //--------------------------------------------------------
          results[ (*grit).first ] = pair<HistogramBins,size_t>(result,(*grit).second.Cells());
          for ( HistogramBins::iterator
                ritt=result.begin(); ritt!=result.end(); ritt++ ) (*ritt).second = 0.;

     } // end for all regions

 } // end RegionPropertyHistogramsElement


    /**
     
     Goes through all Regions of current Model and calculates the area/volume of the region with property values
     within the ranges defined in 'bins'. Normalizes the computed values by the total area/volume of the Region.
     Returns a map of bins and area/volume normalized values for each region
     using the region name as key / identifier.
     
     */
    
    template<uint32_t dim>
    void StatisticalAnalyzer<dim>::RegionPropertyHistogramsIntegrationPoint( const char* prop,
                                                                            const HistogramBins& bins,
                                                                            //  region_name  up.bin.lt, value, n-samples
                                                                            map<std::string,pair<HistogramBins,size_t> >& results, const std::string& flow_domain )
    const
    {
        csmp::Index  prop_key = pref.StorageKey(prop);
        HistogramBins        result( bins.size() );
        double             total_volume(0.), volume, element_vol(0.);
        VectorVariable<dim>  vc;
        double             val;
        int count(0);
        
        // for tubes
        double patm(100325.);
        const csmp::Index p_key(pref.StorageKey("fluid pressure"));
        
        double max(0.), min(1e300);
        
        cout<<"RegionPropertyHistogramsIntegrationPoint: "<<prop<<endl;
        
        // defining upper bin limits in result vector (lowest limit is always 0.)
        HistogramBins::iterator rit=result.begin();
        for ( HistogramBins::const_iterator
             vit=bins.begin(); vit!=bins.end(); vit++, rit++ ) {
            (*rit).first  = (*vit).second;
            (*rit).second = 0.;
        }
        
        if ( !results.empty() ) results.clear();
        
        // 1. go through all groups and do the binning
        // -------------------------------------------
        const Region<dim>& subdomain(sref.Region(flow_domain));
        //    cout<<"velocity"<<endl;
        for ( auto it=subdomain.CellsBegin(); it!=subdomain.CellsEnd(); it++ )
        {
            // will work for BCC but not for tubes
            //total_volume += fabs( volume = (*it)->Volume() );
            // 1.1 reading the property value, taking length of vectors, and determinant of tensors
            // ------------------------------------------------------------------------------------
            val = 0.;
            element_vol = 0.;
            //volume = (*it)->Volume();
            for ( uint32_t n=0; n<(*it)->IntegrationPoints(); ++n ) {
                if ((*it)->PropertyValueAtIntegrationPoint( p_key, n) >= patm) {
                    (*it)->PropertyValueAtIntegrationPoint( prop_key, n, vc );
                    double det_J((*it)->det_J_AtIntegrationPoint(n));
                    //total_volume += fabs( volume );
                    volume = det_J * (*it)->WeightAtIntegrationPoint(n);
                    val += vc.Length() * volume;
                    element_vol += volume;
                } //if ((*it)->PropertyValueAtIntegrationPoint( p_key, n) >= patm) {
            } // for ( auto n=0; n<(*it)->IntegrationPoints(); ++n ) {
            //val /= static_cast<double>((*it)->FE()->IntegrationPoints());
            if (element_vol > 0.) {
                val /= element_vol;
                total_volume += element_vol;
                count++;
                if (val < min) min = val;
                if (val > max) max = val;
                
                // 1.2 binning the value
                // ---------------------
                size_t i{0U};
                for ( typename HistogramBins::const_iterator
                     vit=bins.begin(); vit!=bins.end(); vit++, i++ ) {
                    if (vit == bins.begin()) {
                        if ( val >= (*vit).first && val <= (*vit).second )
                            result[i].second += fabs( element_vol );
                        //result[i].second += 1.;
                    }
                    else {
                        if ( val > (*vit).first && val <= (*vit).second )
                            result[i].second += fabs( element_vol );
                        //result[i].second += 1.;
                    }
                }
            } // if (element_vol > 0.) {
            
        }
        //    cout<<endl;
        cout<<"min= "<<min<<" max= "<<max<<" count= "<<count<<" total_volume= "<<total_volume<<endl;
        // normalizing by volume, i.e. how much of total volume has this characteristic
        // ------------------------------------------------------------------------
        //for ( rit=result.begin(); rit !=result.end(); rit++ ) (*rit).second /= total_volume;
        for ( rit=result.begin(); rit !=result.end(); rit++ ) {
            //(*rit).second /= count;
            (*rit).second /= total_volume;
            cout<<(*rit).second<<endl;
        }
        cout<<endl;
        // 2. storing result map and zeroing vector for next group
        //--------------------------------------------------------
        results[ flow_domain ] = pair<HistogramBins,size_t>(result,subdomain.Cells());
        for ( auto r=result.begin(); r!=result.end(); r++ ) (*r).second = 0.;
        
    } // end RegionPropertyHistogramsIntegrationPoint





/**
    Whoever wrote this, should explain what it does !
*/
template<uint32_t dim>
void StatisticalAnalyzer<dim>::RegionPropertyHistogramsElementProperty2BinningBasedOnProperty1(
                                                       const char* prop1,
                                                       const char* prop2,
                                                       const HistogramBins& bins,
                                                       std::map<std::string,std::pair<HistogramBins,size_t> >& results1,
                                                       std::map<std::string,std::pair<HistogramBins,size_t> >& results2,
                                                       bool weighted_by_porosity ) const
 {
     csmp::Index  prop1_key(pref.StorageKey(prop1));
     csmp::Index  prop2_key(pref.StorageKey(prop2));
     csmp::Index  poro_key( (weighted_by_porosity==true) ? pref.StorageKey("porosity") : csmp::Index() );
     HistogramBins        result1( bins.size() );
     HistogramBins        result2( bins.size() );
     double             total_volume, volume;
     ScalarVariable       sc;
     VectorVariable<dim>  vc;
     TensorVariable<dim>  ts;
     double             val1(0.), val2(0.);
     double             porosity(1.);

     // defining upper bin limits in result vector (lowest limit is always 0.)
     HistogramBins::iterator rit1=result1.begin();
     HistogramBins::iterator rit2=result2.begin();
     for ( HistogramBins::const_iterator
           vit=bins.begin(); vit!=bins.end(); vit++, rit1++, rit2++ ) {
          (*rit1).first  = (*vit).second;
          (*rit1).second = 0.;
          (*rit2).first  = (*vit).second;
          (*rit2).second = 0.;
       }

     // 1. go through all groups and do the binning
     // -------------------------------------------
     for ( typename map<std::string,Region<dim> >::const_iterator
           grit=sref.UniqueRegionsBegin(); grit!=sref.UniqueRegionsEnd(); grit++ )
       {
          switch ( prop1_key.place )
            {
                case NODE:
                     total_volume = 0.;
                     for ( auto it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++ )
                       {
                          if ( weighted_by_porosity ) porosity = (*it)->Read( poro_key );
                          total_volume += fabs( volume = (*it)->Volume() ) * porosity;

                          // 1.0 reading the second property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          switch ( prop2_key.place )
                            {
                                case NODE:
                                        if ( prop2_key.type == SCALAR ) {
                                              (*it)->PropertyValueAtBaryCenter( prop2_key, sc );
                                              val2 = sc();
                                           }
                                        else if ( prop2_key.type == VECTOR ) {
                                              (*it)->PropertyValueAtBaryCenter( prop2_key, vc );
                                              val2 = vc.Length();
                                           }
                                        else { // TENSOR
                                              (*it)->PropertyValueAtBaryCenter( prop2_key, ts );
                                              val2 = ts.Determinant();
                                           }
                                   break;

                                case ELEMENT_INTEGRATION_POINT:
                                        if ( prop2_key.type == SCALAR ) {
                                              val2 = 0.;
                                              for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                                   val2 += (*it)->PropertyValueAtIntegrationPoint( prop2_key, n );
                                                }
                                              // averaging the ip values
                                              val2 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                                           }
                                        else if ( prop2_key.type == VECTOR ) {
                                              val2 = 0.;
                                              for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                                   (*it)->PropertyValueAtIntegrationPoint( prop2_key, n, vc );
                                                   val2 += vc.Length();
                                                }
                                              val2 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                                           }
                                        else { // TENSOR
                                              val2 = 0.;
                                              for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                                   (*it)->PropertyValueAtIntegrationPoint( prop2_key, n, ts );
                                                   val2 += ts.Determinant();
                                                }
                                              val2 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                                           }
                                   break;

                                case ELEMENT:
                                        if ( prop2_key.type == SCALAR ) {
                                              val2 = (*it)->Read( prop2_key );
                                           }
                                        else if ( prop2_key.type == VECTOR ) {
                                              (*it)->Read( prop2_key, vc );
                                              val2 = vc.Length();
                                           }
                                        else { // TENSOR
                                              (*it)->Read( prop2_key, ts );
                                              val2 = ts.Determinant();
                                           }
                                   break;

                                 default:
                                   cout <<"\nStatisticalAnalyzer<dim>::RegionPropertyHistogramsElement: ";
                                   cout <<" placement of property analyzed not recognized: "<< prop2 << endl;

                            } // end of switch ( prop2_key.place )

                          // 1.1 reading the first property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          if ( prop1_key.type == SCALAR ) {
                                (*it)->PropertyValueAtBaryCenter( prop1_key, sc );
                                val1 = sc();
                             }
                          else if ( prop1_key.type == VECTOR ) {
                                (*it)->PropertyValueAtBaryCenter( prop1_key, vc );
                                val1 = vc.Length();
                             }
                          else { // TENSOR
                                (*it)->PropertyValueAtBaryCenter( prop1_key, ts );
                                val1 = ts.Determinant();
                             }
                          // 1.2 binning the value
                          // ---------------------
                          size_t i{0U};
                          for ( typename HistogramBins::const_iterator
                                vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                          {
                              if (vit == bins.begin())
                              {
                                  if ( val1 >= (*vit).first && val1 <= (*vit).second )
                                  {
                                      result1[i].second += fabs( volume*porosity );
                                      result2[i].second += fabs( val2*volume*porosity );
                                  }
                              }
                              else
                              {
                                  if ( val1 > (*vit).first && val1 <= (*vit).second )
                                  {
                                      result1[i].second += fabs( volume*porosity );
                                      result2[i].second += fabs( val2*volume*porosity );
                                  }
                              }
                          }

                       }
                     // normalizing by volume, i.e. how much of total volume has this characteristic
                     // ------------------------------------------------------------------------
                     for ( rit1=result1.begin(); rit1 !=result1.end(); rit1++ ) (*rit1).second /= total_volume;
                     for ( rit2=result2.begin(); rit2 !=result2.end(); rit2++ ) (*rit2).second /= total_volume;
                   break;

                case ELEMENT_INTEGRATION_POINT:
                     total_volume = 0.;
                     for ( auto it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++ )
                       {
                          if ( weighted_by_porosity ) porosity = (*it)->Read( poro_key );
                          total_volume += fabs( volume = (*it)->Volume() ) * porosity;
                          // 1.0 reading the second property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          switch ( prop2_key.place )
                            {
                                case NODE:
                                        if ( prop2_key.type == SCALAR ) {
                                              (*it)->PropertyValueAtBaryCenter( prop2_key, sc );
                                              val2 = sc();
                                           }
                                        else if ( prop2_key.type == VECTOR ) {
                                              (*it)->PropertyValueAtBaryCenter( prop2_key, vc );
                                              val2 = vc.Length();
                                           }
                                        else { // TENSOR
                                              (*it)->PropertyValueAtBaryCenter( prop2_key, ts );
                                              val2 = ts.Determinant();
                                           }
                                   break;

                                case ELEMENT_INTEGRATION_POINT:
                                        if ( prop2_key.type == SCALAR ) {
                                              val2 = 0.;
                                              for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                                   val2 += (*it)->PropertyValueAtIntegrationPoint( prop2_key, n );
                                                }
                                              // averaging the ip values
                                              val2 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                                           }
                                        else if ( prop2_key.type == VECTOR ) {
                                              val2 = 0.;
                                              for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                                   (*it)->PropertyValueAtIntegrationPoint( prop2_key, n, vc );
                                                   val2 += vc.Length();
                                                }
                                              val2 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                                           }
                                        else { // TENSOR
                                              val2 = 0.;
                                              for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                                   (*it)->PropertyValueAtIntegrationPoint( prop2_key, n, ts );
                                                   val2 += ts.Determinant();
                                                }
                                              val2 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                                           }
                                   break;

                                case ELEMENT:
                                        if ( prop2_key.type == SCALAR ) {
                                              val2 = (*it)->Read( prop2_key );
                                           }
                                        else if ( prop2_key.type == VECTOR ) {
                                              (*it)->Read( prop2_key, vc );
                                              val2 = vc.Length();
                                           }
                                        else { // TENSOR
                                              (*it)->Read( prop2_key, ts );
                                              val2 = ts.Determinant();
                                           }
                                   break;

                                 default:
                                   cout <<"\nStatisticalAnalyzer<dim>::RegionPropertyHistogramsElement: ";
                                   cout <<" placement of property analyzed not recognized: "<< prop2 << endl;

                            } // end of switch ( prop2_key.place )

                          // 1.1 reading the property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          if ( prop1_key.type == SCALAR ) {
                                val1 = 0.;
                                for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                     val1 += (*it)->PropertyValueAtIntegrationPoint( prop1_key, n );
                                  }
                                // averaging the ip values
                                val1 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                             }
                          else if ( prop1_key.type == VECTOR ) {
                                val1 = 0.;
                                for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                     (*it)->PropertyValueAtIntegrationPoint( prop1_key, n, vc );
                                     val1 += vc.Length();
                                  }
                                val1 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                             }
                          else { // TENSOR
                                val1 = 0.;
                                for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                     (*it)->PropertyValueAtIntegrationPoint( prop1_key, n, ts );
                                     val1 += ts.Determinant();
                                  }
                                val1 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                             }

                          // 1.2 binning the value
                          // ---------------------
                          size_t i{0U};
                          for ( typename HistogramBins::const_iterator
                                vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                          {
                              if (vit == bins.begin())
                              {
                                  if ( val1 >= (*vit).first && val1 <= (*vit).second )
                                  {
                                      result1[i].second += fabs( volume*porosity );
                                      result2[i].second += fabs( val2*volume*porosity );
                                   }
                              }
                              else
                              {
                                  if ( val1 > (*vit).first && val1 <= (*vit).second )
                                  {
                                      result1[i].second += fabs( volume*porosity );
                                      result2[i].second += fabs( val2*volume*porosity );
                                  }
                              }

                          }

                       }
                     // normalizing by volume, i.e. how much of total volume has this characteristic
                     // ------------------------------------------------------------------------
                     for ( rit1=result1.begin(); rit1 !=result1.end(); rit1++ ) (*rit1).second /= total_volume;
                     for ( rit2=result2.begin(); rit2 !=result2.end(); rit2++ ) (*rit2).second /= total_volume;
                   break;

                case ELEMENT:
                     total_volume = 0.;
                     for ( auto it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++ )
                       {
                          if ( weighted_by_porosity ) porosity = (*it)->Read( poro_key );
                          total_volume += fabs( volume = (*it)->Volume() ) * porosity;

                          // 1.0 reading the second property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          switch ( prop2_key.place )
                            {
                                case NODE:
                                        if ( prop2_key.type == SCALAR ) {
                                              (*it)->PropertyValueAtBaryCenter( prop2_key, sc );
                                              val2 = sc();
                                           }
                                        else if ( prop2_key.type == VECTOR ) {
                                              (*it)->PropertyValueAtBaryCenter( prop2_key, vc );
                                              val2 = vc.Length();
                                           }
                                        else { // TENSOR
                                              (*it)->PropertyValueAtBaryCenter( prop2_key, ts );
                                              val2 = ts.Determinant();
                                           }
                                   break;

                                case ELEMENT_INTEGRATION_POINT:
                                        if ( prop2_key.type == SCALAR ) {
                                              val2 = 0.;
                                              for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                                   val2 += (*it)->PropertyValueAtIntegrationPoint( prop2_key, n );
                                                }
                                              // averaging the ip values
                                              val2 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                                           }
                                        else if ( prop2_key.type == VECTOR ) {
                                              val2 = 0.;
                                              for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                                   (*it)->PropertyValueAtIntegrationPoint( prop2_key, n, vc );
                                                   val2 += vc.Length();
                                                }
                                              val2 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                                           }
                                        else { // TENSOR
                                              val2 = 0.;
                                              for ( auto n=0U; n<(*it)->FE()->IntegrationPoints(); n++ ) {
                                                   (*it)->PropertyValueAtIntegrationPoint( prop2_key, n, ts );
                                                   val2 += ts.Determinant();
                                                }
                                              val2 /= static_cast<double>((*it)->FE()->IntegrationPoints());
                                           }
                                   break;

                                case ELEMENT:
                                        if ( prop2_key.type == SCALAR ) {
                                              val2 = (*it)->Read( prop2_key );
                                           }
                                        else if ( prop2_key.type == VECTOR ) {
                                              (*it)->Read( prop2_key, vc );
                                              val2 = vc.Length();
                                           }
                                        else { // TENSOR
                                              (*it)->Read( prop2_key, ts );
                                              val2 = ts.Determinant();
                                           }
                                   break;

                                 default:
                                   cout <<"\nStatisticalAnalyzer<dim>::RegionPropertyHistogramsElement: ";
                                   cout <<" placement of property analyzed not recognized: "<< prop2 << endl;

                            } // end of switch ( prop2_key.place )
                          // 1.1 reading the property value, taking length of vectors, and determinant of tensors
                          // ------------------------------------------------------------------------------------
                          if ( prop1_key.type == SCALAR ) {
                                val1 = (*it)->Read( prop1_key );
                             }
                          else if ( prop1_key.type == VECTOR ) {
                                (*it)->Read( prop1_key, vc );
                                val1 = vc.Length();
                             }
                          else { // TENSOR
                                (*it)->Read( prop1_key, ts );
                                val1 = ts.Determinant();
                             }

                          // 1.2 binning the value
                          // ---------------------
                          size_t i{0U};
                          for ( typename HistogramBins::const_iterator
                                vit=bins.begin(); vit!=bins.end(); vit++, i++ )
                          {
                              if (vit == bins.begin())
                              {
                                  if ( val1 >= (*vit).first && val1 <= (*vit).second )
                                  {
                                      result1[i].second += fabs( volume*porosity );
                                      result2[i].second += fabs( val2*volume*porosity );
                                  }
                              }
                              else
                              {
                                  if ( val1 > (*vit).first && val1 <= (*vit).second )
                                  {
                                      result1[i].second += fabs( volume*porosity );
                                      result2[i].second += fabs( val2*volume*porosity );
                                   }
                              }

                          }

                       }
                     // normalizing by volume, i.e. how much of total volume has this characteristic
                     // ------------------------------------------------------------------------
                     for ( rit1=result1.begin(); rit1 !=result1.end(); rit1++ ) (*rit1).second /= total_volume;
                     for ( rit2=result2.begin(); rit2 !=result2.end(); rit2++ ) (*rit2).second /= total_volume;
                   break;

                 default:
                   cout <<"\nStatisticalAnalyzer<dim>::RegionPropertyHistogramsElement: ";
                   cout <<" placement of property analyzed not recognized: "<< prop1 << endl;
            }

          // 2. storing result map and zeroing vector for next group
          //--------------------------------------------------------
          results1[ (*grit).first ] = pair<HistogramBins,size_t>(result1,(*grit).second.Cells());
          for ( HistogramBins::iterator
                rit=result1.begin(); rit!=result1.end(); rit++ ) (*rit).second = 0.;

          results2[ (*grit).first ] = pair<HistogramBins,size_t>(result2,(*grit).second.Cells());
          for ( HistogramBins::iterator
                rit=result2.begin(); rit!=result2.end(); rit++ ) (*rit).second = 0.;

     } // end for all regions

 } // end RegionPropertyHistogramsElementProperty2BinningBasedOnProperty1




template<uint32_t dim>
void StatisticalAnalyzer<dim>::WriteHistogramToTextfile( const char* fname,
                                                         const HistogramBins& hist,
                                                         size_t points )
const
 {
    char file[200];
    strcpy( file, fname );
    strcat( file, ".txt" );
    ofstream  ofs( file );

    ofs <<"StatisticalAnalyzer::WriteHistogramToTextfile: "<< fname <<", N samples: "<< points << endl;
    ofs <<"upper bound of bin,\t (normalized) counts" << endl;
    ofs.setf( ios::scientific );
    for ( HistogramBins::const_iterator it=hist.begin(); it!=hist.end(); it++ )
      ofs << (*it).first <<"\t "<< (*it).second << endl;
    ofs.close();
    cout <<"\nStatisticalAnalyzer::WriteHistogramToTextfile: file '"<< file;
    cout <<"' written successfully..." << endl;

 } // end WriteHistogram






/**  Output to Maple (old statistics package), format:

     writing: datasetname := [ Weight(1..3, 5), Weight(3..5, 10), Weight(5..7, 8) ]:
*/
template<uint32_t dim>
void StatisticalAnalyzer<dim>::WriteHistogramToMapleTextfile( const char* fname,
                                                              const char* datasetname,
                                                              const HistogramBins& bins,
                                                              const HistogramBins& hist,
                                                              size_t points,
                                                              bool log10_of_bin_data )
const
 {
    typename HistogramBins::const_iterator  it  = hist.begin();
    typename HistogramBins::const_iterator  itb = bins.begin();
    char file[200], dataset[200];
    strcpy( dataset, datasetname );
    strcpy( file, fname );
    strcat( file, ".txt" );
    ofstream  ofs( file );
    double binval1, binval2;

    ofs <<"StatisticalAnalyzer::WriteHistogramToMapleTextfile: "<< fname <<", N samples: "<< points << endl;
    ofs << endl << endl;

//cout <<"\nbins:"<< endl;
//for ( itb=bins.begin(); itb!= bins.end(); itb++ )
//  cout <<"\n"<< (*itb).first <<" - "<< (*itb).second;
//cout << endl;
//itb=bins.begin();

    // removing potential hyphens from the dataset name
    int i(0);
    while ( dataset[i++] != '\0' ) { if ( dataset[i] == '-' ) dataset[i] = '_'; }
    // writing a weighted list for maple
    // listname and assignment
    ofs << dataset <<" := [ ";

    // writing the entries
//    ofs.setf( ios::scientific );
    while ( it!=hist.end() )
      {
         if ( log10_of_bin_data ) {
              binval1 = log10( (*itb).first );
              binval2 = log10( (*itb).second );
           }
         else  {
              binval1 = (*itb).first;
              binval2 = (*itb).second;
           }
         // the weight                                         column height
         ofs <<"Weight("<< binval1 <<".."<< binval2 <<","<< (*it).second <<")";
         // the item separator
         if ( ++it != hist.end() ) ofs <<", ";
         itb++;
      }
    ofs <<" ]:"<< endl;
    ofs.close();
    cout <<"\nStatisticalAnalyzer::WriteHistogramToMapleTextfile: file '"<< file;
    cout <<"' written successfully..." << endl;

 } // end WriteHistogramToMapleTextfile





/// writing: datasetname := [ [bin_bound1+bin_bound2/2., value1], [bin_bound2+bin_bound3/2., value2],... ]:
template<uint32_t dim>
void StatisticalAnalyzer<dim>::WriteCurvePointsToMapleTextfile( const char* fname,
                                                               const char* datasetname,
                                                               const HistogramBins& bins,
                                                               const HistogramBins& hist,
                                                               size_t points,
                                                               bool log10_of_bin_data )
const
 {
    typename HistogramBins::const_iterator  it  = hist.begin();
    typename HistogramBins::const_iterator  itb = bins.begin();
    char file[200], dataset[200];
    strcpy( dataset, datasetname );
    strcpy( file, fname );
    strcat( file, ".txt" );
    ofstream  ofs( file );
    double binval1, binval2, bin_center;

    ofs <<"StatisticalAnalyzer::WriteCurvePointsToMapleTextfile: "<< fname <<", N samples: "<< points << endl;
    ofs << endl << endl;

    // removing potential hyphens from the dataset name
    int i(0);
    while ( dataset[i++] != '\0' ) { if ( dataset[i] == '-' ) dataset[i] = '_'; }
    // writing a weighted list for maple
    // listname and assignment
    ofs << dataset <<" := [ ";

    // writing the entries
//    ofs.setf( ios::scientific );
    while ( it!=hist.end() )
      {
         if ( log10_of_bin_data ) {
              binval1 = log10( (*itb).first );
              binval2 = log10( (*itb).second );
           }
         else  {
              binval1 = (*itb).first;
              binval2 = (*itb).second;
           }
         // finding the bin center
         bin_center = (binval1 + binval2) / 2.;
         // the weight                  column height
         ofs <<"["<< bin_center <<","<< (*it).second <<"]";
         // the item separator
         if ( ++it != hist.end() ) ofs <<", ";
         itb++;
      }
    ofs <<" ]:"<< endl;
    ofs.close();
    cout <<"\nStatisticalAnalyzer::WriteCurvePointsToMapleTextfile: file '"<< file;
    cout <<"' written successfully..." << endl;

 } // end WriteCurvePointsToMapleTextfile





template<uint32_t dim>
void StatisticalAnalyzer<dim>::OutputRegionPropertyHistograms( const char* prop,
                                                             // bin ranges from < to <=
                                    const map<std::string,pair<HistogramBins,size_t> >& results )
const
 {
     std::string  file_name;

     // 1. Retrieving histograms from Model
     // ----------------------------------------

     typename map<std::string,pair<HistogramBins,size_t> >::const_iterator  it;

     if ( results.empty() ) {
          cout <<"\nStatisticalAnalyzer::OutputRegionPropertyHistograms: No result values yet."<< endl;
          return;
       }

     // 2. Writing histograms to a list of files
     // ----------------------------------------
     for ( it=results.begin(); it!=results.end(); it++ )
       {
           file_name  = (*it).first;
           file_name += "-";
           file_name += prop;
           for ( string::iterator t=file_name.begin(); t!=file_name.end(); t++ ) if ( *t == ' ' ) *t = '_';
           WriteHistogramToTextfile( file_name.c_str(), (*it).second.first, (*it).second.second );
       }
 } // end OutputRegionPropertyHistograms





/**

In Maple paste the data into the list as shown below:

 with(stats):
> with(stats[statplots]):
> data1:=[ Weight(1..3, 5), Weight(3..5, 10), Weight(5..7, 8)]:
> histogram(data1, color=cyan);
*/
template<uint32_t dim>
void StatisticalAnalyzer<dim>::OutputRegionPropertyHistogramsMaple( const char* prop,
                                                                    const HistogramBins& bins,
                                                                    // bin ranges from < to <=
                                                                    const map<std::string,pair<HistogramBins,size_t> >& results,
                                                                    bool log10_of_bin_values )
const
 {
     std::string  file_name;

     // 1. Retrieving histograms from Model
     // ----------------------------------------
     if ( results.empty() ) {
          cout <<"\nStatisticalAnalyzer::OutputRegionPropertyHistogramsMaple: No result values yet."<< endl;
          return;
       }

     // 2. Writing histograms to a list of files
     // ----------------------------------------
     for ( map<std::string,pair<HistogramBins,size_t> >::const_iterator
           it=results.begin(); it!=results.end(); it++ )
       {
           // name of the group
           file_name  = (*it).first;
           file_name += "_";
           // name of the property
           file_name += prop;
           for ( string::iterator t=file_name.begin(); t!=file_name.end(); t++ ) if ( *t == ' ' ) *t = '_';
           //                             filename            datasetname
           WriteHistogramToMapleTextfile( file_name.c_str(), file_name.c_str(), bins,
                                         (*it).second.first, (*it).second.second,
                                          log10_of_bin_values );
       }

 } // end OutputRegionPropertyHistogramsMaple








/// outputs a list of bin-center, column-height pairs to plot histogram-like curves
template<uint32_t dim>
void StatisticalAnalyzer<dim>::OutputRegionPropertyAbundancePolygonsMaple( const char* prop,
                                    const HistogramBins& bins,
                                                             // bin ranges from < to <=
                                    const map<std::string,pair<HistogramBins,size_t> >& results,
                                    bool log10_of_bin_values )
const
 {
     std::string  file_name;

     if ( results.empty() ) {
          cout <<"\nStatisticalAnalyzer::OutputRegionPropertyAbundancePolygonsMaple: No result values yet."<< endl;
          return;
       }

     // 1. Retrieving histograms from Model
     // ----------------------------------------
     for ( auto it=results.begin(); it!=results.end(); it++ )
       {
           // name of the group
           file_name  = (*it).first;
           file_name += "_";
           // name of the property
           file_name += prop;

           // 2. Writing histograms to a list of files
           // ----------------------------------------
           for ( string::iterator t=file_name.begin(); t!=file_name.end(); t++ ) if ( *t == ' ' ) *t = '_';
             //                               filename            datasetname
             WriteCurvePointsToMapleTextfile( file_name.c_str(), file_name.c_str(), bins,
                                             (*it).second.first, (*it).second.second,
                                              log10_of_bin_values );
       }

 } // end OutputRegionPropertyAbundancePolygonsMaple






/// outputs a list of bin-center, column-height pairs to plot histogram-like curves
template<uint32_t dim>
void StatisticalAnalyzer<dim>::OutputRegionPropertyAbundancePolygonsMaple( const char* prop,
                                    const char* file_name_prefix,
                                    const HistogramBins& bins,
                                                             // bin ranges from < to <=
                                    const map<std::string,pair<HistogramBins,size_t> >& results,
                                    bool log10_of_bin_values )
const
 {
     std::string  file_name;

     if ( results.empty() ) {
          cout <<"\nStatisticalAnalyzer::OutputRegionPropertyAbundancePolygonsMaple: No result values yet."<< endl;
          return;
       }

     // 1. Retrieving histograms from Model
     // ----------------------------------------
     for ( auto it=results.begin(); it!=results.end(); it++ )
       {
           file_name  = file_name_prefix;
           file_name +="_";
           // name of the group
           file_name += (*it).first;
           file_name += "_";
           // name of the property
           file_name += prop;

           // 2. Writing histograms to a list of files
           // ----------------------------------------
           for ( string::iterator t=file_name.begin(); t!=file_name.end(); t++ ) if ( *t == ' ' ) *t = '_';
             //                               filename            datasetname
             WriteCurvePointsToMapleTextfile( file_name.c_str(), file_name.c_str(), bins,
                                             (*it).second.first, (*it).second.second,
                                              log10_of_bin_values );
       }

 } // end OutputRegionPropertyAbundancePolygonsMaple (with name prefix)







/**

Function defines the bins for the histogram on the basis of the user
input as specified in file which is read by this function. The file
should contain the following information:

title
number of bin-dilimiter values (integer)
bin bound values (floating point)
EOF

Example
'myfile.bins'
8
1.0e-5
2.0e-5
4.0e-5
6.0e-5
8.0e-5
1.0e-4
2.0e-4
6.0e-4
EOF
*/
template<uint32_t dim>
void StatisticalAnalyzer<dim>::DefineBins( const char* bin_file,
                                           HistogramBins& data ) const
 {
    data.erase( data.begin(), data.end() );
    double  first_val, val;

    ifstream ifs( bin_file );

    if ( !ifs.is_open() ) {
         throw csmp::Exception( ERROR, "StatisticalAnalyzer<dim>::DefineBins",
                                bin_file, "Histogram definition file could not be openend. Create it if missing.");
         return;
      }

    // file header
    char title[250];
    ifs.getline( title, 250 );
    cout <<"\nStatisticalAnalyzer::DefineBins: File header: "<< title << endl;

    // number of bin boundaries to read
    uint32_t bins;
    ifs >> bins;

    if ( bins > 100 )
      cout <<"\nStatisticalAnalyzer::DefineBins: '"<< bin_file
           <<"' contains more than 100 bins ? n="<< bins << endl;

    data.reserve(bins);

    // bin boundaries
    for ( auto n=0U; n<bins; n++ ) {
         if ( ifs.eof() ) break;
         ifs >> val;
         if ( (fabs(val) < 1.0e-30 && val != 0.0) || fabs(val) > 1.0e+30 ) {
              cout <<"\nStatisticalAnalyzer::DefineBins: Encountered anomalous bin value: "<< val << endl;
              cout <<" at position: "<< n <<" in file."<< endl;
           }
         if ( n == 0 ) first_val = val;
         else data.push_back( make_pair(first_val,val) );
         first_val = val;
      }

    // feedback
    cout <<"\nStatisticalAnalyzer::DefineBins: defined the bins (ranges): "<< endl;
    cout.setf( ios::scientific );
    typename HistogramBins::const_iterator  it;
    int i=1;
    for ( it=data.begin(); it!=data.end(); it++ )
     cout << i++ <<": "<< (*it).first <<" to "<< (*it).second << endl;
    cout << endl;

 } // end defineBins






template<uint32_t dim>
void StatisticalAnalyzer<dim>::DefineBins( HistogramBins& bins, double first_val, ... ) const
 {
    //bins.erase( bins.begin(), bins.end() );
    bins.clear();
    double  val;

    // unamed argument pointer
    va_list  ap;

    // point ap to last named argument in function
    va_start( ap, first_val );
    size_t isize(0);
    while ( (val=va_arg(ap,double)) != 0. )
      {
         isize++;
      }
//cout<<"StatisticalAnalyzer test isize"<<isize<<endl;
    bins.reserve( isize );
    // go through un-named argument list
    va_start( ap, first_val );
    while ( (val=va_arg(ap,double)) != 0. )
      {
//cout<<"Val ."<<val<<endl;
         // put here the minimum and maximum values allowed in Property database
         if ( log10(val) < -30.0 || log10(val) > 30.0 ) break;         
         bins.push_back( make_pair(first_val,val) );
         first_val = val;
      }
    va_end( ap );

    cout <<"\nStatisticalAnalyzer::DefineBins: defined the bins (ranges): "<< endl;
    cout.setf( ios::scientific );
    int32_t i(1);
    for ( typename HistogramBins::const_iterator it=bins.begin(); it!=bins.end(); it++ )
     cout << i++ <<": "<< (*it).first <<" to "<< (*it).second << endl;
    cout << endl;

 } // end defineBins


template<uint32_t dim>
void StatisticalAnalyzer<dim>::DefineBins( const double minimum, 
                                           const double maximum, 
                                           const size_t number_of_bins,
                                           HistogramBins& data ) const
 {
    
    if ( ( minimum >= maximum ) || ( number_of_bins <= 0 ) )
      {
         throw csmp::Exception( ERROR, "StatisticalAnalyzer<dim>::DefineBins",
                                    "Unacceptable range or number of bins");
         return;
      }
      
    data.erase( data.begin(), data.end() );
    data.reserve( number_of_bins + 1 );
    double bin_size ( ( maximum - minimum ) / static_cast<double>(number_of_bins) );

    double  first_val, val;


    for ( size_t n=0U; n < (number_of_bins + 1); n++ ) {
         val = minimum + n * bin_size;
         if ( n == 0 ) first_val = val;
         else data.push_back( make_pair(first_val,val) );
         first_val = val;
      }

    // feedback
    cout <<"\nStatisticalAnalyzer::DefineBins: defined the bins (ranges): "<< endl;
    cout.setf( ios::scientific );
    typename HistogramBins::const_iterator  it;
    int i=1;
    for ( it=data.begin(); it!=data.end(); it++ )
     cout << i++ <<": "<< (*it).first <<" to "<< (*it).second << endl;
    cout << endl;

 } // end defineBins


template class StatisticalAnalyzer<1U>;
template class StatisticalAnalyzer<2U>;
template class StatisticalAnalyzer<3U>;

} // csmp
