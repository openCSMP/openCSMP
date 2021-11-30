#include "TextInterface.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Matrix.h"
#include "ArrayVariable.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

TextInterface::TextInterface()
 {
    cout <<"\nTextInterface(constructor): Building module for text I/O."<< endl;
 }




/**
 
Saves the values of a physical variable into a text output file.
The first column in the output file gives you the unique ID number of 
each node, edge, face, interface or element depending where the specified 
variable has been defined. For node variables the coordinates of the nodes 
are supplied in the second and third column of the output file. For
face, interface or element variables coordinates are also supplied 
in the second and third column. These coordinates refer to the centerpoint 
(barycentre) of each element. For integration points coordinates
are given for each point. 

@section arguments Input Arguments 

Outputs into a file named by the first argument string plus extension 
".txt" the variable identified by the second string argument. An optional 
third argument allows you to append the timestep to the filename. 
This timestep information must be an integer. 
  */
template<size_t dim>
void TextInterface::OutputDataAsTextColumns( const Model<dim>& sg,
                                             const char* fname, const char* var ) const
 {
    csmp::Index  prop_key = sg.Database().StorageKey(var);
    char         file_name[200]; strcpy( file_name, fname );
    FILE*        fp(0);

     // 0. Opening file but only if the property is placed on model subdomain
     if (  prop_key.place == MODEL or prop_key.place == REGION or
           prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
         if ((fp = fopen( file_name,"wt")) == NULL ) {
              strcat( file_name, ": could not be created !"); 
              throw csmp::Exception( ERROR, "TextInterface::OutputDataAsTextColumns:", file_name );
              return;
           }
         fprintf( fp, "%s  data (SI units)\n", var );
       }
    
    // 1. Dealing with Region properties
    if (  prop_key.place == MODEL ) {
           fprintf( fp, "Model\tproperty value\n" );
           switch (prop_key.type)
             {
                case SCALAR:
                    fprintf( fp, "%E\n", sg.Read( prop_key ) );
                  break;
                case VECTOR: {
                     VectorVariable<dim>  vc;
                     sg.Read( prop_key, vc );
                     for ( size_t j=0U; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                     fprintf( fp, "\n");
                    }
                  break;
                case TENSOR: { 
                     TensorVariable<dim>  ts;
                     sg.Read( prop_key, ts );
                     for ( size_t j=0U; j<dim; j++ ) 
                       for ( size_t k=0U; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                     fprintf( fp, "\n");
                  }
                  break;
                case ARRAY: {
                     ArrayVariable ary;
                     sg.Read( prop_key, ary );
                     for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                     fprintf( fp, "\n");
                  }
                  break;
                case FLAGGEDARRAY: {
                     FlaggedArrayVariable ary;
                     sg.Read( prop_key, ary );
                     for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                     fprintf( fp, "\n");
                  }
                  break;
                default:
                  throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                "variable type not supported yet.");
           }
       }
   
    // 2. Dealing with Region properties
    else if (  prop_key.place == REGION ) {
          fprintf( fp, "Region\tproperty value\n" );
          // unique regions
          for ( typename map<string,csmp::Region<dim> >::const_iterator
                it=sg.UniqueRegionsBegin(); it!=sg.UniqueRegionsEnd(); it++ )
           {
               fprintf( fp, "%s \t", (*it).first.c_str() );
               switch (prop_key.type)
                 {
                    case SCALAR:
                        fprintf( fp, "%E\n", (*it).second.Read( prop_key ) );
                      break;
                    case VECTOR: {
                         VectorVariable<dim>  vc;
                         (*it).second.Read( prop_key, vc );
                         for ( size_t j=0U; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                         fprintf( fp, "\n");
                        }
                      break;
                    case TENSOR: { 
                         TensorVariable<dim>  ts;
                         (*it).second.Read( prop_key, ts );
                         for ( size_t j=0U; j<dim; j++ ) 
                           for ( size_t k=0U; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case ARRAY: {
                         ArrayVariable ary;
                         (*it).second.Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case FLAGGEDARRAY: {
                         FlaggedArrayVariable ary;
                         (*it).second.Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    default:
                      throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                    "variable type not supported yet.");
                 }
           }
          // other regions
          for ( typename map<string,csmp::Region<dim> >::const_iterator
                it=sg.RegionsBegin(); it!=sg.RegionsEnd(); it++ )
           {
               fprintf( fp, "%s \t", (*it).first.c_str() );
               switch (prop_key.type)
                 {
                    case SCALAR:
                        fprintf( fp, "%E\n", (*it).second.Read( prop_key ) );
                      break;
                    case VECTOR: {
                         VectorVariable<dim>  vc;
                         (*it).second.Read( prop_key, vc );
                         for ( size_t j=0U; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                         fprintf( fp, "\n");
                        }
                      break;
                    case TENSOR: { 
                         TensorVariable<dim>  ts;
                         (*it).second.Read( prop_key, ts );
                         for ( size_t j=0U; j<dim; j++ ) 
                           for ( size_t k=0U; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case ARRAY: {
                         ArrayVariable ary;
                         (*it).second.Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case FLAGGEDARRAY: {
                         FlaggedArrayVariable ary;
                         (*it).second.Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    default:
                      throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                    "variable type not supported yet.");
                 }
           }
          fclose(fp);
          return;  
       }

    // 3. Dealing with Boundary and SplitBoundary properties
    if (  prop_key.place == BOUNDARY ) {
          fprintf( fp, "Boundary\tproperty value\n" );
          // boundaries
          for ( typename map<std::string,csmp::Boundary<dim> >::const_iterator
                it=sg.BoundariesBegin(); it!=sg.BoundariesEnd(); it++ )
           {
               string  bname( (*it).first );
               fprintf( fp, "%s \t", bname.c_str() );
               switch (prop_key.type)
                 {
                    case SCALAR:
                        fprintf( fp, "%E\n", (*it).second.Read( prop_key ) );
                      break;
                    case VECTOR: {
                         VectorVariable<dim>  vc;
                         (*it).second.Read( prop_key, vc );
                         for ( size_t j=0U; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                         fprintf( fp, "\n");
                        }
                      break;
                    case TENSOR: { 
                         TensorVariable<dim>  ts;
                         (*it).second.Read( prop_key, ts );
                         for ( size_t j=0U; j<dim; j++ ) 
                           for ( size_t k=0U; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case ARRAY: {
                         ArrayVariable ary;
                         (*it).second.Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case FLAGGEDARRAY: {
                         FlaggedArrayVariable ary;
                         (*it).second.Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    default:
                      throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                    "variable type not supported yet.");
                 }
           }
          // split boundaries
          fprintf( fp, "SplitBoundary\tproperty value\n" );
          for ( typename map<std::string,csmp::SplitBoundary<dim> >::const_iterator
                it=sg.SplitBoundariesBegin(); it!=sg.SplitBoundariesEnd(); it++ )
           {
               string  bname( (*it).first );
               fprintf( fp, "%s \t", bname.c_str() );
               switch (prop_key.type)
                 {
                    case SCALAR:
                        fprintf( fp, "%E\n", (*it).second.Read( prop_key ) );
                      break;
                    case VECTOR: {
                         VectorVariable<dim>  vc;
                         (*it).second.Read( prop_key, vc );
                         for ( size_t j=0U; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                         fprintf( fp, "\n");
                        }
                      break;
                    case TENSOR: { 
                         TensorVariable<dim>  ts;
                         (*it).second.Read( prop_key, ts );
                         for ( size_t j=0U; j<dim; j++ ) 
                           for ( size_t k=0U; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case ARRAY: {
                         ArrayVariable ary;
                         (*it).second.Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case FLAGGEDARRAY: {
                         FlaggedArrayVariable ary;
                         (*it).second.Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    default:
                      throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                    "variable type not supported yet.");
                 }
           }
          fclose(fp);
          return;  
     }
 
    // 4. Output of any other properties
    OutputDataAsTextColumns( "Model", sg, fname, var );
 
 } // end





/**
    Output of data from one specific region of interest.
*/
template<size_t dim>
void TextInterface::OutputDataAsTextColumns( const char* region,
                                             const Model<dim>& sg,
                                             const char* fname, const char* s ) const
 {
     const Region<dim>  super_group(sg.Region(region));
     char  file_name[200];
     strcpy( file_name, fname );
     strcat( file_name, ".txt");
     FILE                                           *fp;
     size_t                                          j, k;
     csmp::Index                                     prop_key = sg.Database().StorageKey(s);
     ScalarVariable                                  sc;
     VectorVariable<dim>                             vc, center;
     TensorVariable<dim>                             ts;
     vector<double>                                xyz;

     // 2. Getting the file ready
     if ((fp = fopen ( file_name,"wt")) == NULL )
       {
          strcat( file_name, ": could not be created !"); 
          throw csmp::Exception( ERROR, "TextInterface::OutputDataAsTextColumns", file_name );
          return;
       }

     // 3. Output to file
     fprintf( fp, "%s  data (SI units)\n", s );

     switch( prop_key.place ) 
      {
        case NODE: { 
            if ( dim == 1U ) fprintf( fp, "Node\tX\t %s\n", s );
            else if ( dim == 2U ) fprintf( fp, "Node\tX\tY\t %s\n", s );
            else if ( dim == 3U ) fprintf( fp, "Node\tX\tY\tZ\t %s\n", s );
            for ( auto nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
              {
                fprintf( fp, "%E\t", (*nit)->x() );
                if ( dim != 1U ) fprintf( fp, "%E\t", (*nit)->y() );
                if ( dim == 3U ) fprintf( fp, "%E\t", (*nit)->z() );
                switch (prop_key.type)
                  {
                     case SCALAR: 
                          (*nit)->Read( prop_key, sc ); 
                          fprintf( fp, "%E\n", sc() );
                       break;
                     case VECTOR: 
                          (*nit)->Read( prop_key, vc );
                          for ( j=0; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                          fprintf( fp, "\n" );
                       break;
                     case TENSOR: 
                          (*nit)->Read( prop_key, ts );
                          for ( j=0; j<dim; j++ ) 
                            for ( k=0; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                          fprintf( fp, "\n" );
                      break;
                    case ARRAY: {
                         ArrayVariable ary;
                         (*nit)->Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case FLAGGEDARRAY: {
                         FlaggedArrayVariable ary;
                         (*nit)->Read( prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    default:
                      throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                    "variable type not supported yet.");
                  }
              }
            break;
          }
        case ELEMENT_INTEGRATION_POINT: { 
            if ( dim == 1U ) fprintf( fp, "Element ID, IntegrationPoint\tX\t %s \n", s );
            else if ( dim == 2U ) fprintf( fp, "Element ID, IntegrationPoint\tX\tY\t %s \n", s );
            else if ( dim == 3U ) fprintf( fp, "Element ID, IntegrationPoint\tX\tY\tZ\t %s \n", s );
            for ( auto eit=super_group.ElementsBegin(); eit!=super_group.ElementsEnd(); eit++ )
              for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
              {
                // printing the element id first
                fprintf( fp, "%u\t", static_cast<uint32_t>((*eit)->Idx()) );
                
                // printing the constraint point coordinates
                Point<dim> xyz((*eit)->IntegrationPoint(i));
                fprintf( fp, "%E\t", xyz[0] );
                if ( dim != 1U ) fprintf( fp, "%E\t", xyz[1] );
                if ( dim == 3U ) fprintf( fp, "%E\t", xyz[2] );

                switch (prop_key.type)
                  {
                     case SCALAR: 
                          (*eit)->Read( i, prop_key, sc ); 
                          fprintf( fp, "%E\n", sc() );
                       break;
                     case VECTOR: 
                          (*eit)->Read( i, prop_key, vc );
                          for ( j=0; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                          fprintf( fp, "\n" );
                       break;
                     case TENSOR: 
                          (*eit)->Read( i, prop_key, ts );
                          for ( j=0; j<dim; j++ ) 
                            for ( k=0; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                          fprintf( fp, "\n" );
                      break;
                    case ARRAY: {
                         ArrayVariable ary;
                         (*eit)->Read( i, prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    case FLAGGEDARRAY: {
                         FlaggedArrayVariable ary;
                         (*eit)->Read( i, prop_key, ary );
                         for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                         fprintf( fp, "\n");
                      }
                      break;
                    default:
                      throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                    "variable type not supported yet.");
                  }
              }
             break; 
          }
        case SECTOR_INTEGRATION_POINT: {
            if ( dim == 1U ) fprintf( fp, "Element ID, SectorIntegrationPoint\tX\t %s \n", s );
            else if ( dim == 2U ) fprintf( fp, "Element ID, SectorIntegrationPoint\tX\tY\t %s \n", s );
            else if ( dim == 3U ) fprintf( fp, "Element ID, SectorIntegrationPoint\tX\tY\tZ\t %s \n", s );
            for ( auto eit=super_group.ElementsBegin(); eit!=super_group.ElementsEnd(); eit++ ) {
                assert( (*eit)->FV() != NULL );
                // node numbering is equivalent to sector numbering
                for ( size_t i=0U; i<(*eit)->Nodes(); ++i )
                {
                  // printing the element id first
                  fprintf( fp, "%u\t", static_cast<uint32_t>((*eit)->Idx()) );
                  
                  // printing the sector integration point coordinates in global coordinates
                  for ( size_t j=1U; j<(*eit)->FV()->SectorPoints(i); ++j ) {
                       fprintf( fp, " | " );
                       Point<dim> xyz((*eit)->RstToXYZ((*eit)->FV()->SectorPoint(i,j)));
                       fprintf( fp, "%E\t", xyz[0] );
                       if ( dim != 1U ) fprintf( fp, "%E\t", xyz[1] );
                       if ( dim == 3U ) fprintf( fp, "%E\t", xyz[2] );
                    }
                  
                  switch (prop_key.type)
                    {
                       case SCALAR:
                             for ( size_t j=0U; j<(*eit)->FV()->SectorPoints(i); ++j )
                              (*eit)->Read( i, prop_key, sc );
                            fprintf( fp, "%E\n", sc() );
                         break;
                       case VECTOR: 
                            for ( size_t j=0U; j<(*eit)->FV()->SectorPoints(i); ++j ) {
                                 (*eit)->Read( i, j, prop_key, vc );
                                 fprintf( fp, " | " );
                                 for ( k=0; k<dim; ++k ) fprintf( fp, "%E\t", vc(k) );
                              }
                            fprintf( fp, "|\n" );
                         break;
                       case TENSOR: 
                            for ( size_t j=0U; j<(*eit)->FV()->SectorPoints(i); ++j ) {
                                 (*eit)->Read( i, j, prop_key, ts );
                                 fprintf( fp, " | " );
                                 for ( k=0; k<dim; ++k )
                                   for ( size_t l=0; l<dim; ++l ) fprintf( fp, "%E\t", ts(k,l) );
                              }
                            fprintf( fp, "|\n" );
                        break;
                      case ARRAY: {
                            ArrayVariable ary;
                            for ( size_t j=0U; j<(*eit)->FV()->SectorPoints(i); ++j ) {
                                 (*eit)->Read( i, j, prop_key, ary );
                                 fprintf( fp, " | " );
                                 for ( k=0U; k<ary.Size(); ++k ) fprintf( fp, "%E\t", ary(k) );
                              }
                           fprintf( fp, " |\n");
                        }
                        break;
                      case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            for ( size_t j=0U; j<(*eit)->FV()->SectorPoints(i); ++j ) {
                                 (*eit)->Read( i, j, prop_key, ary );
                                 fprintf( fp, "|" );
                                 for ( k=0U; k<ary.Size(); ++k ) fprintf( fp, "%E\t", ary(k) );
                              }
                           fprintf( fp, " |\n");
                        }
                        break;
                      default:
                        throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                      "variable type not supported yet.");
                    }
                }
               break; 
            } // end for loop
          }

        case FACET_INTEGRATION_POINT: {
            if ( dim == 1U ) fprintf( fp, "Element ID, SectorIntegrationPoint\tX\t %s \n", s );
            else if ( dim == 2U ) fprintf( fp, "Element ID, SectorIntegrationPoint\tX\tY\t %s \n", s );
            else if ( dim == 3U ) fprintf( fp, "Element ID, SectorIntegrationPoint\tX\tY\tZ\t %s \n", s );
            for ( auto eit=super_group.ElementsBegin(); eit!=super_group.ElementsEnd(); eit++ ) {
                assert( (*eit)->FV() != NULL );
                // node numbering is equivalent to sector numbering
                for ( size_t i=0U; i<(*eit)->FV()->Facets(); ++i )
                {
                  // printing the element id first
                  fprintf( fp, "%u\t", static_cast<uint32_t>((*eit)->Idx()) );
                  
                  // printing the facet integration point locations in global coordinates
                  for ( size_t j=1U; j<(*eit)->FV()->FacetPoints(i); ++j ) {
                       fprintf( fp, " | " );
                       Point<dim> xyz((*eit)->RstToXYZ((*eit)->FV()->FacetPoint(i,j)));
                       fprintf( fp, "%E\t", xyz[0] );
                       if ( dim != 1U ) fprintf( fp, "%E\t", xyz[1] );
                       if ( dim == 3U ) fprintf( fp, "%E\t", xyz[2] );
                    }

                  switch (prop_key.type)
                    {
                       case SCALAR:
                             for ( size_t j=0U; j<(*eit)->FV()->FacetPoints(i); ++j )
                              (*eit)->Read( i, prop_key, sc );
                            fprintf( fp, "%E\n", sc() );
                         break;
                       case VECTOR: 
                            for ( size_t j=0U; j<(*eit)->FV()->FacetPoints(i); ++j ) {
                                 (*eit)->Read( i, j, prop_key, vc );
                                 fprintf( fp, " | " );
                                 for ( k=0; k<dim; ++k ) fprintf( fp, "%E\t", vc(k) );
                              }
                            fprintf( fp, "|\n" );
                         break;
                       case TENSOR: 
                            for ( size_t j=0U; j<(*eit)->FV()->FacetPoints(i); ++j ) {
                                 (*eit)->Read( i, j, prop_key, ts );
                                 fprintf( fp, " | " );
                                 for ( k=0; k<dim; ++k )
                                   for ( size_t l=0; l<dim; ++l ) fprintf( fp, "%E\t", ts(k,l) );
                              }
                            fprintf( fp, " |\n" );
                        break;
                      case ARRAY: {
                            ArrayVariable ary;
                            for ( size_t j=0U; j<(*eit)->FV()->FacetPoints(i); ++j ) {
                                 (*eit)->Read( i, j, prop_key, ary );
                                 fprintf( fp, " | " );
                                 for ( k=0U; k<ary.Size(); ++k ) fprintf( fp, "%E\t", ary(k) );
                              }
                           fprintf( fp, "|\n");
                        }
                        break;
                      case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            for ( size_t j=0U; j<(*eit)->FV()->FacetPoints(i); ++j ) {
                                 (*eit)->Read( i, j, prop_key, ary );
                                 fprintf( fp, " | " );
                                 for ( k=0U; k<ary.Size(); ++k ) fprintf( fp, "%E\t", ary(k) );
                              }
                           fprintf( fp, " |\n");
                        }
                        break;
                      default:
                        throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                      "variable type not supported yet.");
                    }
                }
               break; 
            } // end for loop
          }
          
        case ELEMENT: { 
             if ( dim == 1U ) fprintf( fp, "Element\tX\t %s \n", s );
             else if ( dim == 2U ) fprintf( fp, "Element\tX\tY\t %s \n", s );
             else if ( dim == 3U ) fprintf( fp, "Element\tX\tY\tZ\t %s \n", s );
             for ( auto eit=super_group.ElementsBegin(); eit!=super_group.ElementsEnd(); eit++ )
              {
                 Point<dim>  center((*eit)->BaryCenter());
                 fprintf( fp, "%E\t", center[0] );
                 if ( dim != 1U ) fprintf( fp, "%E\t", center[1] );
                 if ( dim == 3U ) fprintf( fp, "%E\t", center[2] );

                 switch (prop_key.type)
                   {
                      case SCALAR: 
                           (*eit)->Read( prop_key, sc ); 
                           fprintf( fp, "%E\n", sc() );
                        break;
                      case VECTOR: 
                           (*eit)->Read( prop_key, vc );
                           for ( j=0; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                           fprintf( fp, "\n");
                        break;
                      case TENSOR: 
                           (*eit)->Read( prop_key, ts );
                           for ( j=0; j<dim; j++ ) 
                             for ( k=0; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                           fprintf( fp, "\n");
                        break;
                      case ARRAY: {
                           ArrayVariable ary;
                           (*eit)->Read( prop_key, ary );
                           for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                           fprintf( fp, "\n");
                        }
                        break;
                      case FLAGGEDARRAY: {
                           FlaggedArrayVariable ary;
                           (*eit)->Read( prop_key, ary );
                           for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                           fprintf( fp, "\n");
                        }
                        break;
                      default:
                        throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                      "variable type not supported yet.");
                   }
              }
             break; 
          }
        default:
             throw csmp::Exception( WARNING, "TextInterface::OutputDataAsTextColumns", 
                                             "variable placement not handled yet (Face, InterFace, and their IPs)"); 		       
      }
    fclose( fp );
    
    cout <<"\nModel::OutputDataAsTextColumns: '"<< file_name <<"' has been written successfully."<< endl;

 } // end OutputDataAsTextColumns
 






/// output of a suite of variables with same placement for a specific region; user can transform coordinates
template<size_t dim>
void TextInterface::OutputDataAsTextColumns( const Model<dim>& model,
                                             const char* file_name, const char* region,
                                             const CoordinateTransformer<dim>& xyzt,
                                             const list<string>& outvars ) const
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( !model.ContainsRegion(region) ) {
         csmp_error.notice( ERROR, "TextInterface::OutputDataAsTextColumns(multiple variables):",
                            region, " not found, nothing was done.");
         return;
      }
   
    // 1. checking that all input variables have the same placement, else no output is made
    csmp::Index  key = model.Database().StorageKey( (*outvars.begin()).c_str() );
    for ( list<string>::const_iterator it=outvars.begin(); it!=outvars.end(); it++ )
      if ( model.Database().Placement( (*it).c_str() ) != key.place ) {
           csmp_error.notice( ERROR, "TextInterface::OutputDataAsTextColumns(multiple variables):",
                             (*it).c_str(), "has different placement than other variables; this cannot be handled.");
           return;
        }

    if ( key.place != NODE and key.place != ELEMENT ) {
         csmp_error.notice( ERROR, "TextInterface::OutputDataAsTextColumns(multiple variables):",
                           "method works only for NODE and ELEMENT variable placements; nothing was done.");
         return;
      }

     const Region<dim> output_region(model.Region(region));

     // 2. Getting the output file ready
     string  outfile_prefix(file_name);
     outfile_prefix +="_";
     outfile_prefix += region;
     outfile_prefix +="_";
     outfile_prefix += parsePlacement(key.place);
     outfile_prefix +=".txt";
     FILE*  fp(0);
     if ((fp = fopen ( outfile_prefix.c_str(),"wt")) == NULL ) {
          csmp_error.notice( ERROR, "TextInterface::OutputDataAsTextColumns",
                             outfile_prefix.c_str(), "could not be created." );
          return;
       }

     // 3. Output to file
     fprintf( fp, "%s  data (SI units).\n", outfile_prefix.c_str() );
     if ( key.place == NODE or
          (parsePlacement(key.place)).find("POINT") != string::npos ) {
          if ( dim == 1U )      fprintf( fp, "Point data: X\t" );
          else if ( dim == 2U ) fprintf( fp, "Point data: X\tY\t" );
          else if ( dim == 3U ) fprintf( fp, "Point data: X\tY\tZ\t" );
       }
     // any other properties are considered cell properties
     else {
          if ( dim == 1U )      fprintf( fp, "Cell data: X\t" );
          else if ( dim == 2U ) fprintf( fp, "Cell data: X\tY\t" );
          else if ( dim == 3U ) fprintf( fp, "Cell data: X\tY\tZ\t" );
       }
     // printing variable names, appending indices if these are vector or tensor variables
     vector<csmp::Index>  var_keys(outvars.size());
     size_t               var_count(0);
     for ( list<string>::const_iterator pit=outvars.begin(); pit!=outvars.end(); pit++ ) {
          var_keys[var_count] = model.Database().StorageKey( (*pit).c_str() );
          if ( var_keys[var_count].type == SCALAR ) fprintf( fp, "%s\t", (*pit).c_str() );
          else if ( var_keys[var_count].type == VECTOR ) {
                 for ( size_t i=0U; i<dim; i++ )
                   fprintf( fp, "%s%zu\t", (*pit).c_str(), i );
            }
          else if ( var_keys[var_count].type == TENSOR ) {
                 for ( size_t i=0U; i<dim; i++ )
                   for ( size_t j=0U; j<dim; j++ )
                     fprintf( fp, "%s%zu%zu\t", (*pit).c_str(), i, j );
            }
          else if ( var_keys[var_count].type == ARRAY or var_keys[var_count].type == FLAGGEDARRAY ) {
                 for ( size_t i=0U; i<var_keys[var_count].dataDepth; i++ )
                   fprintf( fp, "%s%zu\t", (*pit).c_str(), i );
            }
          var_count++;
       }
     fprintf( fp, "\n" );
   
     // =============================================================================
     //   writing the variables to file
     // =============================================================================
     VectorVariable<dim>  vc;
     TensorVariable<dim>  ts;
     ArrayVariable        ary;
     FlaggedArrayVariable fary;
     vector<double>     xyz;

     switch( key.place )
       {
         case NODE:
             for ( auto nit=output_region.NodesBegin(); nit!=output_region.NodesEnd(); nit++ )
               {
                  // 1. coordinates
                  Point<dim>  p((*nit)->Coordinate());
                  xyzt.Transform(p);
                  fprintf( fp, "%E\t", p[0] );
                  if ( dim != 1U ) fprintf( fp, "%E\t", p[1] );
                  if ( dim == 3U ) fprintf( fp, "%E\t", p[2] );

                  // 2. all variables placed on this node
                  for ( vector<csmp::Index>::const_iterator
                        vt=var_keys.begin(); vt!=var_keys.end(); vt++ )
                    {
                      switch ( (*vt).type )
                        {
                           case SCALAR:
                                fprintf( fp, "%E\t", (*nit)->Read(*vt) );
                             break;
                           case VECTOR: 
                                (*nit)->Read( (*vt), vc );
                                for ( size_t j=0; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                             break;
                           case TENSOR: 
                                (*nit)->Read( (*vt), ts );
                                for ( size_t j=0; j<dim; j++ )
                                  for ( size_t k=0; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                            break;
                          case ARRAY:
                               (*nit)->Read( (*vt), ary );
                               for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                            break;
                          case FLAGGEDARRAY:
                               (*nit)->Read( (*vt), fary );
                               for ( size_t j=0U; j<fary.Size(); j++ ) fprintf( fp, "%E\t", fary(j) );
                            break;
                          default:
                            throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                          "variable type not supported yet.");
                        }
                    }
                 fprintf( fp, "\n" );
               }
            break;
          
         case ELEMENT:
             for ( auto eit=output_region.ElementsBegin(); eit!=output_region.ElementsEnd(); eit++ )
               {
                   // 1. coordinate of element barycenter
                   Point<dim>  center((*eit)->BaryCenter());
                   xyzt.Transform(center);
                   fprintf( fp, "%E\t", center[0] );
                   if ( dim != 1U ) fprintf( fp, "%E\t", center[1] );
                   if ( dim == 3U ) fprintf( fp, "%E\t", center[2] );

                   // 2. all variables placed on this node
                   for ( vector<csmp::Index>::const_iterator
                         vt=var_keys.begin(); vt!=var_keys.end(); vt++ )
                     {
                        switch ((*vt).type)
                          {
                             case SCALAR:
                                  fprintf( fp, "%E\t", (*eit)->Read((*vt)) );
                               break;
                             case VECTOR:
                                  (*eit)->Read( (*vt), vc );
                                  for ( size_t j=0; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                               break;
                             case TENSOR:
                                  (*eit)->Read( (*vt), ts );
                                  for ( size_t j=0; j<dim; j++ )
                                    for ( size_t k=0; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                               break;
                             case ARRAY:
                                  (*eit)->Read( (*vt), ary );
                                  for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                               break;
                             case FLAGGEDARRAY:
                                  (*eit)->Read( (*vt), ary );
                                  for ( size_t j=0U; j<ary.Size(); j++ ) fprintf( fp, "%E\t", ary(j) );
                               break;
                             default:
                               throw csmp::Exception( ERROR, "TextInterface<dim>::OutputDataAsTextColumns:",
                                                             "variable type not supported yet.");
                           }
                     }
                 fprintf( fp, "\n" );
               }
             break;
         default:
            throw csmp::Exception( WARNING, "TextInterface::OutputDataAsTextColumns",
                                  "variable placement not handled yet (Face, InterFace, and their IPs)");
      }

    fclose( fp );
    
    cout <<"\nModel::OutputDataAsTextColumns: '"<< outfile_prefix <<"' has been written successfully."<< endl;

 } // end OutputDataAsTextColumns


template void TextInterface::OutputDataAsTextColumns( const Model<1U>&, const char*, const char*,
                                                      const CoordinateTransformer<1U>&, const list<string>& ) const;
template void TextInterface::OutputDataAsTextColumns( const Model<2U>&, const char*, const char*,
                                                      const CoordinateTransformer<2U>&, const list<string>& ) const;
template void TextInterface::OutputDataAsTextColumns( const Model<3U>&, const char*, const char*,
                                                      const CoordinateTransformer<3U>&, const list<string>& ) const;









template<size_t dim>
void TextInterface::OutputDataAsTextColumnsNumbered( const Model<dim>& sg,
                                                     const char* fname, const char* s ) const
 {
    OutputDataAsTextColumnsNumbered( "Model", sg, fname, s );
 }








template<size_t dim>
void TextInterface::OutputDataAsTextColumnsNumbered( const char* region, const Model<dim>& sg,
                                                     const char* fname, const char* s ) const
 {
     const Region<dim>  super_group(sg.Region(region));
     csmp::Index                prop_key = sg.Database().StorageKey(s);
     char  file_name[200];
     strcpy( file_name, fname );
     strcat( file_name, ".txt");
     FILE*                fp;
     size_t               j, k;
     ScalarVariable       sc;
     VectorVariable<dim>  vc;
     TensorVariable<dim>  ts;
     vector<double>     xyz;

     // 2. Getting the file ready
     if ((fp = fopen ( file_name,"wt")) == NULL ) {
          strcat( file_name, ": could not be created !"); 
          throw csmp::Exception( ERROR, "Model<dim>::OutputDataAsTextColumnsNumbered", file_name );
          return;
       }

     // 3. Output to file
     fprintf( fp, "%s  data (SI units)\n", s );

     switch( prop_key.place ) 
      {
        case NODE: { 
            if      ( dim == 1U ) fprintf( fp, "Node\tX\t %s\n", s );
            else if ( dim == 2U ) fprintf( fp, "Node\tX\tY\t %s\n", s );
            else if ( dim == 3U ) fprintf( fp, "Node\tX\tY\tZ\t %s\n", s );
            for ( auto nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
              {
                fprintf( fp, "%6.0lu\t",  (*nit)->Idx() );
                fprintf( fp, "%E\t", (*nit)->x() );
                if ( dim != 1U ) fprintf( fp, "%E\t", (*nit)->y() );
                if ( dim == 3U ) fprintf( fp, "%E\t", (*nit)->z() );
                switch (prop_key.type)
                  {
                     case SCALAR: 
                          (*nit)->Read( prop_key, sc ); 
                          fprintf( fp, "%E\n", sc() );
                       break;
                     case VECTOR: 
                          (*nit)->Read( prop_key, vc );
                          for ( j=0U; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                          fprintf( fp, "\n" );
                       break;
                     case TENSOR: 
                          (*nit)->Read( prop_key, ts );
                          for ( j=0U; j<dim; j++ ) 
                            for ( k=0U; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                          fprintf( fp, "\n" );
                     default:
                       throw csmp::Exception( WARNING, "TextInterface::OutputDataAsTextColumnsNumbered",
                                             "variable type not handled yet (Array & FixedArray)");
                  }
              }
            break;
          }
        case ELEMENT_INTEGRATION_POINT: { 
            if ( dim == 1U ) fprintf( fp, "Element ID, IntegrationPoint\tX\t %s \n", s );
            else if ( dim == 2U ) fprintf( fp, "Element ID, IntegrationPoint\tX\tY\t %s \n", s );
            else if ( dim == 3U ) fprintf( fp, "Element ID, IntegrationPoint\tX\tY\tZ\t %s \n", s );
            for ( auto eit=super_group.ElementsBegin(); eit!=super_group.ElementsEnd(); eit++ )
              for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
              {
                // printing the element id first
                fprintf( fp, "%u\t", static_cast<uint32_t>((*eit)->Idx()) );
                
                Point<dim>  xyz((*eit)->IntegrationPoint(i));

                fprintf( fp, "%u\t", static_cast<uint32_t>((*eit)->Idx()) );
                fprintf( fp, "%E\t", xyz[0] );
                if ( dim != 1U ) fprintf( fp, "%E\t", xyz[1] );
                if ( dim == 3U ) fprintf( fp, "%E\t", xyz[2] );

                switch (prop_key.type)
                  {
                     case SCALAR: 
                          (*eit)->Read( i, prop_key, sc ); 
                          fprintf( fp, "%E\n", sc() );
                       break;
                     case VECTOR: 
                          (*eit)->Read( i, prop_key, vc );
                          for ( j=0; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                          fprintf( fp, "\n" );
                       break;
                     case TENSOR: 
                          (*eit)->Read( i, prop_key, ts );
                          for ( j=0; j<dim; j++ ) 
                            for ( k=0; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                          fprintf( fp, "\n" );
                     default:
                       throw csmp::Exception( WARNING, "TextInterface::OutputDataAsTextColumnsNumbered",
                                             "variable type not handled yet (Array & FixedArray)");
                  }
              }
             break; 
          }
        case ELEMENT: { 
             if      ( dim == 1U ) fprintf( fp, "Element\tX\t %s \n", s );
             else if ( dim == 2U ) fprintf( fp, "Element\tX\tY\t %s \n", s );
             else if ( dim == 3U ) fprintf( fp, "Element\tX\tY\tZ\t %s \n", s );
             for ( auto eit=super_group.ElementsBegin(); eit!=super_group.ElementsEnd(); eit++ )
              {
                 fprintf( fp, "%u\t", static_cast<uint32_t>((*eit)->Idx()) );
                 xyz = ((*eit)->BaryCenter()).Coordinates();
                 fprintf( fp, "%E\t", xyz[0] );
                 if ( dim != 1U ) fprintf( fp, "%E\t", xyz[1] );
                 if ( dim == 3U ) fprintf( fp, "%E\t", xyz[2] );

                 switch (prop_key.type)
                   {
                      case SCALAR: 
                           (*eit)->Read( prop_key, sc ); 
                           fprintf( fp, "%E\n", sc() );
                        break;
                      case VECTOR: 
                           (*eit)->Read( prop_key, vc );
                           for ( j=0; j<dim; j++ ) fprintf( fp, "%E\t", vc(j) );
                           fprintf( fp, "\n");
                        break;
                      case TENSOR: 
                           (*eit)->Read( prop_key, ts );
                           for ( j=0; j<dim; j++ ) 
                             for ( k=0; k<dim; k++ ) fprintf( fp, "%E\t", ts(j,k) );
                           fprintf( fp, "\n");
                     default:
                       throw csmp::Exception( WARNING, "TextInterface::OutputDataAsTextColumnsNumbered",
                                             "variable type not handled yet (Array & FixedArray)");
                   }
              }
             break; 
          }
        case MODEL: { 
             fprintf( fp, "%s is a scalar (global) variable. ", s );
             throw csmp::Exception( WARNING, "TextInterface::OutputDataAsTextColumnsNumbered", 
                                             "no Model vars so far"); 		       
             break; 
         }
        default:
             throw csmp::Exception( ERROR, "TextInterface::OutputDataAsTextColumnsNumbered", 
                                           "placement of output variable could not be identified"); 		       
      }
    fclose( fp );
    
    cout <<"\nModel::OutputDataAsTextColumnsNumbered: '"<< file_name;
    cout <<"' has been written successfully."<< endl;

 } // end OutputDataAsTextColumnsNumbered





template<size_t dim>
void TextInterface::OutputDataAsTextColumns( const Model<dim>& sg,
                                             const char* file_name, const char* s, 
                                             long timestep, bool numbered ) const
 {
    char  fname[INFO_STRING], step[NAME_STRING];

    // the timestep is appended to the name of the file that is created
    strcpy( fname, file_name );

    sprintf( step, "%ld", timestep );
    strcat( fname, step );

    // outputting to text file
    if ( numbered ) OutputDataAsTextColumnsNumbered( sg, fname, s );
    else OutputDataAsTextColumns( sg, fname, s );
    
 } // end OutputDataAsTextColumns









/**
 
Outputs a model property (distributed physical variable) into a series of
text files with names corresponding to sub-regions of the model identified
as groups. The resulting filenames are extended by a timestep integer and 
the '.txt' extension. 

@section arguments Input Arguments 

The name of the output property and the timestep which denotes the current
state of the model. 

@section implementation Implementation

Uses the method OutputRegionDataAsTextColumns(). 

@section application Application

The output text-column format is useful for a variety of plotting programs 
such as Texplot, Sigma Plot, SpyGlass Transform, Excel etc. 

@section messages Messages 

The method will report an error and return without creating any output if 
no groups were defined. 
 
*/
template<size_t dim>
void TextInterface::OutputRegionsToTextFiles( const Model<dim>& sg, const char* property ) const
 {
    typename map<string,Region<dim> >::const_iterator  it = sg.RegionsBegin();
    char    file_name[INFO_STRING];
    string  prop(property);
    for ( string::iterator i=prop.begin(); i!=prop.end(); i++ ) if ( *i == ' ' ) *i = '_';

    while ( it!=sg.RegionsEnd() )
      {
         strcpy( file_name, (*it).first.CharPointer() );
         strcat( file_name, "_" );
         strcat( file_name, prop.c_str() );
         OutputDataAsTextColumns( (*it).first.c_str(), sg, file_name, property );
         it++;
      }

 } // end OutputRegionsToTextFiles








/**
 
Appends a line to a textfile. This line will contain the current model
time, the object-host ID of a physical variable (e.g., the node number of 
the node for which the variable value is output), and the variable value at 
specifically flagged host-object locations inside of a subregion of a
CSMP model identified as a group object.  

@section arguments Input Arguments 

AppendRegionToTextFileWhere() expects four arguments, the name of the 
textfile to which the output line shall be appended to, the name of the
subregion of the model, the name of the property which shall be output 
to file, and the group-internal flag of the host-object for which the 
variable shall be output (either PLAIN or BOUNDARY). 

@section implementation Implementation

Calls the Region class interface AppendToTextFileWhere().  

@section application Application

AppendRegionToTextFileWhere() allows the user to monitor property values 
in small subregions of a model over time. The resulting data can be 
visualized directly with Spyglass Plot, Excel or other programs.  

@section messages Messages 

The method will report an error and exit without completing its task, if 
the target group does not exist. 

 */
/*
void TextInterface::AppendRegionToTextFileWhere( const Model<2U>& model,
                                                 const char* fname,
                                                 const char* groupname,
                                                 const char* property,
                                                 VARIABLE_FLAG flag ) const
 {
    map<string,Region,less<string> >::const_iterator  iter;
    
    // finding the region in the region list
    iter = group_list.find(groupname);
    if ( iter != group_list.end() ) 
      (*iter).second.AppendToTextFileWhere( fname, property, flag );
                             
    else throw csmp::Exception( ERROR, "Model<dim>::AppendRegionToTextFileWhere", 
                                "region does not exist: ", groupname );

 } // end AppendRegionToTextFileWhere
*/




/**
 
Reads (text) pixel color output file from NIH Image vs 1.53 and later. 
The method determines the number of rows (m) and columns (n) in this 
input file. 

The matrix size, i. e. rows and columns.

*/
void TextInterface::SizeofPixelTextImage256( const char* fname, size_t& m, size_t& n )
 {
    ifstream  ifs;
    char      c;
    bool      state(true);

    // opening the output file 
    ifs.open ( fname, ios::in );
    if ( !ifs )
      throw csmp::Exception( FATAL_ERROR, "TextInterface::SizeofPixelTextImage256", 
                                          "Input file could not be opened");
    noskipws( ifs ); 

    // 1. count integers in the first line = number of columns
    //    the file starts with the first number
    n = 0; 
    do 
       {
          ifs.get( c );
          if ( c != ' ' && c != '\t' ) state = true;
          if ( (state == true && c == ' ') || (state == true && c == '\t') )
            {
               state = false;
               n++;
            }
      } while ( c != '\n' ); 
    n++; // count the last number 

    // 2. reading all rows of pixel data 
    m = 0; // the first row "\n' has not been read yet
    while ( ifs ) 
      {
        ifs.get( c );
        if ( c == '\n' ) m++;
      }
    ifs.close ();

    cout <<"\n\nTextInterface::SizeofPixelTextImage256:  ";
    cout << m <<" rows,  "<< n <<" columns." << endl;

 } // end SizeofPixelTextImage256
 






/** Initialises the supplied matrix from file.
 */
void  TextInterface::ReadPixelTextImage256( const char* fname, Matrix& data )
 {
    ifstream  ifs;
    double  value;
   
    // opening the output file 
    ifs.open ( fname, ios::in );
    if ( !ifs )
      throw csmp::Exception( FATAL_ERROR, "TextInterface::ReadPixelTextImage256", 
                                          "Input file could not be opened");
      
     for ( size_t i=0U; i<data.Rows(); i++ )
       for ( size_t j=0U; j<data.Cols(); j++ )
         {   
           ifs >> value;
           data(i, j) = value ;
         }

    ifs.close ();
  
 } // end ReadPixelTextImage256







/** Output string to textfile either in overwrite or append mode.
 */
void TextInterface::WriteStringToTextFile( const char* fname, 
                                           const string& str, 
                                           bool overwrite )
 {
    ofstream ofs;
    if ( overwrite == true )
      {
         ofs.open ( fname, ios::out|ios::trunc );
         if ( !ofs )
           throw csmp::Exception( FATAL_ERROR, "TextInterface::WriteStringToTextFile", 
                                        "Output file could not be opened");
         else ofs << str << endl;
      }
    else 
      {
         ofs.open ( fname, ios::out|ios::app );
         if ( !ofs )
           throw csmp::Exception( FATAL_ERROR, "TextInterface::WriteStringToTextFile", 
                                        "Output file could not be opened");
         else ofs << str << endl;
      }
      
 } // end WriteStringToTextFile 1




/** Output string to textfile either in overwrite or append mode.
 */
void TextInterface::WriteStringToTextFile( const char* fname, const char* str, bool overwrite )
 {
    ofstream ofs;
    if ( overwrite == true )
      {
         ofs.open ( fname, ios::out|ios::trunc );
         if ( !ofs )
           throw csmp::Exception( FATAL_ERROR, "TextInterface::WriteStringToTextFile", 
                                        "Output file could not be opened");
         else ofs << str << endl;
      }
    else 
      {
         ofs.open ( fname, ios::out|ios::app );
         if ( !ofs )
           throw csmp::Exception( FATAL_ERROR, "TextInterface::WriteStringToTextFile", 
                                        "Output file could not be opened");
         else ofs << str << endl;
      }
      
 } // end WriteStringToTextFile 2






/** Dump Matrix to File using scientific notation first line gives m x n.
  */
void TextInterface::WriteMatrixToTextfile( const char* fname, 
                                           DenseMatrix<DM_MIN>& mtrx )
 {
    ofstream ofs;
    size_t i, j, m = mtrx.Rows(), n = mtrx.Cols();
        
    // opening the output file 
    ofs.open ( fname, ios::out|ios::trunc );
    if ( !ofs )
      throw csmp::Exception( FATAL_ERROR, "TextInterface::WriteMatrixToTextfile", 
                                   "Output file could not be opened");
   
    // writing numbers and coords into the file as text 
    ofs.setf( ios::scientific );
    ofs << m <<"  "<< n << endl; 
    for ( i=0; i<m; i++ )
      { 
        for ( j=0; j<n; j++ ) 
          ofs << mtrx(i,j) << "\t";
        ofs << endl;
      }

    ofs.close ();

 } // end WriteTextBoundaryCondition





void  TextInterface::AppendDataToText( const char* fname, double time, long idx, double value )
 {
    ofstream  ofs;

    // opening the output file for appending 
    ofs.open ( fname, ios::out|ios::app );
    if ( !ofs )
      throw csmp::Exception( FATAL_ERROR, "TextInterface::AppendDataToText", 
                                   "Output file could not be opened");
     // writing out the timestep
     ofs << time <<"\t"<< idx <<"\t"<< value << endl;
     // writing out the properties
 }



/** Peaks into next file of file to see whether it contains the search string
of interest.  */
bool  TextInterface::IsInNextLine( ifstream& ifs, const char* search_string ) const
 {
    char text_line[256];
   
    // 1. remembering stream positions
    // -------------------------------
    streampos pos = ifs.tellg();
    ifs.get( text_line, 256 );
    ifs.seekg( pos );

    // seeking for word in text string
    if ( strstr( text_line, search_string ) == NULL ) return false;
    return true;
}


template void
TextInterface::OutputDataAsTextColumns<1U>( const Model<1U>& sg,
                                                         const char* fname, const char* s ) const;

template void
TextInterface::OutputDataAsTextColumns<2U>( const Model<2U>& sg, 
                                                         const char* fname, const char* s ) const;

template void
TextInterface::OutputDataAsTextColumns<3U>( const Model<3U>& sg,
                                                         const char* fname, const char* s ) const;

template void
TextInterface::OutputDataAsTextColumnsNumbered<1U>( const Model<1U>& sg,
                                                         const char* fname, const char* s ) const;

template void
TextInterface::OutputDataAsTextColumnsNumbered<2U>( const Model<2U>& sg, 
                                                         const char* fname, const char* s ) const;

template void
TextInterface::OutputDataAsTextColumnsNumbered<3U>( const Model<3U>& sg,
                                                         const char* fname, const char* s ) const;

template void
TextInterface::OutputDataAsTextColumns<1U>( const Model<1U>& sg, 
                                                         const char* fname, const char* s,
                                                         long timestep, bool ) const;

template void
TextInterface::OutputDataAsTextColumns<2U>( const Model<2U>& sg, 
                                                         const char* fname, const char* s,
                                                         long timestep, bool ) const;

template void
TextInterface::OutputDataAsTextColumns<3U>( const Model<3U>& sg, 
                                                         const char* fname, const char* s,
                                                         long timestep, bool  ) const;

} // csmp
