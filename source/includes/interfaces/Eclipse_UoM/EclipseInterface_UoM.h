#ifndef ECLIPSE_INTERFACE_UOM_H
#define ECLIPSE_INTERFACE_UOM_H

#include "CornerPointGrid_UoM.h"

#include "CSMP_highLevelUtilities.h"

namespace csmp {

namespace eclipse {
/// @file EclipseInterface.h - glue that creates a file-based interface between CSMP and Eclipse

/**
Full specifications of the run are in the RUNSPEC (run specification data) (first) section of the input file.
Some properties will appear only if invoked by keywords

RUNSPEC
GRID
PROPS(material and fluid properties, PVT data etc.), PHASES, UNIT SYSTEM .... grid definition
it is followed by:

REGIONS satnum, PVTNUM etc.
SOLUTION methods like IMPES
SUMMARY
SCHEDULE well management, time stepping, convergence control
OUTPUT
*/

/// Eclipse integer tags needed to interpret properties and constitutive relations applied to grid cells
enum ECLIPSE_PROPERTY
{
  ECLIPSE_FIPNUM = -5, ///< how many fluid-in-place regions the model contains
  ECLIPSE_ROCKNUM = -4, ///< number of different rocktypes assigned to regions (defined in conjuction with ROCKTAB and ROCKTABH)
  ECLIPSE_EQLNUM = -3, ///< number of PVT equilibrium regions to be equilibrated separately
  ECLIPSE_PVTNUM = -2, ///< number of regions with different PVT property data
  ECLIPSE_SATNUM = -1, ///< number of saturation functions tabulated in RUNSPEC section
  ECLIPSE_NTG = 1,  ///< value of Net-to-Gross ratio for specific simulation cell
  ECLIPSE_PORO = 2,  ///< porosity value of cell
  ECLIPSE_PERM = 3,  ///< permeability
  ECLIPSE_PRESSURE = 4,  ///< pressure
  ECLIPSE_SWAT = 5,  ///< water saturation
  ECLIPSE_SOIL = 6,  ///< oil saturation
  ECLIPSE_SGAS = 7,  ///< gas saturation
  ECLIPSE_MULTFLT = 8,  ///< multi-fault segment definition (FAULTS block)
  ECLIPSE_MULT = 9,  ///< generic multipliers for individual cells
  ECLIPSE_TRAN = 10  ///< transmissibility multiplier (RUNSPEC file)
};


/// ACTNUM active region number; in definition of regions?


/// phase definitions for RUNSPEC section (how many and which phases are present in run)
enum ECLIPSE_PHASE
{
  ECLIPSE_OIL = 0, ///<
  ECLIPSE_WATER = 1, ///<
  ECLIPSE_GAS = 2, ///<
  ECLIPSE_LIQ = 3  ///<
};



struct EclipseWellCompletion
{
  EclipseWellCompletion();
  ~EclipseWellCompletion();

  /// JC: not clear - COMPDAT connection specification?
  size_t ic_;      ///<
  size_t jc_;      ///<
  size_t k_top_;   ///<
  size_t k_bot_;   ///<
  size_t w_start_; ///<
  size_t w_end_;   ///<
  bool is_open_;   ///< connected with the reservoir or not
};



struct EclipseWell
{
  EclipseWell();
  ~EclipseWell();

  std::string well_name_;
  std::string group_name_;
  size_t i_start_;
  size_t j_start_;
  double z_bhp_;
  ECLIPSE_PHASE phase_;
  double drainage_radius_;
  size_t inflow_type_;
  bool stop_or_shut_in_; // true - stop, false - shut in
  bool cross_flow_;      // true - crossflow allowed, false - crossflow is not allowed
  std::vector<EclipseWellCompletion> well_data_;
};


struct EclipseWellPathEntry
{
  ijk cell;
  CORNER_POINT_CELL_FACE_INDEX from, to;

  EclipseWellPathEntry( const ijk& cell, CORNER_POINT_CELL_FACE_INDEX from, CORNER_POINT_CELL_FACE_INDEX to )
    : cell( cell ), from( from ), to( to )
  {
  }
};


struct EclipseWellPath
{
  std::vector<EclipseWellPathEntry> path;
};

struct EclipseFault
{
  std::vector< std::pair<ijk, CORNER_POINT_CELL_FACE_INDEX> > fault;

  void Reserve( size_t count )
  {
    fault.reserve( count );
  }

  void Add( size_t i, size_t j, size_t k, CORNER_POINT_CELL_FACE_INDEX face )
  {
    fault.emplace_back( ijk( i, j, j ), face );
  }
};



/**
Contains the CornerPointGrid.
*/
class  EclipseInterface
{
public:
  EclipseInterface();
  ~EclipseInterface();

  /// MASTER METHOD not only reads the Eclipse input files, but also creates grid objects and converts them into CSMP cells / mesh
  bool ReadFile( csmp::VSet<3U>& vset,
                 csmp::ModelTopology& model_topology,
                 const std::string& fname,
                 bool exclude_inactive_cells,
                 bool tetra_mesh );

  /// configures reader for the case where INCLUDE is set so that instructions are spread across multiple files
  bool ReadFile( std::ifstream& ifs, size_t line_length );

  /// assigns CSMP parameter specifications to the interface so that corresponding variables can be read
  void SetProperties( const std::map<int, csmp::Parameter>& );

  /// inserts wells in the Eclipse data
  void AddWell();

  /// inserts a well from its start and end points
  void AddWell( const std::string& well_name, const Point<3U>& well_start_point, const Point<3U>& well_end_point );

#if 0
  void AddWellFacePath( const std::string& well_name, const std::vector<size_t>& cell_ids );

  void AddWellFacePath( const std::string& well_name, const std::vector<size_t>& cell_ids,
                        const std::vector<std::pair<size_t, size_t> >& face_ids );

  /// inserts a well path that penetrates the edges following the sides of the supplied cells
  void AddWellEdgePath( const std::string& well_name, const std::vector<size_t>& cell_ids,
                        const std::vector<std::pair<size_t, size_t> >& edge_ids );
#endif

  template<class Container>  void GetRegions( Container& data );

  template<class Container>  void GetFaults( Container& data );

  template<class Container>  void GetWells( Container& wells );

  // SKM's replacement functions created in the debugging process
public:
  /// ReadPillarCoordinates -> Read_COORD
  bool Read_COORD( std::ifstream& ifs, char* text_line, size_t line_length );

  /// ReadCornerDepths -> Read_ZCORN
  bool Read_ZCORN( std::ifstream& ifs, char* text_line, size_t line_length );

  std::multimap<ijk, size_t>& IJKMap();


private:

  // Read different blocks of data from file
  bool ReadBlock( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadIncludeFile( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadUnknownBlock( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadDimensions( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadCellSizes( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadCellDepths( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadGridSpecs( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadPillarCoordinates( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadCornerDepths( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadActiveCells( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadWellSpecs( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadWellCompletionsData( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadExplicitFaceWellCompletionsData( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadExplicitNodeWellCompletionsData( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadFaultsData( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadFaultTransmissibilityMultipliers( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadBoxData( std::ifstream& ifs, char* text_line, size_t line_length );
  bool ReadScalarProperty( std::ifstream& ifs, char* text_line, size_t line_length, std::vector<csmp::ScalarVariable>& prop_data );
  bool ReadTensorProperty( std::ifstream& ifs, char* text_line, size_t line_length, std::vector<csmp::TensorVariable<3U> >& prop_data,
                           const std::string& compx_name, const std::string& compy_name, const std::string& compz_name );

  // Properties operations

  void WritePropertiesToVSet();

  void WriteScalarPropertyToVSet( csmp::VSet<3U>& vset,
                                  const CornerPointGrid_UoM& grid,
                                  const std::vector<csmp::ScalarVariable>& scalar_data,
                                  const std::string& property_name,
                                  const csmp::PLACEMENT& place );

  void WriteScalarPropertyToVSet( csmp::VSet<3U>& vset,
                                  const CornerPointGrid_UoM& grid,
                                  const std::vector<csmp::VectorVariable<3U> >& vector_data,
                                  const std::string& property_name,
                                  const csmp::PLACEMENT& place );

  void WriteScalarPropertyToVSet( csmp::VSet<3U>& vset,
                                  const CornerPointGrid_UoM& grid,
                                  const std::vector<csmp::TensorVariable<3U> >& tensor_data,
                                  const std::string& property_name,
                                  const csmp::PLACEMENT& place );

  void WriteVectorPropertyToVSet( csmp::VSet<3U>& vset,
                                  const CornerPointGrid_UoM& grid,
                                  const std::vector<csmp::VectorVariable<3U> >& vector_data,
                                  const std::string& property_name,
                                  const csmp::PLACEMENT& place );

  void WriteTensorPropertyToVSet( csmp::VSet<3U>& vset,
                                  const CornerPointGrid_UoM& grid,
                                  const std::vector<csmp::TensorVariable<3U> >& tensor_data,
                                  const std::string& property_name,
                                  const csmp::PLACEMENT& place );

  void SaveVectorProperty( size_t component,
                           const std::vector<csmp::ScalarVariable>& scalar_data,
                           std::vector<csmp::VectorVariable<3U> >&  vector_data );

  void SaveTensorProperty( size_t component,
                           const std::vector<csmp::ScalarVariable>& scalar_data,
                           std::vector<csmp::TensorVariable<3U> >&  vector_data );
  void ClearBefore();
  void ClearAfter();

  void AssignGridDimensions();

private:

  // reading process related
  // text file interface
  std::string                 keyword_;              // current keyword
  std::vector<std::string>    keyword_parameters_;   // current keyword parameters

                                                     // csmp model
  std::string                 model_name_;
  csmp::ModelTopology*        model_topology_;
  csmp::VSet<3U>*             vset_;

  // specific regions
  std::set<std::string>       regions_;
  std::set<std::string>       faults_;
  std::set<std::string>       wells_;

  // temporal grid data
  size_t                      NX_;
  size_t                      NY_;
  size_t                      NZ_;
  std::vector<double64>       zcorn_;
  CornerPointGrid_UoM         grid_;
  CellCenteredGrid            block_grid_;

#if 0
  // local mesh block
  std::vector<size_t>         box_;
#endif

  // properties data
  std::map<int, csmp::Parameter>           properties_;
  std::vector<csmp::ScalarVariable>       ntg_;
  std::vector<csmp::ScalarVariable>       poro_;
  std::vector<csmp::TensorVariable<3U> >  permxyz_;
  std::map<std::string, double>            multflt_;
  std::vector<csmp::ScalarVariable>       pressure_;
  std::vector<csmp::ScalarVariable>       swat_;
  std::vector<csmp::ScalarVariable>       soil_;
  std::vector<csmp::ScalarVariable>       sgas_;

  // regions data
  std::vector<csmp::ScalarVariable>       satnum_;
  std::vector<csmp::ScalarVariable>       pvtnum_;
  std::vector<csmp::ScalarVariable>       rocknum_;
  std::vector<csmp::ScalarVariable>       eqlnum_;
  std::vector<csmp::ScalarVariable>       fipnum_;

  std::map<std::string, EclipseFault>      faults_data_;

  // well specific data
  std::map<std::string, EclipseWell>       well_data_;
  std::map<std::string, EclipseWellPath>   well_face_path_;
};




/**
@brief PROPS and several other specifications from the RUNSPECS file record

@author R. Manasipov
@date 2014-2015

ASCII file interface to read corner point grid from ( *.GRDECL, *.grdecl ) file.
In order to build a native CSMP Model the provided grid information is assigned to
VSet and ModelTopology objects.
The regions that are constructed by default are:
% regions which include all the cells
MATRIX
MATRIX_ACTIVE
MATRIX_INACTIVE
% and specific fault regions

In addition to grid information the properties such as:
'permeability', 'porosity', 'pressure', 'saturation' can be readed as well.
User is able to control property type and placement.
For example one can create tensor or scalar permeability variables or
to place 'pressure' at the nodes or elements.

@code
KEYWORDS that can be currently handled:

--- Grid
DIMENS
DX
DY
DZ
SPECGRID
TOPS
COORD
ZCORN
FAULTS

--- Regions
ACTNUM
SATNUM
PVTNUM
ROCKNUM
EQLNUM

--- Props
NTG
PORO
PERMX            /// 3. assign well data
well.well_data_.push_back( wellcomp );
PERMY
PERMZ
PRESSURE
SWAT
SOIL
SGAS

--- Runspec
WELSPECS
COMPDAT (relies on WELLSPECS and can therefore not be defined in isolation)
@endcode

*/
class EclipseModelSettings {
public:

  explicit EclipseModelSettings( const std::string& mesh_file_prefix );
  ~EclipseModelSettings();

  void MeshSetup( bool exclude_inactive_cells,
                  bool tetra_mesh,
                  bool create_boundaries = true );

  void MeshSetup( const std::string& regions_file_prefix,
                  bool exclude_inactive_cells,
                  bool tetra_mesh,
                  bool create_boundaries = true );

  void MeshSetup( const std::set<std::string>& regions,
                  bool exclude_inactive_cells,
                  bool tetra_mesh,
                  bool create_boundaries = true );

  /// properties setup
  void SatNumPropertySetup( const std::string& prop_name = "satnum",
                            const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                            const csmp::PLACEMENT& prop_place = csmp::ELEMENT );

  void PvtNumPropertySetup( const std::string& prop_name = "ptvnum",
                            const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                            const csmp::PLACEMENT& prop_place = csmp::ELEMENT );

  /// usually this is called "rocknum"
  void RockNumPropertySetup( const std::string& prop_name = "rocktype",
                             const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                             const csmp::PLACEMENT& prop_place = csmp::ELEMENT );

  void EqlNumPropertySetup( const std::string& prop_name = "eqlnum",
                            const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                            const csmp::PLACEMENT& prop_place = csmp::ELEMENT );

  void FipNumPropertySetup( const std::string& prop_name = "fipnum",
                            const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                            const csmp::PLACEMENT& prop_place = csmp::ELEMENT );

  void NtgPropertySetup( const std::string& prop_name = "ntg",
                         const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                         const csmp::PLACEMENT& prop_place = csmp::ELEMENT );

  void PoroPropertySetup( const std::string& prop_name = "porosity",
                          const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                          const csmp::PLACEMENT& prop_place = csmp::ELEMENT );

  void PermPropertySetup( const std::string& prop_name = "permeability",
                          const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                          const csmp::PLACEMENT& prop_place = csmp::ELEMENT,
                          const std::string& prop_unit = "m2" );

  void PressurePropertySetup( const std::string& prop_name = "pressure",
                              const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                              const csmp::PLACEMENT& prop_place = csmp::NODE,
                              const std::string& prop_unit = "Pa" );

  void SwatPropertySetup( const std::string& prop_name = "saturation water",
                          const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                          const csmp::PLACEMENT& prop_place = csmp::NODE );

  void SoilPropertySetup( const std::string& prop_name = "saturation oil",
                          const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                          const csmp::PLACEMENT& prop_place = csmp::NODE );

  void SgasPropertySetup( const std::string& prop_name = "saturation gas",
                          const csmp::VARIABLE_TYPE& prop_type = csmp::SCALAR,
                          const csmp::PLACEMENT& prop_place = csmp::NODE );

  /// mesh setting
  std::string           mesh_file_prefix_;       ///< name of model or main mesh file
  std::string           regions_file_prefix_;    ///< regions file prefix
  bool                  create_boundaries_;      ///< true = create boundaries around model, false = do not create boundaries
  std::set<std::string> regions_;                ///< name of regions to be used
  bool                  tetra_mesh_;             ///< tetrahedral mesh
  bool                  exclude_inactive_cells_; ///< do not consider inactive cells

                                                 /// property settings
  std::map<int, csmp::Parameter>   properties_;  // eclipse properties ( porosity, permeability, transmissibility, etc. )

private:
  EclipseModelSettings() {}; // should not be used
};



// common reading functinality
bool isEclipseCommentLine( char* text_line );
bool isEclipseEndOfBlock( char* text_line );
bool isEclipseLineWithEndOfBlock( char* text_line );

bool readEclipseKeyword( std::ifstream& ifs, char* text_line, size_t line_length,
                         std::string& keyword, std::vector<std::string>& keyword_parameters );

bool readEclipseNonBlankLine( std::ifstream& ifs, char* text_line, size_t line_length,
                              bool& endOfblock );

int  readEclipseFirstLineInBlock( std::ifstream& ifs, char* text_line, size_t line_length, bool& endOfblock );

int  readEclipseIncludeFileName( std::ifstream& ifs, char* text_line, size_t line_length,
                                 std::string& fname );

template<class VarType>
bool readEclipseValue( const std::string& text_line, bool& default_value, size_t& num, VarType& value );

// reading specific block's of data
int readEclipseDimensions( size_t& NX, size_t& NY, size_t& NZ,
                           std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

int readEclipseGridSpecs( size_t& NX, size_t& NY, size_t& NZ,
                          std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

int readEclipsePillarCoordinates( size_t& NX, size_t& NY,
                                  CornerPointGrid_UoM& grid,
                                  std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                                  );

int readEclipseCornerDepths( size_t NX, size_t NY, size_t& NZ,
                             CornerPointGrid_UoM& grid,
                             std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                             );

int readEclipseActiveCells( std::vector<uint8_t>& cell_activity,
                            std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                            );

int readEclipseCellData( std::vector<csmp::ScalarVariable>& values,
                         std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                         );

int readEclipseWellSpecs( std::map<std::string, EclipseWell> &well_data,
                          std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                          );

int readEclipseWellCompletionsData( size_t NX, size_t NY, size_t NZ,
                                    std::map<std::string, EclipseWell>& well_data,
                                    std::map<std::string, EclipseWellPath>& well_path,
                                    std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                                    );

int readEclipseExplicitFaceWellCompletionsData( size_t NX, size_t NY, size_t NZ,
                                                std::map<std::string, EclipseWell>& well_data,
                                                std::map<std::string, EclipseWellPath>& well_path,
                                                std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                                                );

int readEclipseExplicitNodeWellCompletionsData( size_t NX, size_t NY, size_t NZ,
                                                std::map<std::string, EclipseWell>& well_data,
                                                std::map<std::string, EclipseWellPath>& well_path,
                                                std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                                                );

int readEclipseFaultData( size_t NX, size_t NY, size_t NZ,
                          std::map<std::string, EclipseFault>& faults_data,
                          std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                          );

int readEclipseFaultTransmissibilityMultipliers( size_t NX, size_t NY, size_t NZ,
                                                 std::map<std::string, double>& multflt,
                                                 std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                                                 );

int readEclipseBoxData( std::vector<size_t>& box_data,
                        std::ifstream& ifs, char* text_line, size_t line_length, bool verbose
                        );

int skipEclipseBlock( std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

void removeSymbolsFromString( std::string &str, const char* symbolsToRemove );


void addWellPath( size_t NX, size_t NY, size_t NZ,
                  const std::string& well_name,
                  const std::vector<ijk>& cell_ids,
                  std::map<std::string, EclipseWellPath>& well_path );

void addWellPath( const std::string& well_name,
                  const std::vector<ijk>& cell_ids,
                  const std::vector<std::pair<ijk, CORNER_POINT_CELL_FACE_INDEX> >& face_ids,
                  std::map<std::string, EclipseWellPath>& well_path );

void addWellPath( const std::string& well_name,
                  const std::vector<ijk>& cell_ids,
                  std::pair<ijk, CORNER_POINT_CELL_FACE_INDEX> face_ids,
                  std::map<std::string, EclipseWellPath>& well_path );

void addWellPath( const std::string& well_name,
                  const std::vector<ijk>& cell_ids,
                  std::pair<CORNER_POINT_CELL_FACE_INDEX, CORNER_POINT_CELL_FACE_INDEX> face_ids,
                  std::map<std::string, EclipseWellPath>& well_path );

void addWellPath( const std::string& well_name,
                  const ijk& cell_ids,
                  std::pair<ijk, CORNER_POINT_CELL_FACE_INDEX> face_id,
                  std::map<std::string, EclipseWellPath>& well_path );



} // eclipse

}// end namespace csmp

#endif

