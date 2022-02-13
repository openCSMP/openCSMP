#include "CornerPointGrid.h"
#include "PropertyData.h"
#include "CSMP_highLevelUtilities.h"
#include "STL_utilities.h"
#include "Pillar.h"
#include "CellGenerator.h"
#include "Box.h"

#include "ErrorHandler.h"
#include "EclipseInterface.h"
#include "PolygonGrid.h"
#include "VSet.h"
#include "ModelTopology.h"


#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearPrism.h"
#include "IsoparametricLinearLineElement.h"

#define CSMP_SIZE_OF_ARRAY(a)  (sizeof(a) / sizeof(a[0]))

using namespace std;

namespace csmp {

// CORNER POINT GRID

CornerPointGrid::CornerPointGrid()
  : NX_( 0 ),
  NY_( 0 ),
  NZ_( 0 ),
  NX_x_NY_( 0 )
{
}

CornerPointGrid::~CornerPointGrid()
{
}


void CornerPointGrid::AssignDimensions( size_t nx, size_t ny, size_t nz )
{
  NX_ = nx;
  NY_ = ny;
  NZ_ = nz;
}


void CornerPointGrid
::InitializeGridSpecs()
{
  std::cout << "InitializeGridSpecs\n";
  NX_x_NY_ = NX_ * NY_;
  elements_ = NX_x_NY_ * NZ_;
  if ( cell_activity_.empty() ) {
    cell_activity_.resize( elements_, 1 );
    active_elements_ = elements_;
  }
  else {
    active_elements_ = 0;
    for ( auto activity : cell_activity_ ) {
      if ( activity ) ++active_elements_;
    }
  }
}

/**
x stays, y = -z and z=y.
*/
void CornerPointGrid
::ConvertFromReservoirToCSMPcoordinateSystem( csmp::Point<3U>& pt ) const
{
  const double y( pt[1U] );
  pt[0U] = pt[0U];
  pt[1U] = -pt[2U];
  pt[2U] = y;
};


uint8_t&
CornerPointGrid::CellActivity( size_t i, size_t j, size_t k )
{
  assert( i < NX_ && j < NY_ && k < NZ_ );
  return cell_activity_[i + j * NX_ + k * NX_x_NY_];
}


std::vector<uint8_t>&
CornerPointGrid::GetCellActivity()
{
  return cell_activity_;
}



/**
Builds the pillars and columns.
*/
void CornerPointGrid::ConstructPillarsAndColumns( const std::vector<double>& zcorn )
{
  const size_t NXxNY = NX_ * NY_;
  size_t node_count = 0;

  {
    // 1. Build pillars
    std::vector<double> zcoord;
    zcoord.reserve( 4 * (NZ_ + 1) );
    for ( auto i = 0; i <= NX_; ++i )
    {
      for ( size_t j = 0; j <= NY_; ++j )
      {
        for ( size_t k = 0; k < NZ_; ++k ) {
          if ( i < NX_ && j < NY_ && CellActivity( i, j, k ) ) {
            zcoord.push_back( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 0] ); // t_nw
            zcoord.push_back( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 4] ); // b_nw
          }

          if ( i > 0 && j < NY_ && CellActivity( i - 1, j, k ) ) {
            zcoord.push_back( zcorn[((i - 1) + j * NX_ + k * NXxNY) * 8 + 1] ); // t_ne
            zcoord.push_back( zcorn[((i - 1) + j * NX_ + k * NXxNY) * 8 + 5] ); // b_ne
          }

          if ( i < NX_ && j > 0 && CellActivity( i, j - 1, k ) ) {
            zcoord.push_back( zcorn[(i + (j - 1) * NX_ + k * NXxNY) * 8 + 2] ); // t_sw
            zcoord.push_back( zcorn[(i + (j - 1) * NX_ + k * NXxNY) * 8 + 6] ); // b_sw
          }

          if ( i > 0 && j > 0 && CellActivity( i - 1, j - 1, k ) ) {
            zcoord.push_back( zcorn[((i - 1) + (j - 1) * NX_ + k * NXxNY) * 8 + 3] ); // t_se
            zcoord.push_back( zcorn[((i - 1) + (j - 1) * NX_ + k * NXxNY) * 8 + 7] ); // b_se
          }
        }

        sortAndUnique( zcoord );
        auto& p = (*this)(i, j);
        p.SetZCoords( zcoord );
        p.SetFirstNodeNum( node_count );
        node_count += zcoord.size();
        zcoord.clear();
      }
    }
  }
  ordinaryNodes_ = node_count;
  std::cerr << " " << node_count << " unique nodes detected\n";

  // 2. Build columns
  size_t cell_count = 0;
  size_t fully_degenerate_cells = 0;
  {
    std::vector<ColumnCell> cells;
    cells.reserve( NZ_ );
    for ( auto i = 0; i < NX_; ++i )
    {
      for ( size_t j = 0; j < NY_; ++j )
      {
        std::pair<size_t, size_t> index( i, j );
        Pillar& p0 = (*this)(i + 0, j + 0);
        Pillar& p1 = (*this)(i + 1, j + 0);
        Pillar& p2 = (*this)(i + 1, j + 1);
        Pillar& p3 = (*this)(i + 0, j + 1);

        for ( size_t k = 0; k < NZ_; ++k )
        {
          if ( !CellActivity( i, j, k ) ) {
            continue;
          }

          ColumnCell cell;
          cell.k = k;
          cell.z[0][0] = p0.FindPoint( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 0] ); // t_nw
          cell.z[0][1] = p0.FindPoint( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 4] ); // b_nw
          cell.z[1][0] = p1.FindPoint( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 1] ); // t_ne
          cell.z[1][1] = p1.FindPoint( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 5] ); // b_ne
          cell.z[2][0] = p2.FindPoint( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 3] ); // t_se
          cell.z[2][1] = p2.FindPoint( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 7] ); // b_se
          cell.z[3][0] = p3.FindPoint( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 2] ); // t_sw
          cell.z[3][1] = p3.FindPoint( zcorn[(i + j * NX_ + k * NXxNY) * 8 + 6] ); // b_sw

          // If all four corners are degenerate, the cell is fully degenerate.
          if ( cell.z[0][0] == cell.z[0][1]
               && cell.z[1][0] == cell.z[1][1]
               && cell.z[2][0] == cell.z[2][1]
               && cell.z[3][0] == cell.z[3][1] ) {
            ++fully_degenerate_cells;
            continue;
          }
          ++cell_count;
          cells.push_back( cell );
        }

        columns_.emplace( index, cells );
        cells.clear();
      }
    }
  }
  std::cerr << " " << cell_count << " unique active cells detected\n";
  if ( fully_degenerate_cells > 0 ) {
    std::cerr << " " << fully_degenerate_cells << " fully degenerate cells detected\n";
  }

  std::cerr << "Pillars and columns built\n";

} // end ConstructPillarsAndColumns




void CornerPointGrid::ConstructFiniteElementsFromColumns( VSet<3U>& vset )
{
  vset.HybridElementTypeMesh( true );

  // 1. Classify the cells

  size_t skewCells = 0;
  for ( auto& index_column : columns_ ) {
    auto& column = index_column.second;
    auto i = index_column.first.first;
    size_t j = index_column.first.second;

    const size_t iNrCells = column.cells_.size();
    for ( size_t iCell = 0; iCell < iNrCells; ++iCell ) {
      auto& cell = column.cells_[iCell];

      Pillar& p0 = (*this)(i + 0, j + 0); //nw
      Pillar& p1 = (*this)(i + 1, j + 0); //ne
      Pillar& p2 = (*this)(i + 1, j + 1); //se
      Pillar& p3 = (*this)(i + 0, j + 1); //sw

      double z[4][2];
      z[0][0] = p0.GetZCoord( cell.z[0][0] );
      z[0][1] = p0.GetZCoord( cell.z[0][1] );
      z[1][0] = p1.GetZCoord( cell.z[1][0] );
      z[1][1] = p1.GetZCoord( cell.z[1][1] );
      z[2][0] = p2.GetZCoord( cell.z[2][0] );
      z[2][1] = p2.GetZCoord( cell.z[2][1] );
      z[3][0] = p3.GetZCoord( cell.z[3][0] );
      z[3][1] = p3.GetZCoord( cell.z[3][1] );

      //AB: Classify the degeneracy
      uint8_t classification = 0;
      for ( size_t v = 0; v < 4; ++v ) {
        if ( cell.z[v][0] == cell.z[v][1] ) {
          classification |= (1 << (3 - v));
        }
      }
      cell.classification = static_cast<ECLIPSE_CELL_CLASSIFICATION>(classification);
    }
  }


  // 2. Constructing Elements
  generator_ = new CellGenerator( *this );
  generator_->ordinaryNodes.reserve( ordinaryNodes_ );

  std::vector<Point<3U>> tp_pts;
  std::vector<Point<3U>> bt_pts;
  for ( auto i = 0; i <= NX_; ++i ) {
    for ( size_t j = 0; j <= NY_; ++j ) {
      Pillar& pillar = (*this)(i, j);      
      for ( size_t k = 0; k < pillar.GetNumPoints(); ++k ) {
        generator_->ordinaryNodes.push_back( pillar.GetPoint( k ) );        
      }

      // keep the top & the bottom points
      Point<3U> bt = pillar.EndPoint();
      Point<3U> tp = pillar.StartPoint();
      ConvertFromReservoirToCSMPcoordinateSystem( bt );
      ConvertFromReservoirToCSMPcoordinateSystem( tp );
      bt_pts.push_back( bt );
      tp_pts.push_back( tp );
    }
  }
    
  badHexahedra_ = 0;
  badPyramids_ = 0;
  badTetrahedra_ = 0;
  badPrisms_ = 0;
  badLines_ = 0;

  //Degeneration process starts!
  cout << "\nCornerPointGrid::ConstructFiniteElementsFromColumns: creating elements...";  
  for ( auto column = columns_.begin(); column != columns_.end(); column++ ) {
    Column& Col = column->second;
    for ( size_t k = 0; k < Col.cells_.size(); k++ ) {
      ColumnCell&  cell = Col.cells_[k];
      pair<size_t, size_t> index = column->first;
      auto i = index.first;
      size_t j = index.second;

      generator_->setIJ( i, j );

      ECLIPSE_CELL_CLASSIFICATION cellType = cell.classification;
      
      /// referred to the classifications of the cell above and the cell beneath  
      ColumnCell* cellAbove = nullptr;
      ColumnCell* cellBeneath = nullptr;

      if ( k > 0 ) {
        cellAbove = &Col.cells_[k - 1];
      }

      if ( k < Col.cells_.size() - 1 ) {
        cellBeneath = &Col.cells_[k + 1];
      }

      
      assert( generator_->elementID == generator_->plist.size() );
      assert( generator_->elementID == generator_->fem_types.size() );

      switch ( cellType ) {
        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_HEXAHEDRON:				// 0000
          if ( !ConstructEclipseCell0000( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_HEXAHEDRON(" << i << "," << j << "," << k << ")";            
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_310_312:		// 0001
          if ( !ConstructEclipseCell0001( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PYRAMIDS_310_312(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_201_203:		// 0010
          if ( !ConstructEclipseCell0010( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PYRAMIDS_201_203(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PRISM_23:				// 0011
          if ( !ConstructEclipseCell0011( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PRISM_23(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_130_132:		// 0100
          if ( !ConstructEclipseCell0100( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PYRAMIDS_130_132(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_TETRAHEDRONS_130_132:	// 0101
          if ( !ConstructEclipseCell0101( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_TETRAHEDRONS_130_132(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PRISM_12:				// 0110
          if ( !ConstructEclipseCell0110( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PRISM_12(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_0:				// 0111
          if ( !ConstructEclipseCell0111( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PYRAMID_0(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_021_023:		// 1000
          if ( !ConstructEclipseCell1000( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PYRAMIDS_021_023(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PRISM_03:				// 1001
          if ( !ConstructEclipseCell1001( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PRISM_03(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_TETRAHEDRONS_021_023:	// 1010
          if ( !ConstructEclipseCell1010( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_TETRAHEDRONS_021_023(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_1:				// 1011
          if ( !ConstructEclipseCell1011( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PYRAMID_1(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PRISM_01:				// 1100
          if ( !ConstructEclipseCell1100( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PRISM_01(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_2:				// 1101
          if ( !ConstructEclipseCell1101( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PYRAMID_2(" << i << "," << j << "," << k << ")";
          break;

        case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_3:				// 1110
          if ( !ConstructEclipseCell1110( cell, i, j, k, cellAbove, cellBeneath ) )
            cout << "\n detected a broken cell: ECLIPSE_CELL_PYRAMID_3(" << i << "," << j << "," << k << ")";
          break;
      }
    }
  }

  std::cerr << "\n\n";
  std::cerr << " removed invalid hexahedra: " << badHexahedra_ << '\n';
  std::cerr << " removed invalid pyramids: " << badPyramids_ << '\n';
  std::cerr << " removed invalid tetrahedra: " << badTetrahedra_ << '\n';
  std::cerr << " removed invalid prisms: " << badPrisms_ << '\n';
  std::cerr << " removed invalid lines: " << badLines_ << '\n';

  // 3. store node coordinates
  deque<double> x, y, z;
  for ( auto i = 0; i <= NX_; ++i ) {
    for ( size_t j = 0; j <= NY_; ++j ) {
      Pillar& pillar = (*this)(i, j);
      for ( size_t k = 0; k < pillar.GetNumPoints(); ++k ) {
        csmp::Point<3u> point = pillar.GetPoint( k );
        // convert to CSMP coordinate
        ConvertFromReservoirToCSMPcoordinateSystem( point );
        x.push_back( point[0] );
        y.push_back( point[1] );
        z.push_back( point[2] );
      }
    }
  }

  for ( auto p : generator_->extraNodes ) {
    ConvertFromReservoirToCSMPcoordinateSystem( p );
    x.push_back( p[0] );
    y.push_back( p[1] );
    z.push_back( p[2] );
  }
  vset.AddXYZ( x, y, z );

  // 4. assigning top & bottom boundary flags on nodes  
  cout << "\nCornerPointGrid::ConstructFiniteElementsFromColumns: assigning top & bottom boundary flags on nodes...\n";
  struct isEqual {
    isEqual( const Point<3U>& pt ) : m_pt( pt ) {};
    bool operator()( const Point<3U>& lpt )
    {
      if ( Point<3U>( lpt - m_pt ).Length() > 0.02 ) return false; // the threshold 0.02 was optimized for the Otway model since there are broken elements in the top or bottom surfaces
      else return true;
    };

    Point<3U> m_pt;
  };
  
  vector<std::int8_t> pbflags( vset.Vertices(), 0 ); // boundary flags
  for ( auto i = 0U; i < vset.Vertices(); ++i )
    {
      vector<double> coord( 3U );
      for ( size_t j = 0U; j<3U; ++j ) coord[j] = vset.P( j, i );
      Point<3U> pt( coord );
      if ( std::find_if( tp_pts.begin(), tp_pts.end(), isEqual( pt ) ) != tp_pts.end() )
        pbflags[i] = BOX_BOUNDARY::TOP;
      else if ( std::find_if( bt_pts.begin(), bt_pts.end(), isEqual( pt ) ) != bt_pts.end() )
        pbflags[i] = BOX_BOUNDARY::BOTTOM;
      else
        pbflags[i] = BOX_BOUNDARY::IRREGULAR;
    }

  // 5. Set up the rest of the vset
  // There is no neighbor information in the Eclipse data(*.grdecl). The information will be created later.
  vset.RemovePfverts();
  vset.ResizePlist( generator_->plist.size() );
  vset.AddPlist( generator_->plist.begin(), generator_->plist.end() );
  vset.AddElementTypes( generator_->fem_types.begin(), generator_->fem_types.end() );
  vset.AddBFlags( pbflags.begin(), pbflags.end() );
}




bool CornerPointGrid::ConstructEclipseCell0000( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  bool faceAboveIsQuad = !cellAbove || generator_->getShapeOfBottomFace( *cellAbove ) == FACE_TYPE::FULL_QUAD;
  bool faceBeneathIsQuad = !cellBeneath || generator_->getShapeOfTopFace( *cellBeneath ) == FACE_TYPE::FULL_QUAD;

  size_t invalid_hexs = 0;
  size_t invalid_pyrs = 0;
  size_t invalid_tets = 0;
  size_t invalid_pris = 0;

  std::vector<uint32_t> new_elements;

  if ( faceAboveIsQuad ) {
    if ( faceBeneathIsQuad ) {
      if ( generator_->ConstructHexahedron( cell ) ) {
        size_t eid = generator_->EmitHexahedron( cell );
        addElementToMap( i, j, k, eid );
        new_elements.push_back( eid );
      }
      else {
        invalid_hexs++;
      }

      badHexahedra_ += invalid_hexs;
    }
    else {
      size_t centroid = generator_->generateCellCentroid( cell );
      if ( centroid == 0 ) return false;

      switch ( generator_->getShapeOfTopFace( *cellBeneath ) ) {
        case FACE_TYPE::SPLIT_02_OR_46:
          // Top face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 1, 2, 3, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Left face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 4, 5, 1, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Right face
          if ( generator_->ConstructPyramidOnFace( cell, 3, 2, 6, 7, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Front face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 3, 7, 4, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Back face
          if ( generator_->ConstructPyramidOnFace( cell, 1, 5, 6, 2, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Bottom faces
          if ( generator_->ConstructTetrahedronOnFace( cell, 6, 5, 4, centroid ) ) {
            size_t eid = generator_->EmitTetrahedron( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_tets++;
          }
          if ( generator_->ConstructTetrahedronOnFace( cell, 6, 4, 7, centroid ) ) {
            size_t eid = generator_->EmitTetrahedron( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_tets++;
          }

          badPyramids_ += invalid_pyrs;
          badTetrahedra_ += invalid_tets;

          break;

        case FACE_TYPE::SPLIT_13_OR_57:
          // Top face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 1, 2, 3, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Left face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 4, 5, 1, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Right face
          if ( generator_->ConstructPyramidOnFace( cell, 3, 2, 6, 7, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Front face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 3, 7, 4, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Back face
          if ( generator_->ConstructPyramidOnFace( cell, 1, 5, 6, 2, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Bottom faces
          if ( generator_->ConstructTetrahedronOnFace( cell, 7, 6, 5, centroid ) ) {
            size_t eid = generator_->EmitTetrahedron( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_tets++;
          }
          if ( generator_->ConstructTetrahedronOnFace( cell, 7, 5, 4, centroid ) ) {
            size_t eid = generator_->EmitTetrahedron( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_tets++;
          }

          badPyramids_ += invalid_pyrs;
          badTetrahedra_ += invalid_tets;

          break;

        default:
          throw csmp::Exception( FATAL_ERROR, "CornerPointGrid::ConstructFiniteElementsFromColumns", "Hexahedron with unknown bottom face" );
      }
    }
  }
  else {
    if ( faceBeneathIsQuad ) {
      size_t centroid = generator_->generateCellCentroid( cell );
      if ( centroid == 0 ) return false;

      switch ( generator_->getShapeOfBottomFace( *cellAbove ) ) {
        case FACE_TYPE::SPLIT_02_OR_46:
          // Left face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 4, 5, 1, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Right face
          if ( generator_->ConstructPyramidOnFace( cell, 3, 2, 6, 7, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Front face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 3, 7, 4, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Back face
          if ( generator_->ConstructPyramidOnFace( cell, 1, 5, 6, 2, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Bottom face
          if ( generator_->ConstructPyramidOnFace( cell, 7, 6, 5, 4, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Top faces
          if ( generator_->ConstructTetrahedronOnFace( cell, 0, 1, 2, centroid ) ) {
            size_t eid = generator_->EmitTetrahedron( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_tets++;
          }
          if ( generator_->ConstructTetrahedronOnFace( cell, 0, 2, 3, centroid ) ) {
            size_t eid = generator_->EmitTetrahedron( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_tets++;
          }

          badPyramids_ += invalid_pyrs;
          badTetrahedra_ += invalid_tets;

          break;

        case FACE_TYPE::SPLIT_13_OR_57:
          // Left face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 4, 5, 1, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Right face
          if ( generator_->ConstructPyramidOnFace( cell, 3, 2, 6, 7, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Front face
          if ( generator_->ConstructPyramidOnFace( cell, 0, 3, 7, 4, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Back face
          if ( generator_->ConstructPyramidOnFace( cell, 1, 5, 6, 2, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Bottom face
          if ( generator_->ConstructPyramidOnFace( cell, 7, 6, 5, 4, centroid ) ) {
            size_t eid = generator_->EmitPyramid( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_pyrs++;
          }
          // Top faces
          if ( generator_->ConstructTetrahedronOnFace( cell, 1, 2, 3, centroid ) ) {
            size_t eid = generator_->EmitTetrahedron( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_tets++;
          }
          if ( generator_->ConstructTetrahedronOnFace( cell, 1, 3, 0, centroid ) ) {
            size_t eid = generator_->EmitTetrahedron( cell );
            addElementToMap( i, j, k, eid );
            new_elements.push_back( eid );
          }
          else {
            invalid_tets++;
          }

          badPyramids_ += invalid_pyrs;
          badTetrahedra_ += invalid_tets;

          break;

        default:
          throw csmp::Exception( FATAL_ERROR, "CornerPointGrid::ConstructFiniteElementsFromColumns", "Hexahedron with unknown bottom face" );
      }
    }
    else {
      auto topShape = generator_->getShapeOfBottomFace( *cellAbove );
      auto bottomShape = generator_->getShapeOfTopFace( *cellBeneath );

      size_t centroid = generator_->generateCellCentroid( cell );
      if ( centroid == 0 ) return false;

      switch ( topShape ) {
        case FACE_TYPE::SPLIT_02_OR_46:
          switch ( bottomShape ) {
            case FACE_TYPE::SPLIT_02_OR_46:
              if ( generator_->ConstructPrism( cell, 6, 5, 4, 2, 1, 0 ) ) {
                size_t eid = generator_->EmitPrism( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pris++;
              }
              if ( generator_->ConstructPrism( cell, 0, 2, 3, 4, 6, 7 ) ) {
                size_t eid = generator_->EmitPrism( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pris++;
              }
              badPrisms_ += invalid_pris;

              break;

            case FACE_TYPE::SPLIT_13_OR_57:
              if ( generator_->ConstructPyramidOnFace( cell, 0, 4, 5, 1, centroid ) ) {
                size_t eid = generator_->EmitPyramid( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pyrs++;
              }
              // Right face
              if ( generator_->ConstructPyramidOnFace( cell, 3, 2, 6, 7, centroid ) ) {
                size_t eid = generator_->EmitPyramid( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pyrs++;
              }
              // Front face
              if ( generator_->ConstructPyramidOnFace( cell, 0, 3, 7, 4, centroid ) ) {
                size_t eid = generator_->EmitPyramid( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pyrs++;
              }
              // Back face
              if ( generator_->ConstructPyramidOnFace( cell, 1, 5, 6, 2, centroid ) ) {
                size_t eid = generator_->EmitPyramid( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pyrs++;
              }
              // Top faces
              if ( generator_->ConstructTetrahedronOnFace( cell, 0, 1, 2, centroid ) ) {
                size_t eid = generator_->EmitTetrahedron( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_tets++;
              }
              if ( generator_->ConstructTetrahedronOnFace( cell, 0, 2, 3, centroid ) ) {
                size_t eid = generator_->EmitTetrahedron( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_tets++;
              }
              // Bottom faces
              if ( generator_->ConstructTetrahedronOnFace( cell, 7, 6, 5, centroid ) ) {
                size_t eid = generator_->EmitTetrahedron( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_tets++;
              }
              if ( generator_->ConstructTetrahedronOnFace( cell, 7, 5, 4, centroid ) ) {
                size_t eid = generator_->EmitTetrahedron( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_tets++;
              }

              badPyramids_ += invalid_pyrs;
              badTetrahedra_ += invalid_tets;

              break;

            default:
              throw csmp::Exception( FATAL_ERROR, "CornerPointGrid::ConstructFiniteElementsFromColumns", "Hexahedron with unknown bottom face" );
          }
          break;
        case FACE_TYPE::SPLIT_13_OR_57:
          switch ( bottomShape ) {
            case FACE_TYPE::SPLIT_02_OR_46:
              // Left face
              if ( generator_->ConstructPyramidOnFace( cell, 0, 4, 5, 1, centroid ) ) {
                size_t eid = generator_->EmitPyramid( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pyrs++;
              }
              // Right face
              if ( generator_->ConstructPyramidOnFace( cell, 3, 2, 6, 7, centroid ) ) {
                size_t eid = generator_->EmitPyramid( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pyrs++;
              }
              // Front face
              if ( generator_->ConstructPyramidOnFace( cell, 0, 3, 7, 4, centroid ) ) {
                size_t eid = generator_->EmitPyramid( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pyrs++;
              }
              // Back face
              if ( generator_->ConstructPyramidOnFace( cell, 1, 5, 6, 2, centroid ) ) {
                size_t eid = generator_->EmitPyramid( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pyrs++;
              }
              // Top faces
              if ( generator_->ConstructTetrahedronOnFace( cell, 1, 2, 3, centroid ) ) {
                size_t eid = generator_->EmitTetrahedron( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_tets++;
              }
              if ( generator_->ConstructTetrahedronOnFace( cell, 1, 3, 0, centroid ) ) {
                size_t eid = generator_->EmitTetrahedron( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_tets++;
              }
              // Bottom faces
              if ( generator_->ConstructTetrahedronOnFace( cell, 6, 5, 4, centroid ) ) {
                size_t eid = generator_->EmitTetrahedron( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_tets++;
              }
              if ( generator_->ConstructTetrahedronOnFace( cell, 6, 4, 7, centroid ) ) {
                size_t eid = generator_->EmitTetrahedron( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_tets++;
              }

              badPyramids_ += invalid_pyrs;
              badTetrahedra_ += invalid_tets;

              break;

            case FACE_TYPE::SPLIT_13_OR_57:
              if ( generator_->ConstructPrism( cell, 5, 4, 7, 1, 0, 3 ) ) {
                size_t eid = generator_->EmitPrism( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pris++;
              }
              if ( generator_->ConstructPrism( cell, 7, 6, 5, 3, 2, 1 ) ) {
                size_t eid = generator_->EmitPrism( cell );
                addElementToMap( i, j, k, eid );
                new_elements.push_back( eid );
              }
              else {
                invalid_pris++;
              }
              badPrisms_ += invalid_pris;

              break;

            default:
              throw csmp::Exception( FATAL_ERROR, "CornerPointGrid::ConstructFiniteElementsFromColumns", "Hexahedron with unknown bottom face" );
          }
          break;
      }
    }
  }

  if ( invalid_hexs > 0 || invalid_pyrs > 0 || invalid_tets > 0 || invalid_pris > 0 )
    return false;

  return true;
}


bool CornerPointGrid::ConstructEclipseCell0001( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructPyramidOnFace( cell, 4, 5, 1, 0, generator_->getNodeID( cell, 3 ) ) ) {
    size_t eid = generator_->EmitPyramid( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else
    invalid_elements++;

  if ( generator_->ConstructPyramidOnFace( cell, 1, 5, 6, 2, generator_->getNodeID( cell, 3 ) ) ) {
    size_t eid = generator_->EmitPyramid( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else
    invalid_elements++;

  badPyramids_ += invalid_elements;

  if ( invalid_elements > 0 )
    return false;

  return true;
}

bool CornerPointGrid::ConstructEclipseCell0010( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructPyramidOnFace( cell, 4, 5, 1, 0, generator_->getNodeID( cell, 2 ) ) ) {
    size_t eid = generator_->EmitPyramid( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  if ( generator_->ConstructPyramidOnFace( cell, 3, 7, 4, 0, generator_->getNodeID( cell, 2 ) ) ) {
    size_t eid = generator_->EmitPyramid( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badPyramids_ += invalid_elements;

  if ( invalid_elements > 0 )
    return false;

  return true;
}

bool CornerPointGrid::ConstructEclipseCell0011( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  generator_->plist.emplace( generator_->elementID, generator_->degenerateToOnePrismAtEdge23( cell ) );
  addElementToMap( i, j, k, generator_->elementID );
  generator_->fem_types.push_back( ISOPARAMETRIC_LINEAR_PRISM );
  ++generator_->elementID;

  return true;
}

bool CornerPointGrid::ConstructEclipseCell0100( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructPyramidOnFace( cell, 3, 7, 4, 0, generator_->getNodeID( cell, 1 ) ) ) {
    size_t eid = generator_->EmitPyramid( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else
    invalid_elements++;

  if ( generator_->ConstructPyramidOnFace( cell, 2, 6, 7, 3, generator_->getNodeID( cell, 1 ) ) ) {
    size_t eid = generator_->EmitPyramid( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else
    invalid_elements++;

  badPyramids_ += invalid_elements;

  if ( invalid_elements > 0 )
    return false;

  return true;
}

bool CornerPointGrid::ConstructEclipseCell0101( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  auto elementList = generator_->degenerateToTwoTetrahedraAtEdge13( cell );
  generator_->plist.emplace( generator_->elementID, elementList[0] );
  addElementToMap( i, j, k, generator_->elementID );
  generator_->fem_types.push_back( ISOPARAMETRIC_LINEAR_TETRAHEDRON );
  ++generator_->elementID;
  generator_->plist.emplace( generator_->elementID, elementList[1] );
  addElementToMap( i, j, k, generator_->elementID );
  generator_->fem_types.push_back( ISOPARAMETRIC_LINEAR_TETRAHEDRON );
  ++generator_->elementID;

  return true;
}

bool CornerPointGrid::ConstructEclipseCell0110( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructPrism( cell, 0, 4, 1, 3, 7, 2 ) ) {
    size_t eid = generator_->EmitPrism( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badPrisms_ += invalid_elements;

  if ( invalid_elements > 0 )
    return false;

  return true;
}

//Eclipse cell type: ECLIPSE_CELL_PYRAMID_0
bool CornerPointGrid::ConstructEclipseCell0111( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructTetrahedronOnFace( cell, 7, 5, 4, generator_->getNodeID( cell, 0 ) ) ) {
    size_t eid = generator_->EmitTetrahedron( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badTetrahedra_ += invalid_elements;

  if ( invalid_elements > 0 )
    return false;

  return true;
}

bool CornerPointGrid::ConstructEclipseCell1000( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructPyramidOnFace( cell, 5, 6, 2, 1, generator_->getNodeID( cell, 0 ) ) ) {
    size_t eid = generator_->EmitPyramid( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  if ( generator_->ConstructPyramidOnFace( cell, 7, 3, 2, 6, generator_->getNodeID( cell, 0 ) ) ) {
    size_t eid = generator_->EmitPyramid( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badPyramids_ += invalid_elements;

  if ( invalid_elements > 0 ) 
    return false;

  return true;
}


bool CornerPointGrid::ConstructEclipseCell1001( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  generator_->plist.emplace( generator_->elementID, generator_->degenerateToOnePrismAtEdge30( cell ) );
  addElementToMap( i, j, k, generator_->elementID );
  generator_->fem_types.push_back( ISOPARAMETRIC_LINEAR_PRISM );
  ++generator_->elementID;

  return true;
}

//Eclipse cell type: ECLIPSE_CELL_TETRAHEDRONS_021_023
bool CornerPointGrid::ConstructEclipseCell1010( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  auto elementList = generator_->degenerateToTwoTetrahedrasAtEdge02( cell );
  generator_->plist.emplace( generator_->elementID, elementList[0] );
  addElementToMap( i, j, k, generator_->elementID );
  generator_->fem_types.push_back( ISOPARAMETRIC_LINEAR_TETRAHEDRON );
  ++generator_->elementID;
  generator_->plist.emplace( generator_->elementID, elementList[1] );
  addElementToMap( i, j, k, generator_->elementID );
  generator_->fem_types.push_back( ISOPARAMETRIC_LINEAR_TETRAHEDRON );
  ++generator_->elementID;

  return true;
}

//Eclipse cell type: ECLIPSE_CELL_PYRAMID_1
bool CornerPointGrid::ConstructEclipseCell1011( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructTetrahedronOnFace( cell, 6, 5, 4, generator_->getNodeID( cell, 1 ) ) ) {
    size_t eid = generator_->EmitTetrahedron( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badTetrahedra_ += invalid_elements;

  if ( invalid_elements > 0 )
    return false;

  return true;
}

bool CornerPointGrid::ConstructEclipseCell1100( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructPrism( cell, 0, 3, 7, 1, 2, 6 ) ) {
    size_t eid = generator_->EmitPrism( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badPrisms_ += invalid_elements;

  if ( invalid_elements > 0 )
    return false;

  return true;
}

//Eclipse cell type: ECLIPSE_CELL_PYRAMID_2
bool CornerPointGrid::ConstructEclipseCell1101( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructTetrahedronOnFace( cell, 7, 6, 5, generator_->getNodeID( cell, 2 ) ) ) {
    size_t eid = generator_->EmitTetrahedron( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badTetrahedra_ += invalid_elements;

  if ( invalid_elements > 0 )
    return false;

  return true;
}

//Eclipse cell type: ECLIPSE_CELL_PYRAMID_3
bool CornerPointGrid::ConstructEclipseCell1110( ColumnCell&  cell, size_t& i, size_t& j, size_t& k, ColumnCell* cellAbove, ColumnCell* cellBeneath ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructTetrahedronOnFace( cell, 7, 6, 4, generator_->getNodeID( cell, 3 ) ) ) {
    size_t eid = generator_->EmitTetrahedron( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badTetrahedra_ += invalid_elements;

  if ( invalid_elements > 0 ) {
    //this->elementMap.erase(ijk(i, j, k));
    for ( auto eid : new_elements ) {
      //generator_->plist.erase(eid);
      //generator_->fem_types.pop_back();
      //--generator_->elementID;
    }
    return false;
  }

  return true;
}

bool CornerPointGrid::ConstructLineElement( size_t& i, size_t& j, size_t& k, size_t& new_elemt_idx ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  auto it = columns_.find( make_pair( i, j ) );
  if ( it == columns_.end() ) return false;

  generator_->setIJ( i, j );

  Column& Col = it->second;
  if ( k >= Col.cells_.size() ) return false;

  ColumnCell&  cell = Col.cells_[k];

  if ( generator_->ConstructLine( cell ) ) {
    size_t eid = generator_->EmitLine( cell );
    addElementToMap( i, j, k, eid );
    new_elemt_idx = eid;
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badLines_ += invalid_elements;

  if ( invalid_elements > 0 || new_elements.size() == 0 )
    return false;

  return true;
}


bool CornerPointGrid::ConstructLineElement( ColumnCell& cell, size_t& i, size_t& j, size_t& k ) {
  size_t invalid_elements = 0;
  std::vector<uint32_t> new_elements;

  if ( generator_->ConstructLine( cell ) ) {
    size_t eid = generator_->EmitLine( cell );
    addElementToMap( i, j, k, eid );
    new_elements.push_back( eid );
  }
  else {
    invalid_elements++;
  }

  badLines_ += invalid_elements;

  if ( invalid_elements > 0 )
    return false;

  return true;
}



/**
MASTER METHOD for the creation of corner point grids
from the original cell-centered grids read by the EclipseInterface.

The PolygonGridManager is used to construct the pillars.

@todo this is confused! - why should the corner point grid be responsible for building the CSMP model? - code should be in the EclipseModel

@todo regions are not recognised properly.
*/
void CornerPointGrid::CreateModel( const std::string&     model_name,
                                   csmp::VSet<3U>&        vset,
                                   csmp::ModelTopology&   model_topology,
                                   const std::vector<double>& zcorn,
                                   std::set<std::string>& regions,
                                   std::set<std::string>& faults,
                                   std::set<std::string>&  wells,
                                   bool tetra_mesh, bool exclude_inactive_cells )
{
  // 1. Initialise grid specs
  InitializeGridSpecs();

  // 2. Construct pillars and columns
  ConstructPillarsAndColumns( zcorn );

  // 3. Construct FEs from columns and add to Vset
  ConstructFiniteElementsFromColumns( vset );
}


void CornerPointGrid::addElementToMap( size_t i, size_t j, size_t k, size_t elementID )
{
  this->elementMap.emplace( ijk( i, j, k ), elementID );
}


template<class VarType>
void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&             vset,
                                           const std::vector<VarType>& prop_data,
                                           const std::string&          prop_name,
                                           const csmp::PLACEMENT&      prop_place ) const
{
  /// correcting data cell id's due to existance of embedded cells
  VarType var;
  var = 0.;
  std::vector<VarType> cell_data( vset.Elements(), var );
  for ( auto& entry : elementMap ) {
    auto coord = entry.first;
    size_t idx = coord.i + coord.j * NX_ + coord.k * NX_x_NY_;
    cell_data[entry.second] = prop_data[idx];
  }
  if ( prop_place == csmp::NODE ) {
    // Add nodal data to vset
    std::vector<VarType> nodal_data;
    // TODO: why node property?
    extrapolateElementToNodeProperty<3U, VarType>( vset, cell_data, nodal_data );
    //csmp::FEM_Data<VarType> property_values( prop_place, nodal_data );
    const size_t array_length = (var.Size() > 9U) ? var.Size() : 0U;
    assert( array_length == 1 );
    PropertyData property_values( prop_place, VarType::VariableType, 3U, 0U );
    property_values.Reserve( nodal_data.size() );
    for ( const auto& it : nodal_data ) pushBack( property_values, it );
    vset.AddData( prop_name.c_str(), property_values );
  }
  else {
    // Add cell data to vset
    // csmp::FEM_Data<VarType> property_values( prop_place, cell_data );
    PropertyData property_values( prop_place, VarType::VariableType, 3U, 0U );
    property_values.Reserve( cell_data.size() );
    for ( const auto& it : cell_data ) pushBack( property_values, it );
    vset.AddData( prop_name.c_str(), property_values );
  }

} // end WritePropertyToVSet

template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<ScalarVariable>&, const std::string&, const csmp::PLACEMENT& ) const;
template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<VectorVariable<3U> >&, const std::string&, const csmp::PLACEMENT& ) const;
template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<TensorVariable<3U> >&, const std::string&, const csmp::PLACEMENT& ) const;
template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<ArrayVariable>&, const std::string&, const csmp::PLACEMENT& ) const;
template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<FlaggedArrayVariable>&, const std::string&, const csmp::PLACEMENT& ) const;


void CornerPointGrid::Resize( size_t i_pillar_max, size_t j_pillar_max )
{
  pillars_.resize( i_pillar_max * j_pillar_max );
}


/**
accessing the contained pillars
*/
Pillar & CornerPointGrid::operator()( size_t i, size_t j )
{
  assert( i <= NX_ );
  assert( j <= NY_ );

  return pillars_[j * (NX_ + 1) + i];
}

// CELL CENTERED GRID

CellCenteredGrid::CellCenteredGrid()
  :dx_( 3U )
{
}


CellCenteredGrid::~CellCenteredGrid()
{
}


void CellCenteredGrid::Clear()
{
  /// grid
  tops_.clear();
  dx_.clear();
  dy_.clear();
  dz_.clear();
}


size_t CellCenteredGrid::GetNumCells() const
{
  return NX_*NY_*NZ_;
}


std::vector<csmp::ScalarVariable>& CellCenteredGrid::GetCellDepths()
{
  return tops_;
}


std::vector<ScalarVariable>& CellCenteredGrid::GetCellSizes( size_t i )
{
  if ( i == 0 )
    return dx_;
  if ( i == 1 )
    return dy_;
  return dz_;
}


void CellCenteredGrid::AssignDimensionX( size_t NX )
{
  NX_ = NX;
}


void CellCenteredGrid::AssignDimensionY( size_t NY )
{
  NY_ = NY;
}


void CellCenteredGrid::AssignDimensionZ( size_t NZ )
{
  NZ_ = NZ;
}

} // end namespace csmp
