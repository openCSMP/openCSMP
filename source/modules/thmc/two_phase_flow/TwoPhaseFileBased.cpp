#include "TwoPhaseFileBased.h"
#include "PropertyDatabase.h"
#include "CSMP_mathUtilities.h"
#include "FiniteVolumeStencil.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
TwoPhaseFileBased<dim>::TwoPhaseFileBased( const char* fileName )
 : GRAVITY_(9.80665),
   LOWER_LIMIT_(1.0e-30),
   catchPhrase_( "Pc" ),
   MAX_CAPILLARY_PRESSURE_(1.0e+7),
   MAX_CAPILLARY_PRESSURE_SLOPE_(1.0e+5)
 {
    ReadFile(fileName);
 }

template<uint32_t dim>
TwoPhaseFileBased<dim>::TwoPhaseFileBased()
 : GRAVITY_(9.80665),
   LOWER_LIMIT_(1.0e-30),
   MAX_CAPILLARY_PRESSURE_(1.0e+7),
   MAX_CAPILLARY_PRESSURE_SLOPE_(1.0e+5)
 {
 }

template<uint32_t dim>
TwoPhaseFileBased<dim>::TwoPhaseFileBased( const PropertyDatabase<dim>& database, 
                                           const std::vector<double>& seff,
                                           const std::vector<double>& krn,
                                           const std::vector<double>& krw,
                                           const std::vector<double>& pc  )
 : GRAVITY_(9.80665),
   catchPhrase_( "Pc" ),
   LOWER_LIMIT_(1.0e-30),
   MAX_CAPILLARY_PRESSURE_(1.0e+7),
   MAX_CAPILLARY_PRESSURE_SLOPE_(1.0e+5), // this is a meaningful limit which does not make matrix singular
   TwoPhaseModel<dim>( database, "permeability",
                       "viscosity oil", "viscosity water",
                       "density oil", "density water", "saturation water",
                       "residual saturation oil",
                       "residual saturation water" ),
   seff_( seff ), krw_( krw ), krn_( krn ), pc_( pc )
  {
    krwCurve_.Initialize( seff_, krw_, 0.01, 0.01 );
    krnCurve_.Initialize( seff_, krn_, 0.01, 0.01 );
    pcCurve_.Initialize( seff_, pc_, 0.01, 0.01 );
  }


template<uint32_t dim>
TwoPhaseFileBased<dim>::TwoPhaseFileBased( const PropertyDatabase<dim>& database, const char* fileName )
 : GRAVITY_(9.80665),
   catchPhrase_( "Pc" ),
   LOWER_LIMIT_(1.0e-30),
   MAX_CAPILLARY_PRESSURE_(1.0e+7),
   MAX_CAPILLARY_PRESSURE_SLOPE_(1.0e+5), // this is a meaningful limit which does not make matrix singular
   TwoPhaseModel<dim>( database, "permeability",
                       "viscosity oil", "viscosity water",
                       "density oil", "density water", "saturation water",
                       "residual saturation non-wetting phase",
                       "residual saturation wetting phase" )
 {
   ReadFile(fileName);
 }


template<uint32_t dim>
void TwoPhaseFileBased<dim>::ReadFile( const char* fileName )
 {
   relpermFile_.open( fileName, ios::in | ios::binary);
   if( !relpermFile_ ) {
     throw csmp::Exception( FATAL_ERROR, "TwoPhaseFileBased",
                            "Could not open relative permeability file!",
                            "check input file");
       }

   readData();
   cout << "\n---------------------------------------------------------";
   cout << "\nFile based two phase model setup based on following input:\n";
   cout << "\nSeff\tKrw\tKrn\tPc\n";
   for( uint32_t i=0; i<seff_.size(); ++i )
     cout << seff_.at( i ) << "\t" << krw_.at( i ) << "\t" << krn_.at( i ) << "\t" << pc_.at( i ) << endl;
   cout << "---------------------------------------------------------\n";

   krwCurve_.Initialize( seff_, krw_, 0.01, 0.01 );
   krnCurve_.Initialize( seff_, krn_, 0.01, 0.01 );
   pcCurve_.Initialize( seff_, pc_, 0.01, 0.01 );
 }


template<uint32_t dim>
TwoPhaseFileBased<dim>::~TwoPhaseFileBased()
 {
 }



template<uint32_t dim>
void TwoPhaseFileBased<dim>::Initialize( const Element<dim>& e )
 {
  TwoPhaseModel<dim>::swr_ = e.Read( TwoPhaseModel<dim>::swr_key_ );
  TwoPhaseModel<dim>::snr_ = e.Read( TwoPhaseModel<dim>::snr_key_ );

  if( TwoPhaseModel<dim>::tensor_permeability_){
      e.Read( TwoPhaseModel<dim>::perm_key_, TwoPhaseModel<dim>::K_);
      TwoPhaseModel<dim>::k_ = TwoPhaseModel<dim>::K_.Trace()/static_cast<double>(dim);
  }else{
      TwoPhaseModel<dim>::k_ = e.Read( TwoPhaseModel<dim>::perm_key_ );
      TwoPhaseModel<dim>::K_.operator=( VectorVariable<dim>(PLAIN, TwoPhaseModel<dim>::k_ ) );
  }


 } // end Initialize


template<uint32_t dim>
double TwoPhaseFileBased<dim>::krw_Phase() const
 {
  double se = TwoPhaseModel<dim>::seff_;

  if ( se <= 0.) return static_cast<double>(0.);
  if ( se >= 1.) return static_cast<double>(1.);

   size_t i(0);
    for( i = 1; i < seff_.size(); ++i )
      if( seff_[i] > se )
        break;
    return (krw_[i-1]-krw_[i])/(seff_[i-1]-seff_[i]) * (se - seff_[i-1])+krw_[i-1] ;
 }


template<uint32_t dim>
double TwoPhaseFileBased<dim>::krn_Phase() const
 {
  double se = TwoPhaseModel<dim>::seff_;

  if ( TwoPhaseModel<dim>::seff_ <= 0.) return static_cast<double>(1.);
  if ( TwoPhaseModel<dim>::seff_ >= 1.) return static_cast<double>(0.);

  size_t i(0);
  for( i = 1; i < seff_.size(); ++i )
    if( seff_[i] > se )
      break;
  return (krn_[i-1]-krn_[i])/(seff_[i-1]-seff_[i]) * (se - seff_[i-1])+krn_[i-1] ;
 }

template<uint32_t dim>
double TwoPhaseFileBased<dim>::dkrwds_Phase() const
 {
  double se = TwoPhaseModel<dim>::seff_;

  if ( se <= 0.) return static_cast<double>(0.);
  if ( se >= 1.) return static_cast<double>(0.);

  const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

   size_t i(0);
    for( i = 1; i < seff_.size(); ++i )
      if( seff_[i] > se )
        break;
    return (krw_[i-1]-krw_[i])/(seff_[i-1]-seff_[i])*seff_mult;
 }


template<uint32_t dim>
double TwoPhaseFileBased<dim>::dkrnds_Phase() const
 {
  double se = TwoPhaseModel<dim>::seff_;

  if ( TwoPhaseModel<dim>::seff_ <= 0.) return static_cast<double>(0.);
  if ( TwoPhaseModel<dim>::seff_ >= 1.) return static_cast<double>(0.);

  const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

  size_t i(0);
  for( i = 1; i < seff_.size(); ++i )
    if( seff_[i] > se )
      break;
  return (krn_[i-1]-krn_[i])/(seff_[i-1]-seff_[i])*seff_mult;

 }


template<uint32_t dim>
double TwoPhaseFileBased<dim>::pc_Phase( ) const
{
  double se = TwoPhaseModel<dim>::seff_;

  // only for the min saturation of water precautions are needed
  // the actual saturation is used instead of the effective saturation
  if ( se == 0. ) return TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
  if ( se == 1. ) return static_cast<double>(0.);

  size_t i(0);
  for( i = 1; i < seff_.size(); ++i )
    if( seff_[i] > se )
      break;
  return (pc_[i-1]-pc_[i])/(seff_[i-1]-seff_[i]) * (se - seff_[i-1])+pc_[i-1] ;
}

template<uint32_t dim>
double TwoPhaseFileBased<dim>::dpcds_Phase( ) const
{
  double se = TwoPhaseModel<dim>::seff_;

  if ( se == 0. ) return -MAX_CAPILLARY_PRESSURE_SLOPE_;
  if ( se == 1. ) return -MAX_CAPILLARY_PRESSURE_SLOPE_;

  const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

  size_t i(0);
  for( i = 1; i < seff_.size(); ++i )
    if( seff_[i] > se )
      break;
  return (pc_[i-1]-pc_[i])/(seff_[i-1]-seff_[i])*seff_mult;

}


template<uint32_t dim>
double TwoPhaseFileBased<dim>::Sw_Phase( double pc ) const
{

  // only for the min saturation of water precautions are needed
  // the actual saturation is used instead of the effective saturation
  if ( pc == TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ ){
      TwoPhaseModel<dim>::seff_ = 0.0;
      return TwoPhaseModel<dim>::SeffToSw();
  }
  if ( pc == 0.0 ){
      TwoPhaseModel<dim>::seff_ = 1.0;
      return TwoPhaseModel<dim>::SeffToSw();
  }

  size_t i(0);
  for( i = 1; i < pc_.size() ; ++i )
      if( ( pc <=pc_[i-1] ) && (pc > pc_[i]))
          break;

  TwoPhaseModel<dim>::seff_ = (pc - pc_[i-1])/(pc_[i-1]-pc_[i])*(seff_[i-1]-seff_[i]) + seff_[i-1];
  return TwoPhaseModel<dim>::SeffToSw();

}

template<uint32_t dim>
double TwoPhaseFileBased<dim>::dsdpc_Phase( double pc ) const
{

  if ( pc == TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ )
      return -1.0/MAX_CAPILLARY_PRESSURE_SLOPE_;
  if ( pc == 0.0 )
      return -1.0/MAX_CAPILLARY_PRESSURE_SLOPE_;

  const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

  size_t i(0);
  for( i = 1; i < pc_.size() ; ++i )
      if( ( pc <=pc_[i-1] ) && (pc > pc_[i]))
          break;

  return (seff_[i-1]-seff_[i])/(pc_[i-1]-pc_[i])/seff_mult;

}

template<uint32_t dim>
double TwoPhaseFileBased<dim>::MaxFractionalFlowDerivative() const
 {
    return static_cast<double>(7.); // as computed with dfds_Phase method
 }

template<uint32_t dim>
double TwoPhaseFileBased<dim>::dfds() const
{  
    return TwoPhaseModel<dim>::dfds_numerical();
}  // end dfdS_Phase


template<uint32_t dim>
double TwoPhaseFileBased<dim>::dGds( ) const
{
    return TwoPhaseModel<dim>::dGds_numerical();
}

/// internal helper function to find position in a given file
template<uint32_t dim>
streampos  TwoPhaseFileBased<dim>::findPosition( std::ifstream& file ) const{

  string catchPhrase;

  while( file >> catchPhrase )
    if( catchPhrase == catchPhrase_ )
      return file.tellg();

  return 0;
}

/// reads line values to vector
template<uint32_t dim>
int32_t TwoPhaseFileBased<dim>::readData(){

  int32_t           count = 0;
  double          cache;
  const streampos position( findPosition( relpermFile_ ) );

  relpermFile_.seekg( position, ios::beg );
     while( !relpermFile_.eof() ){
       relpermFile_ >> cache;
       seff_.push_back( cache );
       relpermFile_ >> cache;
       krw_.push_back( cache );
       relpermFile_ >> cache;
       krn_.push_back( cache );
       relpermFile_ >> cache;
       pc_.push_back( cache );
       ++count;
     }

  return count;
}

template<uint32_t dim>
int32_t TwoPhaseFileBased<dim>::writeData(){

  ofstream file( "TwoPhaseFileBased.txt", ios::out );
  file << "Seff\tKrw\tKrn\tPc\n";
  for( auto i = 0; i<=10; ++i )
    file << i*0.1 << "\t" << krwCurve_.Value( i*0.1 ) << "\t" << krnCurve_.Value( i*0.1 ) << "\t" << pcCurve_.Value( i*0.1 ) << endl;

  file.close();

  return 0;
}



template class TwoPhaseFileBased<1U>;
template class TwoPhaseFileBased<2U>;
template class TwoPhaseFileBased<3U>;

} // end namespace csp









