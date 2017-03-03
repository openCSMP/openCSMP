#include "TwoPhaseModel_Test.h"

#include <vector>
#include <iterator>
#include "Model.h"
#include "TwoPhaseModel.h"
#include "VTU_Interface.h"

#include "TwoPhaseModel.h"

#include "BrooksCorey.h"
#include "FourarLenormand.h"
#include "VanGenuchten.h"
#include "TwoPhaseFileBased.h"
#include "Experimental2PhaseModel.h"
#include "ExperimentalRT.h"
#include "FractureMatrixUpscaled.h"
#include "BrooksCoreyWithHysteresis.h"
#include "LinearTwoPhaseModel.h"

using namespace std;

namespace csmp {


TwoPhaseModel_Test::TwoPhaseModel_Test( Model<1U>* model,
                                        TwoPhaseModel<1U>* twoPhaseModel,
                                        const char* testName,
                                        const char* krn,
                                        const char* krw,
                                        const char* pc,
                                        size_t upToElementNumber )
 : model_( model ),
   twoPhaseModel_( twoPhaseModel ),
   upToElementNumber_( upToElementNumber ),
   krw_key_(model->Database().StorageKey(krw)),
   krn_key_(model->Database().StorageKey(krn)),
   pc_key_(model->Database().StorageKey(pc))
{

    setName( testName );

    vtuOutputProps_.push_back( krn );
    vtuOutputProps_.push_back( krw );
    vtuOutputProps_.push_back( pc );


}


TwoPhaseModel_Test::~TwoPhaseModel_Test( )
{
    if( twoPhaseModel_!=NULL)
        delete twoPhaseModel_;
}



// TODO file output, test for more than just nan(test criteria)
void TwoPhaseModel_Test::run()
{
  cout.setf(ios_base::scientific);

  cout << "\nTwoPhaseModel_Test: " << getName() << endl;

  Region<1>& model = model_->Region( "Model" );
  for ( vector<Element<1U>*>::iterator it = model.ElementsBegin();
        it != model.ElementsEnd(); ++it )
  {
  if( upToElementNumber_ != 0 && distance( model.ElementsBegin(), it ) > upToElementNumber_  )
    break;
    // setting up the relative permeability model
    // ---------------------------------------------
    twoPhaseModel_->Initialize( *(*it) );
    twoPhaseModel_->InitializeForNode( *(*it), 1U );
    twoPhaseModel_->EffectiveSaturation();

    const double64 sw( twoPhaseModel_->Saturation(1U) );
    const double64 seff( twoPhaseModel_->EffectiveSaturation() );
    const double64 swr( twoPhaseModel_->Swr() );
    const double64 snr( twoPhaseModel_->Snr() );

    // water saturation & effective water saturation
    // ---------------------------------------------
    _test( !( sw < 0. || sw > 1. ) );
    _test( !( seff < 0. || 1. < seff ) );

    // relative permeability
    // ---------------------
    _test( !( twoPhaseModel_->krw_Phase() < 0. || twoPhaseModel_->krn_Phase() < 0. ) );

    // total mobility
    // --------------

    // fractional flow and its saturation derivative
    // ---------------------------------------------
    _test( !( twoPhaseModel_->f_Phase(2U) < 0. || twoPhaseModel_->f_Phase(2U) > 1. ) );
    //_test( !isnan( twoPhaseModel_->dfds() ) );

    // G and derivative of G-function
    // ------------------------------
    _test( !isnan( twoPhaseModel_->G() ) );
    _test( !isnan( twoPhaseModel_->dGds() ) );

    // 5. capillary pressure and its saturation derivative
    // -------------------------------------------------------------------------------------------
    _test( !isnan( twoPhaseModel_->pc_Phase() ) );
    _test( !( twoPhaseModel_->pc_Phase() > 1.e9 ) );
    _test( !isnan( twoPhaseModel_->dpcds_Phase() ) );


    // the 3 characteristic multipliers: advection, diffusion, gravity
    // -------------------------------------------------------------------------------------------
    _test( !isnan( twoPhaseModel_->AdvectionMultiplier() ) );
    _test( !isnan( twoPhaseModel_->DiffusionMultiplier(2U) ) );
    _test( !isnan( twoPhaseModel_->GravityMultiplier_G() ) );
    _test( !isnan( twoPhaseModel_->GravityMultiplier_dGds() ) );
    _test( !isnan( twoPhaseModel_->CapillaryDiffusionMultiplier() ) );
    _test( !isnan( twoPhaseModel_->ShockSpeed() ) );


    // write data
    twoPhaseModel_->Initialize( *(*it) ); // just in case
    for( size_t i=0U;i<(*it)->Nodes();i++)
    {
        twoPhaseModel_->InitializeForNode( *(*it), i );
        twoPhaseModel_->EffectiveSaturation();
        (*it)->N(i)->Store( krn_key_, makeScalar( PLAIN, twoPhaseModel_->krw_Phase()));
        (*it)->N(i)->Store( krw_key_, makeScalar( PLAIN, twoPhaseModel_->krn_Phase()));
        (*it)->N(i)->Store( pc_key_,  makeScalar( PLAIN, twoPhaseModel_->pc_Phase()));
    }

  }

  VTU_Interface<1U> vtu( *model_ );
  vtu.OmitZeroInFileName( true );
  vtu.OutputDataToVTU(  this->getName().c_str(), vtuOutputProps_,"Model", static_cast<int>(0));


} // run()

} // csmp
