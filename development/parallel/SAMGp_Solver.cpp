#include "SAMGp_Solver.h"
#include "SAMG_Exception.h"
#include "CSP_ErrorHandler.h"
#include "CSP_math.h"

using namespace std;

namespace csp {

template<typename fT, stl_index dim>
  SAMGp_Solver<fT,dim>::SAMGp_Solver() :
    settings_(new SAMG_Settings()),
    u_(0),
    f_(0),
    iscale_(0),
    iu_(0),
    ip_(0),
    ierr_(0),
    output_amg_data_to_text_files_(false)
    {
//    MPI::Group new_group (MPI::COMM_WORLD.Get_group().Excl(0,0));
//    MPI::Intracomm new_comm (MPI::COMM_WORLD.Create(new_group));
//    fortran_communicator = (int)MPI_Comm_c2f(new_comm);
    }
    
    

template<typename fT, stl_index dim>
  SAMGp_Solver<fT,dim>::SAMGp_Solver(SAMG_Settings* settings, SAMGp_CommunicationData<fT,dim>* parameters, bool update_halo) :
    settings_(settings),samgp_parameters_(parameters),update_halo_(update_halo),tag(0),
    u_(0),//    int my_rank = MPI::COMM_WORLD.Get_rank( );
    f_(0),
    iscale_(0),
    iu_(0),
    ip_(0),
    ierr_(0),
    output_amg_data_to_text_files_(false)
    {
      irecdata.resize(samgp_parameters_->ireclist.size());
      isnddata.resize(samgp_parameters_->isndlist.size());
      mpi_req.resize (samgp_parameters_->irankrec.size());
    }



  // SKM addition, default destructor will not work  
template<typename fT, stl_index dim>
  SAMGp_Solver<fT,dim>::~SAMGp_Solver()
   {
/*
        // don't delete settings_, this would create dangling pointer problem
	  delete[] ia_;
	  delete[] ja_;
	  delete[] a_;
	  delete[] u_;
	  delete[] f_;
	  delete[] iscale_;
	  delete[] iu_;
	  delete[] ip_;
          delete[] samgp_parameters_;
*/
   }
   

    
  /*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  int SAMGp_Solver::SetSolverSettings(SAMG_Settings* settings)
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->
Set new solver settings.
<p>
<H4>Input Arguments:</H4><!----------------------------------------------->
Pointer to a SAMG_Settings object.
<p>
<H4>Output Arguments &amp; Return Value</H4><!---------------------------->


<H4>Implementation:</H4><!------------------------------------------------>
The new setting is stored as a pointer variable.
<p>
WARNING: The SAMGp_Solver object DOES NOT delete ANY pointer to SAMG_Settings
objects given to it. The user is therefore solely responsible for guaranteeing
the destruction of SAMG_Settings objects created with new!
<p>
<H4>Application:</H4><!--------------------------------------------------->


<H4>Messages:</H4><!------------------------------------------------------>

<!------------------------------------------------------------------------>
tested: */  
template<typename fT, stl_index dim>
  void SAMGp_Solver<fT,dim>::SetSolverSettings(SAMG_Settings* settings) {
    assert(settings);
    settings_ = settings;
  }
  
template<typename fT, stl_index dim>
  SAMG_Settings* SAMGp_Solver<fT,dim>::GetSolverSettings() const {
    return settings_;
  }
  
template<typename fT, stl_index dim>
  void SAMGp_Solver<fT,dim>::SolveMatrixEquation( SparseMatrix& A, 
                                          vector<csp_float>& b, 
                                          vector<csp_float>& x,
                                          stl_index no_unknowns )
  {
    extern csp::ErrorHandler skm_err;

    nshalo = samgp_parameters_->nshalo();
    npsnd  = samgp_parameters_->npsnd();
    nrhalo = samgp_parameters_->nrhalo();
    nprec  = samgp_parameters_->nprec();

    int new_nnu = static_cast<int>(A.Rows()-nrhalo);
    int new_nsys = static_cast<int>(no_unknowns);

    crmat_.Initialize( A );
    crmat_.RemoveHalo( nrhalo );
    u_.resize (new_nnu+nrhalo);
    f_.resize (new_nnu);
    int new_nna = static_cast<int>(crmat_.a.size());

    if (Verbose()) {
      cout << "\nSAMGp_Solver::SolveMatrixEquation: Sizes NNU(rows=cols): "
           << new_nnu << " NNA(total unknowns): " << new_nna << "   ";
    }


    // set arrays for systems
    if (new_nsys > 1)
     {
      // only rebuild arrays if settings have changed
      if (new_nnu != nnu_ || new_nsys != nsys_)
       {  
        if (new_nnu != nnu_)
         {
	    iu_.resize(new_nnu+nrhalo);
	    ndiu_ = iu_.size(); // update the size indicator for iu vector
	 }
    
        int k = 0; 
        // coupled systems with solution vector [x1, y1, x2, y2, x3, y3, ...., xn, yn]
        if ( settings_->UsePointBasedApproach() )
         {
          for ( int j = 0; j < new_nsys; j++ )
	    for ( int i = 0; i < (iu_.size()/new_nsys); i++ )
	       iu_[k++] = j+1;
         }
        else
         {
          // coupled systems with solution vector [x1, x2, x3,..., xn; y1, y2, y3, ...., yn]
          for ( int i=0; i<(iu_.size()/new_nsys); i++ )
	    for ( int j=0; j<new_nsys; j++ )
	      iu_[k++] = j+1;        
         }
       }
      // use point information as well
      if ( settings_->UsePointBasedApproach() )
       {
	// only rebuild arrays if settings have changed
	if (new_nnu != nnu_ || new_nsys != nsys_)
         {
	  if (new_nnu != nnu_)
           {
            ip_.resize(new_nnu+nrhalo);
	    ndip_ = ip_.size(); // resizing the ip vector size indicator
	   }

	  int k = 0;
	  for ( int i=0; i<(ip_.size()/new_nsys); i++ )
	    for ( int j=0; j<new_nsys; j++ )
	      // must be the corresponding node number 
	      ip_[k++] = i+1;
	 }
       }
    }

   else
    { // new_nsys < 1 or new_nsys=1
      if (ndiu_ != 1)
       {
         iu_.resize(1);
       }
      if (ndip_ != 1)
       {
         ip_.resize(1);
       }
      ndiu_ = 1;
      ndip_ = 1;
    }
    
    // To be done for scalar systems as well
    if (new_nsys != nsys_)
     {
       iscale_.resize(new_nsys);
       // putting zero values into this array switches the scaling off
       fill( iscale_.begin(), iscale_.end(), 0 );
     }
    
    // assign new values to actual ones
    nsys_ = new_nsys;
    nna_  = new_nna;
    nnu_  = new_nnu;
  
    // test for NULL pointers
    if( crmat_.a.empty() || crmat_.ja.empty() || crmat_.ia.empty() || u_.empty() || f_.empty() )
      skm_err.notice(FATAL_ERROR, 
                     "SAMGp_Solver::SolveMatrixEquation", 
                     "Unable to allocate required memory for transfer arrays");

    // initial guess for the solution vector
    for ( stl_index i=0; i<u_.size(); i++ )
      u_[i] = x[i]; 
    // right-hand side
    for ( stl_index i=0; i<f_.size(); i++ )
      f_[i] = b[i]; 

    if (Verbose()) {
      cout << "\nSAMG_Adaptor::SolveSAMG: Calling SAMG..." << endl;
    }

    // 3. initialize settings
    // 3.1 Output parameter
    res_in_         = -1.;
    res_out_        = -1.;
    int ncyc_done = 0;
    
    // 3.2 Paramater from SAMG_Settings object
    int    matrix = settings_->Get_matrix();
    int    nsolve = settings_->Get_nsolve();
    int    ifirst = settings_->Get_ifirst();
    double eps = (ifirst == 0 ? settings_->Get_rel_eps() * std::sqrt(vector_norm2(b, b)) : settings_->Get_eps());
    int    ncyc = settings_->Get_ncyc();
    int    iswtch = settings_->Get_iswtch();
    double a_cmplx = settings_->Get_a_cmplx();
    double g_cmplx = settings_->Get_g_cmplx();
    double p_cmplx = settings_->Get_p_cmplx();
    double w_avrge = settings_->Get_w_avrge();
    double chktol = settings_->Get_chktol();
    int    idump = settings_->Get_idump();
    int    iout = settings_->Get_iout();
    int    ncg = settings_->Get_ncg();

    if (settings_->ExplicitSecondary())
       SAMG_SET_NCG(&ncg);
    else
       SAMG_RESET_SECONDARY();


    int fortran_communicator = (int)MPI_Comm_c2f(MPI::COMM_WORLD);

/*
// ONLY DEBUGGING:
assert(u_.size()==x.size());
assert(x.size()==b.size());
assert(A.Rows()==u_.size());
assert(A.Rows()-nrhalo==f_.size());

// SAMGp DEBUGGING SETTINGS:
int my_rank = MPI::COMM_WORLD.Get_rank( );
int l,ioout,iunf;
iout = 31, l=1, iunf=0, ioout = 8+my_rank, idump = 8;
SAMG_SET_LOGIO(&ioout);
*/
// HERE CALL SAMGP:

cout <<"\nRank "<<MPI::COMM_WORLD.Get_rank()<<": calling SAMGP now........"<< endl;

    // 2. writing SAMG solver input/output data to file
    // -----------------------------------------
    if ( output_amg_data_to_text_files_ )
     {
       int my_rank = MPI::COMM_WORLD.Get_rank();
       string out("csp_solution_data_");
       char rank[200];
       sprintf(rank, "%ld", static_cast<long>(my_rank));
       out += rank;

       cout << "\nSAMGp_Solver): SAMG TEXT FILE OUTPUT HAS BEEN ENABLED ! "
       << "Watch for '.frm', '.amg', '.rhs' and perhaps '.iu' and 'ip' files that will be written."<< endl;
       Write_SAMG_TextInputFile( out.c_str() );
     }


    SAMGP( &nnu_, &nna_, &nsys_,
	   &crmat_.ia[0], &crmat_.ja[0], &crmat_.a[0], &f_[0], &u_[0],
	   &iu_[0], &ndiu_, &ip_[0], &ndip_, &matrix, &iscale_[0],
	   &res_in_, &res_out_, &ncyc_done, &ierr_,
	   &nsolve, &ifirst, &eps, &ncyc, &iswtch,
	   &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
	   &chktol, &idump, &iout,
           &nshalo, &npsnd, &samgp_parameters_->iranksnd[0], &samgp_parameters_->ipts[0], &samgp_parameters_->isndlist[0],
	   &nrhalo, &nprec, &samgp_parameters_->irankrec[0], &samgp_parameters_->iptr[0], &samgp_parameters_->ireclist[0],
	   &fortran_communicator
           );

    // treat errors with exception
    if (ierr_ > 0)
     {
      skm_err.notice(CSP_ERROR,
                     "SAMGp_Solver::SolveMatrixEquation",
                     "SAMGp solver returned with an error; error code:");
      cerr << ierr_ << endl;
     }
    else if (ierr_ < 0)
     {
      skm_err.notice(WARNING,
                    "SAMGp_Solver::SolveMatrixEquation",
                    "SAMGp solver returned with a warning; code:");
      cerr << ierr_ << endl;
     }

    CheckConvergence(eps);

    if (Verbose())
      cout << "\tpassing solution values to 'x' vector<csp_float>..."; 

    // communication fo outerhalo data
    if(update_halo_)
     CommunicateOuterhalo();

    // solver returned ok so lets place contents back into x
    for ( stl_index i=0; i<u_.size(); i++ ) x[i] = u_[i];  

    if (Verbose())
      cout << "OK\n";

/* 
   // 2. writing SAMG solver input/output data to file
    // -----------------------------------------
    if ( output_amg_data_to_text_files_ )
     {
       int my_rank = MPI::COMM_WORLD.Get_rank();
       string out("csp_solution_data_");
       char rank[200];
       sprintf(rank, "%ld", static_cast<long>(my_rank));
       out += rank;

       cout << "\nSAMGp_Solver): SAMG TEXT FILE OUTPUT HAS BEEN ENABLED ! "
       << "Watch for '.frm', '.amg', '.rhs' and perhaps '.iu' and 'ip' files that will be written."<< endl;
       Write_SAMG_TextInputFile( out.c_str() );
     }
*/


  }

template<typename fT, stl_index dim>
void  SAMGp_Solver<fT,dim>::Write_SAMG_TextOutput(bool write) {
    output_amg_data_to_text_files_ = write;
 }

/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  void SAMGp_Solver::OutputVectors( )
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

<p> 

<H4>Input Arguments:</H4><!----------------------------------------------->


<H4>Output Arguments &amp; Return Value</H4><!---------------------------->


<H4>Implementation:</H4><!------------------------------------------------>


<H4>Application:</H4><!--------------------------------------------------->


<H4>Messages:</H4><!------------------------------------------------------>

<!------------------------------------------------------------------------>
*/
template<typename fT, stl_index dim>
void SAMGp_Solver<fT,dim>::OutputVectors()
  {
    int my_rank = MPI::COMM_WORLD.Get_rank();
    string out("dual_output_");
    char rank[200];
    sprintf(rank, "%ld", static_cast<long>(my_rank));
    out += rank;
    out += ".txt";
    ofstream ofs(out.c_str());
    ofs.setf(ios::scientific);
    long prec = ofs.precision(15);
    ofs.precision(prec);

    // output u and f to screen
    ofs<<"\n\nu:"<< endl;
    for (unsigned int i=0;i!=u_.size();i++)
      ofs<<i+1<<", "<<static_cast<long double>(u_[i])<< endl;

    ofs<<"\n\nf:"<< endl;
    for (unsigned int i=0;i!=f_.size();i++)
      ofs<<i+1<<", "<<static_cast<long double>(f_[i])<< endl;
  }



/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  void SAMGp_Solver::CommunicateOuterhalo( )
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

<p> 

<H4>Input Arguments:</H4><!----------------------------------------------->


<H4>Output Arguments &amp; Return Value</H4><!---------------------------->


<H4>Implementation:</H4><!------------------------------------------------>


<H4>Application:</H4><!--------------------------------------------------->


<H4>Messages:</H4><!------------------------------------------------------>

<!------------------------------------------------------------------------>
*/

template<typename fT, stl_index dim>
void SAMGp_Solver<fT,dim>::CommunicateOuterhalo()
 {
    // All process post non-blocking receive to all neighbors
    for (unsigned int i=0;i!=samgp_parameters_->irankrec.size();i++)
      mpi_req[i] = MPI::COMM_WORLD.Irecv(&irecdata[samgp_parameters_->iptr[i]-1],
                                         (samgp_parameters_->iptr[i+1] - samgp_parameters_->iptr[i]),
                                         MPI::DOUBLE,
                                         samgp_parameters_->irankrec[i],
                                         tag);

    // Now update data vector with data to be send:
    for (unsigned int i=0;i!=samgp_parameters_->isndlist.size();i++)
      isnddata[i]=u_[samgp_parameters_->isndlist[i]-1];

    // All process post non-blocking send to all neighbors
    for (unsigned int i=0;i!=samgp_parameters_->iranksnd.size();i++)
      MPI::COMM_WORLD.Isend(&isnddata[samgp_parameters_->ipts[i]-1],
                            (samgp_parameters_->ipts[i+1] -  samgp_parameters_->ipts[i]),
                            MPI::DOUBLE,
                            samgp_parameters_->iranksnd[i],
                            tag);

    // Check if irecdata is received
    for (unsigned int i=0;i!=mpi_req.size();i++)
      mpi_req[i].Wait();

    // Now update variables with received data
    for (unsigned int i=0;i!=irecdata.size();i++)
      u_[samgp_parameters_->ireclist[i]-1]=irecdata[i];
 }


/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  bool SAMG_Adaptor::Write_SAMG_TextInputFile( const char* filename ) const
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

Outputs the input SparseMatrix into the desired format for the 
AMG solver's file-based interface. This entails the generation of
a format file (*.frm), <p> 

<H4>Input Arguments:</H4><!----------------------------------------------->


<H4>Output Arguments &amp; Return Value</H4><!---------------------------->


<H4>Implementation:</H4><!------------------------------------------------>


<H4>Application:</H4><!--------------------------------------------------->


<H4>Messages:</H4><!------------------------------------------------------>

<!------------------------------------------------------------------------>
tested: */
template<typename fT, stl_index dim>
bool  SAMGp_Solver<fT,dim>::Write_SAMG_TextInputFile( const char* file ) const
 {
    string    out_file(file);
    bool      complete_output(true);
    
    out_file += ".frm";
    
    ofstream  ofs(out_file.c_str());
    
    // -----------------------------------------------------
    // .frm  type of data format and sizes of pointer arrays
    // -----------------------------------------------------
    // Hierbei sind nna,nnu,matrix,nsys und der Vektor iscale
    // vom samg(...) interface her bekannt.
    // NPNT ist 1 oder 0 je nachdem, ob Punktinformationen vorhanden
    // sind oder nicht (d.h. ob der Vektor pi_ip die Laenge m_ndip=nnu hat und mit
    // sinnvollen Punktnummern gefuellt ist, oder ob pi_iu ein dummy vector
    // der Laenge m_ndip=1 ist).
    ofs << out_file <<"\nCSP output file created for SAMG test."<< endl;
    ofs <<"f         4"<< endl;
    // in the following matrix is a composite parameter which indicates whether
    // the solution matrix is symmetric (1vs2) and whether the rowsums are 0 (1vs2)
    ofs <<"# NNA  NNU   NSYS  NPNT"<< endl;
    ofs << nna_ <<"  "<< nnu_  <<"  "<< nsys_ <<"  ";
    // if point information is present 
    if ( ndip_ > 1 ) ofs << 1 << endl;
    else             ofs << 0 << endl;
    ofs <<"# ISCALE(1...NSYS)"<< endl;
    // no-scaling of any of the unknowns shall occur
    for ( int i=0; i<nsys_; i++ ) ofs << 0 << endl; 
    ofs.close();
    cout <<"\nSAMG_Adaptor::Write_SAMG_TextInputFile: file '"<< out_file; 
    cout <<" written successfully."<< endl;
    
    
    // -----------------------------------------------------
    // .amg  ia(1...nnu+1), ja(1...nna), a(1...nna) condensed
    //       matrix storage indices and matrix
    // -----------------------------------------------------
    // enthaelt hintereinanderweg die drei Vektoren ia,ja,a
    // (bei Dir pi_ia, pi_ja, pd_a). Immer nur 1 Eintrag pro Zeile,
    // dh. nnu+1 + nna +nna Zeilen.
    assert( !crmat_.ia.empty() );
    assert( !crmat_.ja.empty() );
    assert( !crmat_.a.empty() );
    out_file  = file;
    out_file += ".amg";
    ofs.open(out_file.c_str());
    // this order is O.K. because the indices are not printed
    for ( stl_index i=0; i<nnu_+1; i++ ) ofs << crmat_.ia[i] << endl;
    for ( stl_index i=0; i<nna_;   i++ ) ofs << crmat_.ja[i] << endl;
    ofs.setf(ios::scientific);
    long prec = ofs.precision(15);
    for ( stl_index i=0; i<nna_;   i++ ) ofs << crmat_.a[i] << endl;
    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
    cout <<"\nSAMG_Adaptor::Write_SAMG_TextInputFile: file '"<< out_file; 
    cout <<" written successfully."<< endl;

    // -----------------------------------------------------
    // .rhs		pd_f(1...nnu)
    // -----------------------------------------------------
    // wenn vorhanden enthaelt die right hand side (d.h. nnu Zeilen).
    assert( !f_.empty() );
    out_file  = file;
    out_file += ".rhs";
    ofs.open(out_file.c_str());
    ofs.setf(ios::scientific);
    prec = ofs.precision(15);
    for ( stl_index i=0; i<nnu_; i++ ) ofs << f_[i] << endl;
    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
    cout <<"\nSAMG_Adaptor::Write_SAMG_TextInputFile: file '"<< out_file; 
    cout <<" written successfully."<< endl;

    // -----------------------------------------------------
    // .lhs		pd_f(1...nnu)
    // -----------------------------------------------------
    // wenn vorhanden enthaelt die right hand side (d.h. nnu Zeilen).
    assert( !u_.empty() );
    out_file  = file;
    out_file += ".lhs";
    ofs.open(out_file.c_str());
    ofs.setf(ios::scientific);
    prec = ofs.precision(15);
    for ( stl_index i=0; i<u_.size(); i++ ) ofs << u_[i] << endl;
    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
    cout <<"\nSAMG_Adaptor::Write_SAMG_TextInputFile: file '"<< out_file; 
    cout <<" written successfully."<< endl;

    // -----------------------------------------------------
    // .a		pd_f(1...nnu)
    // -----------------------------------------------------
    // wenn vorhanden enthaelt die right hand side (d.h. nnu Zeilen).
    assert( !crmat_.a.empty() );
    out_file  = file;
    out_file += ".a";
    ofs.open(out_file.c_str());
    ofs.setf(ios::scientific);
    prec = ofs.precision(15);
    stl_index i=0;
    for ( stl_index row=0; row!=nnu_; row++ )
      {
        ofs <<"\n"<<(row+1)<<"\t"; // node number
        for (stl_index j=crmat_.ia[row];j!=crmat_.ia[row+1];j++,i++)
          ofs <<crmat_.a[i]<<"\t["<<crmat_.ja[i]<<"]\t";
      }

    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
    cout <<"\nSAMG_Adaptor::Write_SAMG_TextInputFile: file '"<< out_file; 
    cout <<" written successfully."<< endl;


    // -----------------------------------------------------
    // .iu   iu(1...nnu)	for the ith variable, iu(i) 
    //                      is the identifier of the physical
    //                   	unknown (1...nsys)
    // -----------------------------------------------------
    // notwendig im Falle eines PDE-Systems
    // enthaelt den Vektor iu (bei Dir pi_iu),
    // also (m_)ndiu Zeilen
    // Falls die Matrix nur von einer skalaren PDE
    // herkommt, fallen die Dateien <name>.iu und <name.ip> weg.
    if ( nsys_ > 1 && ndiu_ > 1 ) {
         assert( !iu_.empty() );
	     out_file  = file;
	     out_file += ".iu";
	     ofs.open(out_file.c_str());
	     for ( stl_index i=0; i<ndiu_; i++ ) ofs << iu_[i] << endl;
	     ofs.close();
         cout <<"\nSAMG_Adaptor::Write_SAMG_TextInputFile: file '"<< out_file; 
         cout <<" written successfully."<< endl;
      }

    // -----------------------------------------------------
    // .ip	 ip(1...nnu)	for the ith variable, ip(i) is the identifier of the node where
    //						variable i is located (this can only be used if the variables are
    //						ordered point-wise and, for every point, unknown wise and ascending.
    // -----------------------------------------------------
    // enthaelt den Vektor ip (bei Dir pi_ip),
    // also (m_)ndip Zeilen
    if ( nsys_ > 1 && ndip_ > 1 ) {
         assert( !ip_.empty() );
	     out_file  = file;
	     out_file += ".ip";
	     ofs.open(out_file.c_str());
	     for ( stl_index i=0; i<ndip_; i++ ) ofs << ip_[i] << endl;
	     ofs.close();
       cout <<"\nSAMG_Adaptor::Write_SAMG_TextInputFile: file '"<< out_file; 
       cout <<" written successfully."<< endl;
    }
      
    return complete_output;
    
 } // end Write_SAMG_TextInputFile
 
 
 
 
template<typename fT, stl_index dim>
 csp_float SAMGp_Solver<fT,dim>::LastSolverResidual() const {
  return res_out_;
 }




template<typename fT, stl_index dim>
void SAMGp_Solver<fT,dim>::CheckConvergence(double64 eps) 
 {
    if ( eps <= 1e-25 ) return; // if automatic conversion criterion was set to 0.
 
    extern csp::ErrorHandler skm_err;
    
    if (settings_->Get_ifirst() == 0) {
      if (res_out_ > eps ) {
        skm_err.notice(WARNING,
                       "SAMGp_Solver::SolveMatrixEquation",
                       "Solution criterion was not fulfilled");
      }
    } else {
	    csp_float  signal_to_noise;
	    if (res_in_ == 0.) signal_to_noise = 0.;
	    else               signal_to_noise = res_out_/res_in_;
        
      // make sure that the target residual is strongly violated (best single prec. solution)
      if (signal_to_noise > eps) {
        skm_err.notice(WARNING,
                       "SAMGp_Solver::SolveMatrixEquation",
                       "Solution criterion was not fulfilled");
      }
    }
 }

template class SAMGp_Solver<csp_float,1U>;
template class SAMGp_Solver<csp_float,2U>;
template class SAMGp_Solver<csp_float,3U>;


} // end namespace csp


/*
//#ifndef CSP_PC_WINDOWS_NT_CODE_WARRIOR
void samg_user_coo(int * i,int * ndim, double * x, double * y, double * z) {
    cout << "\nSAMG_Adaptor::SAMG_USER_COO: Setting *ndim to 0 " << endl;
    *ndim = 0;
}
*/

//#endif

/*
// CHECK DATA

    cout<<"\n\nProcessor "<<my_rank<<": Outputting matrix.... "<< endl;
    //Check format:
    //A.Out();
    cout <<"nnu = "<<nnu_<<", nshalo = "<<nshalo<<", nrhalo = "<<nrhalo<< endl;
    cout <<"nna = "<<nna_<< endl;

    cout<<"\nProcessor "<<my_rank<<": Outputting vector x .... "<< endl;
    for (std::vector<csp_float>::const_iterator it=x.begin();it!=x.end();it++) cout << (*it) <<", ";
    cout<<"\nProcessor "<<my_rank<<": Outputting vector b .... "<< endl;
    for (std::vector<csp_float>::const_iterator it=b.begin();it!=b.end();it++) cout << (*it) <<", ";

    // Check samgp_parameters
    cout<<"\n\nisndlist size = "<<samgp_parameters_->nshalo()<<endl;
    for (unsigned int i=0;i!=samgp_parameters_->isndlist.size();i++)
      cout<<samgp_parameters_->isndlist[i]<<", ";
    cout<<"\niranksnd size = "<<samgp_parameters_->npsnd()<<endl;
    for (unsigned int i=0;i!=samgp_parameters_->iranksnd.size();i++)
      cout<<samgp_parameters_->iranksnd[i]<<", ";
    cout<<"\nipts size = "<<samgp_parameters_->ipts.size()<<endl;
    for (unsigned int i=0;i!=samgp_parameters_->ipts.size();i++)
      cout<<samgp_parameters_->ipts[i]<<", ";
    
    cout<<"\n\nireclist size = "<<samgp_parameters_->nrhalo()<<endl;
    for (unsigned int i=0;i!=samgp_parameters_->ireclist.size();i++)
      cout<<samgp_parameters_->ireclist[i]<<", ";
    cout<<"\nirankrec size = "<<samgp_parameters_->nprec()<<endl;
    for (unsigned int i=0;i!=samgp_parameters_->irankrec.size();i++)
      cout<<samgp_parameters_->irankrec[i]<<", ";
    cout<<"\niptr size = "<<samgp_parameters_->iptr.size()<<endl;
    for (unsigned int i=0;i!=samgp_parameters_->iptr.size();i++)
      cout<<samgp_parameters_->iptr[i]<<", ";

    cout<<"\n\ncrmat_.ia (size "<<crmat_.ia.size()<<") = "<<endl;
    for (unsigned int i=0;i!=crmat_.ia.size();i++)
      cout<<crmat_.ia[i]<<", ";
    cout<<"\n\ncrmat_.ja (size "<<crmat_.ja.size()<<") = "<<endl;
    for (unsigned int i=0;i!=crmat_.ja.size();i++)
      cout<<crmat_.ja[i]<<", ";
    cout<<"\n\ncrmat_.a (size "<<crmat_.a.size()<<") = "<<endl;
    for (unsigned int i=0;i!=crmat_.a.size();i++)
      cout<<crmat_.a[i]<<", ";

    cout<<"\n\nProcessor "<<my_rank<<": Calling SAMGp.... "<< endl;
*/

