#include "ANSYS_Model2D.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_Integrator.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include"NumIntegral_BT_op_dV.h"
#include"NumIntegral_PT_op_dV.h"
#include "StressesAndStrains2.h"
#include "ExtractTensorVariableComponent.h"
#include "BoreHole_stability2D_VVCase.h"
#include "Test.h"

using namespace std;


namespace csmp
  {

BoreHole_stability2D_VVCase::BoreHole_stability2D_VVCase(const char* prefix)

{
    this->setName("BoreHole_stability2D_VVCase");
    prefix_=prefix;
}

 void BoreHole_stability2D_VVCase::run()
    {
      // some constants

      // R=13.3 m
      const double SHmax( 90.0e+06 );
      const double SHmin( 51.5e+06);
      const double Pmud( 31.5e+06 ); //mud pressure
      const double Pp( 31.5e+06 ); //pore pressure
      const double length(20.);

     // pore presure = mud pressure

      const ScalarVariable zeroScalar( PLAIN, 0. );
      const VectorVariable<2U> zeroVector( PLAIN, 0. );
   //   const VectorVariable<DIM> zeroVectorDirichlet( DIRICH, 0. );

      // establishing model & output facility
      string input_file_name(prefix_);

      //  string modelName( "WellBore2D" );


      ANSYS_Model2D model( input_file_name.data(),this->getName().c_str(),(this->getName()+".txt").c_str() );

      VTU_Interface<2U> vtu( model );
      vtu.OmitZeroInFileName(true);

      // output properties & initial output
      list<string> outputProps;
      outputProps.push_back( "mean stress" );
      model.InputPropertyValue( "mean stress", zeroScalar );
      model.InputPropertyValue( "displacement", zeroVector );
      model.InputPropertyValue( "force", zeroVector );
      model.InputPropertyValue( "Neumann stress", zeroVector );
      model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN, Pp) );
      vtu.OutputDataToVTU( ((string)(input_file_name.data())+"UnLoaded").data(), outputProps, "Model", static_cast<int>(0) );

      // model configuration

      model.InputPropertyValue( "Poisson's ratio", makeScalar( PLAIN, 0.29 ) );
      model.InputPropertyValue( "Young's modulus", makeScalar( PLAIN, 24.0e+9 ) );

      Boundary<2U>& rightBoundary( model.Boundary( "RIGHT" ) );
      Boundary<2U>& leftBoundary( model.Boundary( "LEFT" ) );
      Boundary<2U>& topBoundary( model.Boundary( "TOP" ) );
      Boundary<2U>& bottomBoundary( model.Boundary( "BOTTOM" ) );
      Region<2U>& wellbore( model.Region("STANDARD") );

      const VectorVariable<2U> SressRight( DIRICH,DIRICH, -SHmax,0. );
      const VectorVariable<2U> SressLeft( DIRICH,DIRICH, SHmax,0. );
      const VectorVariable<2U> SressTop( DIRICH,DIRICH, 0.,-SHmin );
      const VectorVariable<2U> SressBottom( DIRICH,DIRICH, 0.,SHmin );

      const ScalarVariable Pressure( DIRICH,Pmud);

      rightBoundary.InputPropertyValue( "Neumann stress", SressRight );
      leftBoundary.InputPropertyValue( "Neumann stress", SressLeft );
      topBoundary.InputPropertyValue( "Neumann stress", SressTop );
      bottomBoundary.InputPropertyValue( "Neumann stress", SressBottom );

      wellbore.InputPropertyValue( "fluid pressure", Pressure );

      // setting up & solving linear elasticity fea problem
      #ifdef CSMP_WITH_SAMG_SOLVER
      SAMG_Settings settings;
      settings.Set_napproach(2);
      SAMG_Solver solver(&settings);
      PDE_Integrator<2U,Region> deformation( solver );
      #else
      CSMP_DEFAULT_LINEAR_SOLVER solver;
      PDE_Integrator<2U,Region> deformation( solver );
      #endif
      PT_op<2U,Element<2U> > bforces( model.Database(), "force", "displacement" );
      NumIntegral_BT_D_B_dV<2U> stiffness( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );

      NumIntegral_PT_op_dV<2U>  AppliedStress( model.Database(), "Neumann stress", "displacement");

      NumIntegral_BT_op_dV<2U>  WellBorePressure( model.Database(),"fluid pressure", "displacement");

      deformation.Add( &stiffness );
      deformation.Add( &bforces );
      deformation.Add( &AppliedStress );
      deformation.Add( &WellBorePressure );

      StressesAndStrains<2U>  postpro( model, "Young's modulus", "Poisson's ratio", "displacement", true, false );

      postpro.PlaneStress();

      deformation.AddPostProcess( &postpro );
      model.Apply( deformation );


      ExtractTensorVariableComponent<2U>  ystress(   model.Database(), "stress", "stress-y",  1,1 );
      ExtractTensorVariableComponent<2U>  stress_xy( model.Database(), "stress", "stress-xy", 0,1 );
      ExtractTensorVariableComponent<2U>  xstress( model.Database(), "stress", "stress-x", 0,0 );

      model.Apply( ystress );
      model.Apply( stress_xy );
      model.Apply( xstress );

      // apply resulting displacement
      model.MoveNodeCoordinatesBy("displacement");

      // reset connectivity (displacement!) and output
      outputProps.push_back( "displacement" );
      outputProps.push_back( "strain" );
      outputProps.push_back( "sigma1" );
      outputProps.push_back( "sigma2" );
      outputProps.push_back( "stress-x" );
      outputProps.push_back( "stress-xy" );
      outputProps.push_back( "stress-y" );
      outputProps.push_back( "Neumann stress" );
      vtu.DeleteConnectivity();
      vtu.OutputDataToVTU( ((string)(input_file_name.data())+"Loaded").data(), outputProps, "Model", static_cast<int>(0) );

//-----------------------------------------------------------------
// Tests
//-----------------------------------------------------------------


      ifstream fin;
      string str;

      int NumCompPoint(40);// 41 points were used for comparaison
      double tolerance( 10E-6 );

      int test1(0);
      int test2(0);
      int test3(0);
      int test4(0);
      int test5(0);
      int test6(0);


      //mean stress

      //horizontal direction

      fin.open("BoreHole_stability2D_VVCase_ComparisonData_meanstressH.txt");


      for (size_t compit = 0; compit<NumCompPoint; compit++)

      {

       cout << "test Point " << compit << endl;

       std::vector<double> xy(2);

       getline(fin,str,'\t');
       double xcoord=atof(str.c_str());

       xy[0] = xcoord;

//       cout << "xy[0] " << xy[0] <<endl;

       getline(fin,str,'\t');
       double ycoord=atof(str.c_str());

       xy[1] = ycoord;

//       cout << "xy[1] " << xy[1] <<endl;

       getline(fin,str,'\n');

       double valmeanstress(atof(str.c_str()));

       cout << " valmeanstress " << valmeanstress <<endl;

       // checking the element that contains this point

       Region<2> mod_domain(model.Region("Model"));
       for ( auto it=mod_domain.CellsBegin(); it != mod_domain.CellsEnd(); ++it  )
       {

//        cout << "CurrentID " << (*it)-> Idx() <<endl;

        Point<2U>  point;
        point=(*it)->BaryCenter();

//                cout <<  "point.Coordinates(1) " << point.Coordinates()[0] << "point.Coordinates(1)" << point.Coordinates()[1] <<endl;

                double dis(0.) ;
                dis = (point.Coordinates()[0] - xy[0]) * (point.Coordinates()[0] - xy[0]);
                dis += (point.Coordinates()[1] - xy[1]) * (point.Coordinates()[1] - xy[1]);
                dis = sqrt(dis);
                if (dis < 2.)
                {
       (*it)->CoordinateMatrix();

       (*it)->FE()->N( (*it)->FE()->NRST ,xy);

       double area(0.);
       bool bvar=true;
       for ( auto i{0}; i<(*it)->Nodes(); i++ ) {
           if ((*it)->FE()->NRST[i] < -0.001) bvar=false;
           area += (*it)->FE()->NRST[i];
//                    cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
       }
       if (bvar){

                    cout<<"THIS IS THE RIGHT ELEMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

           double ms(0.);

           for ( auto i{0}; i<(*it)->Nodes(); i++ ){
              ms += (*it)->FE()->NRST[i] * (*it)->N(i)->Read( model.Database().StorageKey("mean stress"));
//              cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
           }
           cout <<" ms " << ms<< endl;
          if (abs((valmeanstress-ms)/valmeanstress)<tolerance)
              test1=test1+1;;

        break;
       }

      }
       }

 }
      fin.close(); // end mean stress horizontal test

       //vertical direction

      fin.open("BoreHole_stability2D_VVCase_ComparisonData_meanstressV.txt");


      for (size_t compit = 0; compit<NumCompPoint; compit++)

      {

       cout << "test Point " << compit << endl;

       std::vector<double> xy(2);

       getline(fin,str,'\t');
       double xcoord=atof(str.c_str());

       xy[0] = xcoord;

//       cout << "xy[0] " << xy[0] <<endl;

       getline(fin,str,'\t');
       double ycoord=atof(str.c_str());

       xy[1] = ycoord;

//       cout << "xy[1] " << xy[1] <<endl;

       getline(fin,str,'\n');

       double valmeanstress(atof(str.c_str()));

       cout << " valmeanstress " << valmeanstress <<endl;

       // checking the element that contains this point
       Region<2U>& model_domain(model.Region("Model"));

       for ( auto it( model_domain.CellsBegin() ); it != model_domain.CellsEnd(); ++it  )
       {

//        cout << "CurrentID " << (*it)-> Idx() <<endl;

        Point<2U>  point;
        point=(*it)->BaryCenter();

//                cout <<  "point.Coordinates(1) " << point.Coordinates()[0] << "point.Coordinates(1)" << point.Coordinates()[1] <<endl;

                double dis(0.) ;
                dis = (point.Coordinates()[0] - xy[0]) * (point.Coordinates()[0] - xy[0]);
                dis += (point.Coordinates()[1] - xy[1]) * (point.Coordinates()[1] - xy[1]);
                dis = sqrt(dis);
                if (dis < 2.)
                {

       (*it)->CoordinateMatrix();

       (*it)->FE()->N( (*it)->FE()->NRST ,xy);

       double area(0.);
       bool bvar=true;
       for ( auto i{0}; i<(*it)->Nodes(); i++ ) {
           if ((*it)->FE()->NRST[i] < -0.001 ) bvar=false;
           area += (*it)->FE()->NRST[i];
//                    cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
       }
       if (bvar){

                    cout<<"THIS IS THE RIGHT ELEMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

           double ms(0.);

           for ( auto i{0}; i<(*it)->Nodes(); i++ ){
              ms += (*it)->FE()->NRST[i] * (*it)->N(i)->Read( model.Database().StorageKey("mean stress"));
//              cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
          }

           cout <<" ms " << ms<< endl;
          if (abs((valmeanstress-ms)/valmeanstress)<tolerance)
              test2=test2+1;;

        break;
       }

      }
       }

 }
      fin.close(); // end mean stress vertical check

      // stress-x

      //horizontal direction

      fin.open("BoreHole_stability2D_VVCase_ComparisonData_stressxH.txt");


      for (size_t compit = 0; compit<NumCompPoint; compit++)

      {

       cout << "test Point " << compit << endl;

       std::vector<double> xy(2);

       getline(fin,str,'\t');
       double xcoord=atof(str.c_str());

       xy[0] = xcoord;

//       cout << "xy[0] " << xy[0] <<endl;

       getline(fin,str,'\t');
       double ycoord=atof(str.c_str());

       xy[1] = ycoord;

//       cout << "xy[1] " << xy[1] <<endl;

       getline(fin,str,'\n');

       double valstressx(atof(str.c_str()));

       cout << " valstressx " << valstressx <<endl;


       // checking the element that contains this point
       Region<2>& mod_domain(model.Region("Model"));
       for ( auto it=mod_domain.CellsBegin(); it != mod_domain.CellsEnd(); ++it  )
       {

//        cout << "CurrentID " << (*it)-> Idx() <<endl;

        Point<2U>  point;
        point=(*it)->BaryCenter();

//                cout <<  "point.Coordinates(1) " << point.Coordinates()[0] << "point.Coordinates(1)" << point.Coordinates()[1] <<endl;

                double dis(0.) ;
                dis = (point.Coordinates()[0] - xy[0]) * (point.Coordinates()[0] - xy[0]);
                dis += (point.Coordinates()[1] - xy[1]) * (point.Coordinates()[1] - xy[1]);
                dis = sqrt(dis);
                if (dis < 2.)
                {

       (*it)->CoordinateMatrix();

       (*it)->FE()->N( (*it)->FE()->NRST ,xy);

       double area(0.);
       bool bvar=true;
       for ( auto i{0}; i<(*it)->Nodes(); i++ ) {
           if ((*it)->FE()->NRST[i] <-0.001 ) bvar=false;
           area += (*it)->FE()->NRST[i];
//                    cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
       }
       if (bvar){

                    cout<<"THIS IS THE RIGHT ELEMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

           double sigmax(0.);

           for ( auto i{0}; i<(*it)->Nodes(); i++ ){
              sigmax += (*it)->FE()->NRST[i] * (*it)->N(i)->Read( model.Database().StorageKey("stress-x"));
//              cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
          }

           cout <<" sigmax " << sigmax<< endl;
          if (abs((valstressx-sigmax)/valstressx)<tolerance)
              test3=test3+1;;

        break;
       }
      }
      }

 }      fin.close(); // end stress-x horizontal test

//       //vertical direction

      fin.open("BoreHole_stability2D_VVCase_ComparisonData_stressxV.txt");


      for (size_t compit = 0; compit<NumCompPoint; compit++)

      {

       cout << "test Point " << compit << endl;

       std::vector<double> xy(2);

       getline(fin,str,'\t');
       double xcoord=atof(str.c_str());

       xy[0] = xcoord;

//       cout << "xy[0]  " << xy[0] <<endl;

       getline(fin,str,'\t');
       double ycoord=atof(str.c_str());

       xy[1] = ycoord;

//       cout << "xy[1]  " << xy[1] <<endl;

       getline(fin,str,'\n');

       double valstressx(atof(str.c_str()));

         cout << " valstressx " << valstressx <<endl;


       // checking the element that contains this point
       Region<2> mod_domain(model.Region("Model"));
       for ( auto it=mod_domain.CellsBegin(); it != mod_domain.CellsEnd(); ++it  )
       {

//        cout << "CurrentID " << (*it)-> Idx() <<endl;

        Point<2U>  point;
        point=(*it)->BaryCenter();

//                cout <<  "point.Coordinates(1) " << point.Coordinates()[0] << "point.Coordinates(1)" << point.Coordinates()[1] <<endl;

                double dis(0.) ;
                dis = (point.Coordinates()[0] - xy[0]) * (point.Coordinates()[0] - xy[0]);
                dis += (point.Coordinates()[1] - xy[1]) * (point.Coordinates()[1] - xy[1]);
                dis = sqrt(dis);
                if (dis < 2.)
                {

       (*it)->CoordinateMatrix();

       (*it)->FE()->N( (*it)->FE()->NRST ,xy);

       double area(0.);
       bool bvar=true;
       for ( auto i{0}; i<(*it)->Nodes(); i++ ) {
           if ((*it)->FE()->NRST[i] <-0.001 ) bvar=false;
           area += (*it)->FE()->NRST[i];
//                    cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
       }
       if (bvar){

                    cout<<"THIS IS THE RIGHT ELEMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

           double sigmax(0.);

           for ( auto i{0}; i<(*it)->Nodes(); i++ ){
              sigmax += (*it)->FE()->NRST[i] * (*it)->N(i)->Read( model.Database().StorageKey("stress-x"));
//              cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
          }

           cout << "sigmax " << sigmax << endl;
          if (abs((valstressx-sigmax)/valstressx)<tolerance)
              test4=test4+1;;

        break;
       }

      }
      }
      }

      fin.close(); // end stress-x vertical check

//      // stress-xy

//      //horizontal direction

      fin.open("BoreHole_stability2D_VVCase_ComparisonData_stressxyH.txt");


      for (size_t compit = 0; compit<NumCompPoint; compit++)

      {

       cout << "test Point " << compit << endl;

       std::vector<double> xy(2);

       getline(fin,str,'\t');
       double xcoord=atof(str.c_str());

       xy[0] = xcoord;

//       cout << "xy[0]  " << xy[0] <<endl;

       getline(fin,str,'\t');
       double ycoord=atof(str.c_str());

       xy[1] = ycoord;

//       cout << "xy[1]  " << xy[1] <<endl;

       getline(fin,str,'\n');

       double valstressxy(atof(str.c_str()));

          cout << " valstressxy " << valstressxy <<endl;


       // checking the element that contains this point
       // TODO: refactor this crazy code
       Region<2> mod_domain(model.Region("Model"));

       for ( auto it=mod_domain.CellsBegin(); it != mod_domain.CellsEnd(); ++it  )
       {

//       cout << "CurrentID " << (*it)-> Idx() <<endl;

       Point<2U>  point;
       point=(*it)->BaryCenter();

//               cout <<  "point.Coordinates(1) " << point.Coordinates()[0] << "point.Coordinates(1)" << point.Coordinates()[1] <<endl;

               double dis(0.) ;
               dis = (point.Coordinates()[0] - xy[0]) * (point.Coordinates()[0] - xy[0]);
               dis += (point.Coordinates()[1] - xy[1]) * (point.Coordinates()[1] - xy[1]);
               dis = sqrt(dis);
               if (dis < 2.)
               {

       (*it)->CoordinateMatrix();

       (*it)->FE()->N( (*it)->FE()->NRST ,xy);

       double area(0.);
       bool bvar=true;
       for ( auto i{0}; i<(*it)->Nodes(); i++ ) {
           if ((*it)->FE()->NRST[i] <-0.001 ) bvar=false;
           area += (*it)->FE()->NRST[i];
//                    cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
       }
       if (bvar){

                    cout<<"THIS IS THE RIGHT ELEMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

           double sigmaxy(0.);

           for ( auto i{0}; i<(*it)->Nodes(); i++ ){
              sigmaxy += (*it)->FE()->NRST[i] * (*it)->N(i)->Read( model.Database().StorageKey("stress-xy"));
//              cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
          }

           cout << "sigmaxy " << sigmaxy << endl;

          if (abs((valstressxy-sigmaxy)/valstressxy)<tolerance)
              test5=test5+1;;

        break;
       }

      }
      }
      }
      fin.close(); // end stress-xy horizontal test

       //vertical direction

      fin.open("BoreHole_stability2D_VVCase_ComparisonData_stressxyV.txt");


      for (size_t compit = 0; compit<NumCompPoint; compit++)

      {

       cout << "test Point " << compit << endl;

       std::vector<double> xy(2);

       getline(fin,str,'\t');
       double xcoord=atof(str.c_str());

       xy[0] = xcoord;

//       cout << "xy[0]  " << xy[0] <<endl;

       getline(fin,str,'\t');
       double ycoord=atof(str.c_str());

       xy[1] = ycoord;

//       cout << "xy[1]  " << xy[1] <<endl;

       getline(fin,str,'\n');

       double valstressxy(atof(str.c_str()));

       cout << " valstressxy " << valstressxy <<endl;


       // checking the element that contains this point
       Region<2> mod_domain(model.Region("Model"));
       for ( auto it=mod_domain.CellsBegin(); it != mod_domain.CellsEnd(); ++it  )
       {
//        cout << "CurrentID " << (*it)-> Idx() <<endl;

        Point<2U>  point;
        point=(*it)->BaryCenter();

//                cout <<  "point.Coordinates(1) " << point.Coordinates()[0] << "point.Coordinates(1)" << point.Coordinates()[1] <<endl;

                double dis(0.) ;
                dis = (point.Coordinates()[0] - xy[0]) * (point.Coordinates()[0] - xy[0]);
                dis += (point.Coordinates()[1] - xy[1]) * (point.Coordinates()[1] - xy[1]);
                dis = sqrt(dis);
                if (dis < 2.)
                {

       (*it)->CoordinateMatrix();

       (*it)->FE()->N( (*it)->FE()->NRST ,xy);

       double area(0.);
       bool bvar=true;
       for ( auto i{0}; i<(*it)->Nodes(); i++ ) {
           if ((*it)->FE()->NRST[i] <-0.001 ) bvar=false;
           area += (*it)->FE()->NRST[i];
//                    cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
       }
       if (bvar){

                    cout<<"THIS IS THE RIGHT ELEMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

           double sigmaxy(0.);

           for ( auto i{0}; i<(*it)->Nodes(); i++ ){
              sigmaxy += (*it)->FE()->NRST[i] * (*it)->N(i)->Read( model.Database().StorageKey("stress-xy"));
//              cout << " (*it)->FE()->NRST[i] " << (*it)->FE()->NRST[i] <<endl;
          }

          cout << "sigmaxy " << sigmaxy << endl;

          if (abs((valstressxy-sigmaxy)/valstressxy)<tolerance)
              test6=test6+1;;

        break;
       }

      }
       }
      }
      fin.close(); // end stress-xy vertical check

      if (test1 == NumCompPoint and test2 == NumCompPoint and test3 == NumCompPoint and test4 == NumCompPoint and test5 == NumCompPoint and test6 == NumCompPoint)
      {
          cout << test1  << " Values verified for mean stress Horizontal "<<endl;
          cout << test2  << " Values verified for mean stress vertical" << endl;
          cout << test3  << " Values verified for stress-x Horizontal"  << endl;
          cout << test4  << " Values verified for stress-x vertical" << endl;
          cout << test5  << " Values verified for stress-xy Horizontal" << endl;
          cout << test6  << " Values verified for stress-xy vertical "<< endl;

          _equal( 0, 0., tolerance );
      }
      else
      {
          cout << "test failure: " <<" just "<< test1  << " values verified for first test  "<<endl;
          cout << "test failure: " <<" just "<< test2  << " values verified for second test  "<<endl;
          cout << "test failure: " <<" just "<< test3  << " values verified for third test  "<<endl;
          cout << "test failure: " <<" just "<< test4  << " values verified for fourth test  "<<endl;
          cout << "test failure: " <<" just "<< test5  << " values verified for fifth test  "<<endl;
          cout << "test failure: " <<" just "<< test6  << " values verified for fifth test  "<<endl;
          _equal( -1, 0., tolerance );
      }

      return;
    }

  } // csmp
