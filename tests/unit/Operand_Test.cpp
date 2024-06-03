#include "Operand_Test.h"
#include "Operand.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "Model.h"
#include "vsetMakers.h"

using namespace std;

namespace csmp {

void Operand_Test::run()
{
  // 2D test case
  {
    VSet<2U>    mesh_container;
    create_TrianglePatch_VSet( mesh_container );

    // Building Region object from ANSYS data files
    string mesh_name("triangle_patch");
    cout <<"Operand_Test::run: Building Model..."<<endl;
    Model<2U> model( mesh_container, "CSMP-Operand_Test-variables.txt" );

    ScalarVariable      s(PLAIN,1.0);
    VectorVariable<2U>  v(PLAIN,PLAIN,1.0, 1.0);
    TensorVariable<2U>  t(PLAIN,PLAIN,1.0,1.0,1.0,1.0);
    ScalarVariable      ss(PLAIN,2.0);
    VectorVariable<2U>  vv(PLAIN,PLAIN,2.0, 2.0);
    TensorVariable<2U>  tt(PLAIN,PLAIN,2.0,2.0,2.0,2.0);
    model.InputPropertyValue("nodal scalar variable", s );
    model.InputPropertyValue("nodal vector variable", v );
    model.InputPropertyValue("nodal tensor variable", t );

    //Getting the keys;
    csmp::Index sc_key=model.Database().StorageKey("nodal scalar variable");
    csmp::Index vc_key=model.Database().StorageKey("nodal vector variable");
    csmp::Index ts_key=model.Database().StorageKey("nodal tensor variable");

    //Constructing operands
    //from keys only
    Operand<2U> opsc_from_key(sc_key);
    Operand<2U> opvc_from_key(vc_key);
    Operand<2U> opts_from_key(ts_key);

    //from name and database
    Operand<2U> opsc_from_name_database("nodal scalar variable",model.Database());
    Operand<2U> opvc_from_name_database("nodal vector variable",model.Database());
    Operand<2U> opts_from_name_database("nodal tensor variable",model.Database());

    //copy constructed
    Operand<2U> opsc_copy(opsc_from_key);
    Operand<2U> opvc_copy(opvc_from_key);
    Operand<2U> opts_copy(opts_from_key);

    _test(opsc_from_key.Key()==opsc_from_name_database.Key());
    _test(opvc_from_key.Key()==opvc_from_name_database.Key());
    _test(opts_from_key.Key()==opts_from_name_database.Key());

    _test(opsc_from_key.Key()==opsc_copy.Key());
    _test(opvc_from_key.Key()==opvc_copy.Key());
    _test(opts_from_key.Key()==opts_copy.Key());

    _test(opsc_from_key.Key()==sc_key);
    _test(opvc_from_key.Key()==vc_key);
    _test(opts_from_key.Key()==ts_key);

    double val=1.;
    Operand<2U> opsc_from_sc(sc_key,s);
    Operand<2U> opvc_from_vc(vc_key,v);
    Operand<2U> opts_from_ts(ts_key,t);

    // Testing operators (these require operands constructed from storage and keys.. i.e. a scalar, vector, or tensor space plus corresponding key
    _test(opsc_from_sc+val==ss);
    _test(opvc_from_vc+val==vv);
    _test(opts_from_ts+val==tt);

    val=2.;
    _test(opsc_from_sc*val==ss);
    _test(opvc_from_vc*val==vv);
    _test(opts_from_ts*val==tt);

    _test(opsc_from_sc-val==ss*-0.5);
    _test(opvc_from_vc-val==vv*-0.5);
    _test(opts_from_ts-val==tt*-0.5);

    _test(opsc_from_sc/val==ss*0.25);
    _test(opvc_from_vc/val==vv*0.25);
    _test(opts_from_ts/val==tt*0.25);

    //operands with operands
    _test(opsc_from_sc+opsc_from_sc==ss);
    _test(opvc_from_vc+opvc_from_vc==vv);
    _test(opts_from_ts+opts_from_ts==tt);

    _test(opsc_from_sc*opsc_from_sc==ss/2.);
    _test(opvc_from_vc*opvc_from_vc==vv/2.);
    _test(opts_from_ts*opts_from_ts==tt);

    _test(opsc_from_sc-opsc_from_sc==ss*0.);
    _test(opvc_from_vc-opvc_from_vc==vv*0.);
    _test(opts_from_ts-opts_from_ts==tt*0.);

    _test(opsc_from_sc/opsc_from_sc==ss/2.);
    _test(opvc_from_vc/opvc_from_vc==vv/2.);
    _test(opts_from_ts/opts_from_ts==tt/2.);

    //self assignment
    val=1.;
    opsc_from_sc+=val;
    opvc_from_vc+=val;
    opts_from_ts+=val;
    _test(opsc_from_sc==ss);
    _test(opvc_from_vc==vv);
    _test(opts_from_ts==tt);

    val=2;
    opsc_from_sc*=val;
    opvc_from_vc*=val;
    opts_from_ts*=val;
    _test(opsc_from_sc==ss*2.);
    _test(opvc_from_vc==vv*2.);
    _test(opts_from_ts==tt*2.);

    val=1;
    opsc_from_sc-=val;
    opvc_from_vc-=val;
    opts_from_ts-=val;
    _test(opsc_from_sc==ss+1.);
    _test(opvc_from_vc==vv+1.);
    _test(opts_from_ts==tt+1.);

    val=3.;
    opsc_from_sc/=val;
    opvc_from_vc/=val;
    opts_from_ts/=val;
    _test(opsc_from_sc==ss/2.);
    _test(opvc_from_vc==vv/2.);
    _test(opts_from_ts==tt/2.);

    //self assignment with scalar variables.
    ScalarVariable      sss(PLAIN,1.0);
    VectorVariable<2U>  vvv(PLAIN,PLAIN,1.0, 1.0);
    TensorVariable<2U>  ttt(PLAIN,PLAIN,1.0,1.0,1.0,1.0);
    opsc_from_sc+=sss;
    opvc_from_vc+=vvv;
    opts_from_ts+=ttt;
    _test(opsc_from_sc==ss);
    _test(opvc_from_vc==vv);
    _test(opts_from_ts==tt);

    sss=2;
    vvv=2;
    ttt=2;
    opsc_from_sc*=sss;
    opvc_from_vc*=vvv;
    opts_from_ts*=ttt;
    _test(opsc_from_sc==ss*2.);
    _test(opvc_from_vc==vv*2.);
    _test(opts_from_ts==tt*4.);

    sss=1;
    vvv=1;
    ttt=1;
    opsc_from_sc-=sss;
    opvc_from_vc-=vvv;
    opts_from_ts-=ttt;
    _test(opsc_from_sc==ss+1.);
    _test(opvc_from_vc==vv+1.);
    _test(opts_from_ts==tt+5.);

    sss=3;
    vvv=3;
    ttt=7;
    opsc_from_sc/=sss;
    opvc_from_vc/=vvv;
    opts_from_ts/=ttt;
    _test(opsc_from_sc==ss/2.);
    _test(opvc_from_vc==vv/2.);
    _test(opts_from_ts==tt/2.);

    //self assignment with operands
    opsc_from_sc+=opsc_from_sc;
    opvc_from_vc+=opvc_from_vc;
    opts_from_ts+=opts_from_ts;
    _test(opsc_from_sc==ss);
    _test(opvc_from_vc==vv);
    _test(opts_from_ts==tt);

    opsc_from_sc*=opsc_from_sc;
    opvc_from_vc*=opvc_from_vc;
    opts_from_ts*=opts_from_ts;
    _test(opsc_from_sc==ss*2.);
    _test(opvc_from_vc==vv*2.);
    _test(opts_from_ts==tt*4.);

    opsc_from_sc/=opsc_from_sc;
    opvc_from_vc/=opvc_from_vc;
    opts_from_ts/=opts_from_ts;
    _test(opsc_from_sc==ss/2.);
    _test(opvc_from_vc==vv/2.);
    _test(opts_from_ts==tt/2.);

    opsc_from_sc-=opsc_from_sc;
    opvc_from_vc-=opvc_from_vc;
    opts_from_ts-=opts_from_ts;
    _test(opsc_from_sc==ss*0.);
    _test(opvc_from_vc==vv*0.);
    _test(opts_from_ts==tt*0.);

    //logical tests
    val=1.2;
    opsc_from_sc=0.6;
    opvc_from_vc=0.6;
    opts_from_ts=0.6;
    _test(opsc_from_sc<val);

    val=0.4;
    _test(opsc_from_sc>val);
    _test(opsc_from_sc!=val);

    val=0.6;
    _test(opsc_from_sc==val);
    _test(opsc_from_sc<=val);
    _test(opsc_from_sc>=val);

  }

  // 3D test
  {
    VSet<3U>    mesh_container;
    create_Tetra_VSet( mesh_container );

    // Building Region object from ANSYS data files
    string mesh_name("triangle_patch");
    cout <<"Operand_Test::run: Building Model..."<<endl;
    Model<3U> model( mesh_container, "CSMP-Operand_Test-variables.txt" );

    ScalarVariable      s(PLAIN,1.0);
    VectorVariable<3U>  v(PLAIN,PLAIN,PLAIN,1.0, 1.0,1.0);
    TensorVariable<3U>  t(PLAIN,PLAIN,PLAIN,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0);
    ScalarVariable      ss(PLAIN,2.0);
    VectorVariable<3U>  vv(PLAIN,PLAIN,PLAIN,2.0, 2.0,2.0);
    TensorVariable<3U>  tt(PLAIN,PLAIN,PLAIN,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0);

    model.InputPropertyValue("nodal scalar variable", s );
    model.InputPropertyValue("nodal vector variable", v );
    model.InputPropertyValue("nodal tensor variable", t );

    //Getting the keys;
    csmp::Index sc_key=model.Database().StorageKey("nodal scalar variable");
    csmp::Index vc_key=model.Database().StorageKey("nodal vector variable");
    csmp::Index ts_key=model.Database().StorageKey("nodal tensor variable");

    //Constructing operands
    //from keys only
    Operand<3U> opsc_from_key(sc_key);
    Operand<3U> opvc_from_key(vc_key);
    Operand<3U> opts_from_key(ts_key);

    //from name and database
    Operand<3U> opsc_from_name_database("nodal scalar variable",model.Database());
    Operand<3U> opvc_from_name_database("nodal vector variable",model.Database());
    Operand<3U> opts_from_name_database("nodal tensor variable",model.Database());

    //copy constructed
    Operand<3U> opsc_copy(opsc_from_key);
    Operand<3U> opvc_copy(opvc_from_key);
    Operand<3U> opts_copy(opts_from_key);

    _test(opsc_from_key.Key()==opsc_from_name_database.Key());
    _test(opvc_from_key.Key()==opvc_from_name_database.Key());
    _test(opts_from_key.Key()==opts_from_name_database.Key());

    _test(opsc_from_key.Key()==opsc_copy.Key());
    _test(opvc_from_key.Key()==opvc_copy.Key());
    _test(opts_from_key.Key()==opts_copy.Key());

    _test(opsc_from_key.Key()==sc_key);
    _test(opvc_from_key.Key()==vc_key);
    _test(opts_from_key.Key()==ts_key);

    double val=1.;
    Operand<3U> opsc_from_sc(sc_key,s);
    Operand<3U> opvc_from_vc(vc_key,v);
    Operand<3U> opts_from_ts(ts_key,t);

    // Testing operators (these require operands constructed from storage and keys.. i.e. a scalar, vector, or tensor space plus corresponding key
    _test(opsc_from_sc+val==ss);
    _test(opvc_from_vc+val==vv);
    _test(opts_from_ts+val==tt);

    val=2.;
    _test(opsc_from_sc*val==ss);
    _test(opvc_from_vc*val==vv);
    _test(opts_from_ts*val==tt);

    _test(opsc_from_sc-val==ss*-0.5);
    _test(opvc_from_vc-val==vv*-0.5);
    _test(opts_from_ts-val==tt*-0.5);

    _test(opsc_from_sc/val==ss*0.25);
    _test(opvc_from_vc/val==vv*0.25);
    _test(opts_from_ts/val==tt*0.25);

    //operands with operands
    _test(opsc_from_sc+opsc_from_sc==ss);
    _test(opvc_from_vc+opvc_from_vc==vv);
    _test(opts_from_ts+opts_from_ts==tt);

    _test(opsc_from_sc*opsc_from_sc==ss/2.);
    _test(opvc_from_vc*opvc_from_vc==vv/2.);
    _test(opts_from_ts*opts_from_ts==tt+1.);

    _test(opsc_from_sc-opsc_from_sc==ss*0.);
    _test(opvc_from_vc-opvc_from_vc==vv*0.);
    _test(opts_from_ts-opts_from_ts==tt*0.);

    _test(opsc_from_sc/opsc_from_sc==ss/2.);
    _test(opvc_from_vc/opvc_from_vc==vv/2.);
    _test(opts_from_ts/opts_from_ts==tt/2.);

    //self assignment
    val=1.;
    opsc_from_sc+=val;
    opvc_from_vc+=val;
    opts_from_ts+=val;
    _test(opsc_from_sc==ss);
    _test(opvc_from_vc==vv);
    _test(opts_from_ts==tt);

    val=2;
    opsc_from_sc*=val;
    opvc_from_vc*=val;
    opts_from_ts*=val;
    _test(opsc_from_sc==ss*2.);
    _test(opvc_from_vc==vv*2.);
    _test(opts_from_ts==tt*2.);

    val=1;
    opsc_from_sc-=val;
    opvc_from_vc-=val;
    opts_from_ts-=val;
    _test(opsc_from_sc==ss+1.);
    _test(opvc_from_vc==vv+1.);
    _test(opts_from_ts==tt+1.);

    val=3.;
    opsc_from_sc/=val;
    opvc_from_vc/=val;
    opts_from_ts/=val;
    _test(opsc_from_sc==ss/2.);
    _test(opvc_from_vc==vv/2.);
    _test(opts_from_ts==tt/2.);

    //self assignment with scalar variables.
    ScalarVariable      sss(PLAIN,1.0);
    VectorVariable<3U>  vvv(PLAIN,PLAIN,PLAIN,1.0, 1.0,1.0);
    TensorVariable<3U>  ttt(PLAIN,PLAIN,PLAIN,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0);
    opsc_from_sc+=sss;
    opvc_from_vc+=vvv;
    opts_from_ts+=ttt;
    _test(opsc_from_sc==ss);
    _test(opvc_from_vc==vv);
    _test(opts_from_ts==tt);

    sss=2;
    vvv=2;
    ttt=2;
    opsc_from_sc*=sss;
    opvc_from_vc*=vvv;
    opts_from_ts*=ttt;
    _test(opsc_from_sc==ss*2.);
    _test(opvc_from_vc==vv*2.);
    _test(opts_from_ts==tt*6.);

    sss=1;
    vvv=1;
    ttt=1;
    opsc_from_sc-=sss;
    opvc_from_vc-=vvv;
    opts_from_ts-=ttt;
    _test(opsc_from_sc==ss+1.);
    _test(opvc_from_vc==vv+1.);
    _test(opts_from_ts==tt+9.);

    sss=3;
    vvv=3;
    ttt=11;
    opsc_from_sc/=sss;
    opvc_from_vc/=vvv;
    opts_from_ts/=ttt;
    _test(opsc_from_sc==ss/2.);
    _test(opvc_from_vc==vv/2.);
    _test(opts_from_ts==tt/2.);

    //self assignment with operands
    opsc_from_sc+=opsc_from_sc;
    opvc_from_vc+=opvc_from_vc;
    opts_from_ts+=opts_from_ts;
    _test(opsc_from_sc==ss);
    _test(opvc_from_vc==vv);
    _test(opts_from_ts==tt);

    opsc_from_sc*=opsc_from_sc;
    opvc_from_vc*=opvc_from_vc;
    opts_from_ts*=opts_from_ts;
    _test(opsc_from_sc==ss*2.);
    _test(opvc_from_vc==vv*2.);
    _test(opts_from_ts==tt*6.);

    opsc_from_sc/=opsc_from_sc;
    opvc_from_vc/=opvc_from_vc;
    opts_from_ts/=opts_from_ts;
    _test(opsc_from_sc==ss/2.);
    _test(opvc_from_vc==vv/2.);
    _test(opts_from_ts==tt/2.);

    opsc_from_sc-=opsc_from_sc;
    opvc_from_vc-=opvc_from_vc;
    opts_from_ts-=opts_from_ts;
    _test(opsc_from_sc==ss*0.);
    _test(opvc_from_vc==vv*0.);
    _test(opts_from_ts==tt*0.);

    //logical tests
    val=1.2;
    opsc_from_sc=0.6;
    opvc_from_vc=0.6;
    opts_from_ts=0.6;
    _test(opsc_from_sc<val);

    val=0.4;
    _test(opsc_from_sc>val);
    _test(opsc_from_sc!=val);

    val=0.6;
    _test(opsc_from_sc==val);
    _test(opsc_from_sc<=val);
    _test(opsc_from_sc>=val);
  }

} // end run

} // csmp
