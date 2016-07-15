#include "TotalDensity_Test.h"

using namespace std;

namespace csmp {

TotalDensity_Test::TotalDensity_Test( const string& mesh_name, 
                                      const string& var_file )
{}
  //: tol_(0.001)
/*{
    mi_ = new TRIANGLE_Interface();
   // vset_ = new VSet<csmp_float, 2>();
    
    mi_->ReadTriangle2DMesh( mesh_name.c_str(), *vset_ );
    
  //  sg_ = new SuperGroup<csp_float, 2>(*vset_, var_file.c_str() );
}
*/
TotalDensity_Test::~TotalDensity_Test() {
  /*
    delete mi_;
    delete vset_;
    delete sg_;
    */
}

void TotalDensity_Test::run() {
   // elementTest();
    //nodeTest();
}
/*
void TotalDensity_Test::elementTest() {
 // TotalDensity<csp_float, 2> td(sg_->ReferencePropertyDatabase(),
  //                              "conductivity1",
  //                              "mobility water",
   //                             "mobility oil",
  //                              "pressure");
                                
  vector<csp_float> mw;
  mw.push_back(1.0);
  mw.push_back(2.0);
  mw.push_back(3.0);
  mw.push_back(4.0);
  mw.push_back(5.0);
  setNodeVariable(mw, "mobility water");
  
  vector<csp_float> mo;
  mo.push_back(1.0);
  mo.push_back(1.0);
  mo.push_back(1.0);
  mo.push_back(1.0);
  mo.push_back(1.0);
  setNodeVariable(mo, "mobility oil");
  
  vector<csp_float> p;
  p.push_back(0.0);
  p.push_back(0.0);
  p.push_back(1.0);
  p.push_back(1.0);
  p.push_back(0.5);
  setNodeVariable(p, "pressure");
  
  sg_->Pass(td);
  
  // Set expected values
    vector<csp_float> expected;
    expected.push_back(21.0/6.0);
    expected.push_back(23.0/18.0);
    expected.push_back(13.0/6.0);
    expected.push_back(13.0/6.0);
  
  ScalarVariable<csp_float> sc;
    std::deque<Element<csp_float, 2> >::iterator it;
    csp::Index key(sg_->ReferencePropertyDatabase().StorageKey("conductivity1"));
    
    unsigned int i = 0;
    for (it = sg_->ReferenceMesh().ElementsBegin(); it != sg_->ReferenceMesh().ElementsEnd(); ++it) {
        it->Read(sg_->ReferencePropertyStorage(), key.index, sc);
        _equal(sc(), expected.at(i), tol_);
        ++i;
    }
}

void TotalDensity_Test::nodeTest() {
      TotalDensity<csp_float, 2> td(sg_->ReferencePropertyDatabase(),
                                "mobility",
                                "mobility water",
                                "mobility oil",
                                "pressure");
                                
  vector<csp_float> mw;
  mw.push_back(1.0);
  mw.push_back(2.0);
  mw.push_back(3.0);
  mw.push_back(4.0);
  mw.push_back(5.0);
  setNodeVariable(mw, "mobility water");
  
  vector<csp_float> mo;
  mo.push_back(1.0);
  mo.push_back(1.0);
  mo.push_back(1.0);
  mo.push_back(1.0);
  mo.push_back(1.0);
  setNodeVariable(mo, "mobility oil");
  
  vector<csp_float> p;
  p.push_back(0.0);
  p.push_back(0.0);
  p.push_back(1.0);
  p.push_back(1.0);
  p.push_back(0.5);
  setNodeVariable(p, "pressure");
  
  sg_->Pass(td);
  
  vector<csp_float> expected;
  expected.push_back(1.0);
  expected.push_back(1.0);
  expected.push_back(3.0);
  expected.push_back(4.0);
  expected.push_back(3.0);
  
    ScalarVariable<csp_float> sc;
    std::deque<Node<csp_float, 2> >::iterator it;
    csp::Index key(sg_->ReferencePropertyDatabase().StorageKey("mobility"));
    
    unsigned int i = 0;
    for (it = sg_->ReferenceMesh().NodesBegin(); it != sg_->ReferenceMesh().NodesEnd(); ++it) {
        it->Read(sg_->ReferencePropertyStorage(), key.index, sc);
        _equal(sc(), expected.at(i), tol_);
        ++i;
    }

}

void TotalDensity_Test::setNodeVariable(vector<csp_float>& var, const char* var_name ) {
    std::deque<Node<csp_float, 2> >::iterator it;
    csp::Index key(sg_->ReferencePropertyDatabase().StorageKey(var_name));
        
    unsigned int i = 0;
    for (it = sg_->ReferenceMesh().NodesBegin(); it != sg_->ReferenceMesh().NodesEnd(); ++it) {     
        it->Store(sg_->ReferencePropertyStorage(),
                  key.index,
                  ScalarVariable<csp_float>(PLAIN, static_cast<csp_float>(var[i])));
        
        ++i;
    }
}
*/
} // csmp
