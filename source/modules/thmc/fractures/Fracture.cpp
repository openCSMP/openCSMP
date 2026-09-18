// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#define _USE_MATH_DEFINES
#include <cmath>
#include "Fracture.h"
#include "Exception.h"
#include "Model.h"
#include <regex>

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif


//need to provide and link GSL library if you want to use this functionality
#define gsl_sf_hyperg_2F1(X,Y,Z,C) NAN;


namespace csmp {

template<uint32_t dim>
Fracture<dim>::Fracture(Model<dim>& model, std::string splitboundary, TIP_TYPE tip_type ):
  model_(&model),
  sb_ref_(model.SplitBoundary(splitboundary)),
  tiptype_(tip_type),
  configured_(false)
{

  if ( (*(sb_ref_.CellsBegin()))->HasInterveningElement() ){
        std::string reg_name = std::regex_replace(sb_ref_.Name(), std::regex("SPLITBOUNDARY"), "REGION");
        std::cout << "Fracture class looking for middle region: " << reg_name << std::endl;
        midregion_ = &model.Region(reg_name);
    }

  //0) defining default keys
  disp_key_ = model.Database().StorageKey("displacement");
  nu_key_   = model.Database().StorageKey("Poisson's ratio");
  E_key_    = model.Database().StorageKey("Young's modulus");

  //0.5) Defining propagation regime limits
  K_dominant = 5.0;                       //where toughness dominated regime prevails
  M_small    = 0.1;                       //Up to where small effects of viscosity should be included
  K_small    = 1.0;                       //Up to where toughness should be included
  M_dominant = std::pow(0.5, -4.0);       //Where viscous dominated regime prevails


  //1) Configuring Initial Length
  //SetOldFractureLengthToCurrent();  EP TODO BRING BACK FRACTURE TIP   //needed for old fracture tip construction (not based on dx)

  //2)Creating Fracture tips container
  //CreateFractureTips(); TODO FRACTURE TIP FT

} // end of construction of Fracture



/**
 Destructor for Fracture tips that are dynamically created
 */
template<uint32_t dim>
Fracture<dim>::~Fracture(){

  /* TODO FRACTURE TIP FT
  //Deleting memory allocated in fracture tip
  for (typename std::vector<FractureTip<dim>*>::iterator it = fracturetips_.begin(); it != fracturetips_.end(); ++it){
    delete *it;
  }
  std::cout << "Fracture::Destructor() -- Deleting dynamic memory allocated to fracture tips" << std::endl;

  */

}










/**
  Initialises fracture tips :
    1) Finds tip nods
    2) Find boundary nodes
    3) Get fracture length
    4) Create FractureTip for each tip node
    5) Set Current tip coordinate as old tip coord

    Note: FractureTips are created thinking there is not a situation where they are created and have already propagated! Therefore they do not need old tip coord!

 */
/* TODO FRACTURE TIP FT
template<uint32_t dim>
void Fracture<dim>::CreateFractureTips(){


  //Identifying Tip Nodes in split boundary
  std::vector<Node<dim>*> tipnodes;
  boxboundaryNodes_.clear();
  std::cout << "Fracture::CreateFractureTips()... " << std::endl;
  for ( typename std::vector< InterFace<dim>* >::const_iterator it = sb_ref_.PerimeterCellsBegin();
        it != sb_ref_.CellsEnd(); it++){
      size_t nds_on_face = (*it)->FE()->Nodes();
      for (size_t i = 0; i < nds_on_face; i++ ){
        if ( (*it)->N(i,INSIDE) == (*it)->N(i,OUTSIDE)){
          tipnodes.push_back((*it)->N(i,INSIDE));
          std::cout << "\t...adding tip node " << (*it)->N(i,INSIDE)->Idx() << " from InterFace " << (*it)->Idx() << std::endl;
        }
        else if ( (*it)->N(i,INSIDE)->AtBoundary()  != INTERNAL ){    //only internal or box boundaries should be possible
          assert( (*it)->N(i,INSIDE)->AtBoundary()  != NOT );    //check it really is a box boundary
          assert( (*it)->N(i,OUTSIDE)->AtBoundary() != INTERNAL and (*it)->N(i,OUTSIDE)->AtBoundary() != NOT ); //check outside element the same!
          boxboundaryNodes_.push_back( (*it)->N(i,INSIDE));
          boxboundaryNodes_.push_back( (*it)->N(i,OUTSIDE));
        }
      }
  }

  //checking tip numbers detected (2d only)
  if (tipnodes.size() == 2 ){
    two_tips_ = true;
    assert(boxboundaryNodes_.size() == 0);
  } else {
    two_tips_ = false;
    if (dim == 2)
      assert(boxboundaryNodes_.size() == 2);
  }

  std::cout << "\tTotal tip nodes found: " << tipnodes.size()
            << "\n\tTotal boundary nodes found: " << boxboundaryNodes_.size() << std::endl;


  //getting half length
  double l_current = FractureLength(MIDDLE);
  if (two_tips_)
    l_current *= 0.5;     //get half length


  switch (tiptype_) {
  case COHESIVE_TIP:
    for ( size_t i = 0; i < tipnodes.size(); i++){
      InterFace<dim>* if_obj = FindInterFaceOfNode(tipnodes[i], INSIDE, PERIMETER).begin()->first;
      fracturetips_.push_back( new CohesiveTip<dim>( model_->Database(), if_obj, l_current ));
    }
    break;
  case DC_TIP:
    for ( size_t i = 0; i < tipnodes.size(); i++){
      InterFace<dim>* tip_interface = FindInterFaceOfNode(tipnodes[i], INSIDE, PERIMETER).begin()->first;
      fracturetips_.push_back( new DisplacementCorrelationTip<dim>( model_->Database(), tip_interface, l_current));
    }
    break;
  case JA_HFM_TIP:
    for ( size_t i = 0; i < tipnodes.size(); i++){
      InterFace<dim>* tip_interface = FindInterFaceOfNode(tipnodes[i], INSIDE, PERIMETER).begin()->first;
      fracturetips_.push_back( new J_IntegralViscousTip<dim>( model_->Database(), tip_interface, JA_HFM_TIP, l_current ));
    }
    break;
  case J_HFM_TIP:
    for ( size_t i = 0; i < tipnodes.size(); i++){
      InterFace<dim>* tip_interface = FindInterFaceOfNode(tipnodes[i], INSIDE, PERIMETER).begin()->first;
      fracturetips_.push_back( new J_IntegralViscousTip<dim>( model_->Database(), tip_interface, J_HFM_TIP, l_current, true ));
    }
    break;
  case J_WET_TIP:
    for ( size_t i = 0; i < tipnodes.size(); i++){
      InterFace<dim>* tip_interface = FindInterFaceOfNode(tipnodes[i], INSIDE, PERIMETER).begin()->first;
      fracturetips_.push_back( new J_IntegralTip<dim>( model_->Database(), tip_interface, J_WET_TIP, l_current ));
    }
    break;
  case J_DRY_TIP:
    for ( size_t i = 0; i < tipnodes.size(); i++){
      InterFace<dim>* tip_interface = FindInterFaceOfNode(tipnodes[i], INSIDE, PERIMETER).begin()->first;
      fracturetips_.push_back( new J_IntegralTip<dim>( model_->Database(), tip_interface, J_DRY_TIP, l_current));
    }
    break;
  case HF_TIP:
    for (size_t i = 0; i < tipnodes.size(); i++){
      InterFace<dim>* tip_interface = FindInterFaceOfNode(tipnodes[i], INSIDE, PERIMETER).begin()->first;
      fracturetips_.push_back( new HydraulicFractureTip<dim>( model_->Database(), tip_interface, l_current));
      std::cout << "\t...FractureTip created with tip " << fracturetips_[i]->TipNode()->Idx() << std::endl;
    }
    break;
  } // end of switch


} //end of Create fracture tips

*/







/**
 Method to call when advancing time steps and need to tell Fracture tips their current tip is now the tip at old time step
 */
/* TODO FRACTURE TIP FT
template<uint32_t dim>
void Fracture<dim>::SetOldTipCoordinatesAsCurrent(){

  for (FractureTip<dim>* F : fracturetips_){
    F->SetOldTipCoordinateAsCurrent();
  }

}
*/


/*

// Making tip elements quarter-point:
//This may incur a rotation of the nodes within the tip element to have the tip node as E->N(0).
//attention: this will not work when one element has two tips!
template<uint32_t dim>
void Fracture<dim>::InitializeQuarterPointElements(){

  if (fracturetips_.empty() ){
    throw csmp::Exception(ERROR, "Fracture<dim>::InitializeQuarterPointElements",
                          "Fracture has no tip elements!");
  }

  std::map<Element<dim>*, size_t > rotated_elements;

  for (typename std::vector<FractureTip<dim>*>::iterator it = fracturetips_.begin(); it != fracturetips_.end(); ++it){
    Node<dim>* tipNode = (*it)->TipNode();
    (*it)->InitializeFirstSet();
    std::set<Element<dim>*> tipElements = (*it)->TipElements();

    for (typename std::set<Element<dim>*>::const_iterator eit = tipElements.begin(); eit != tipElements.end(); ++eit ){
      if ( (*eit)->Interpolation() != 2)
        throw csmp::Exception(ERROR,"Fracture::InitializeQuarterPointElements()",
                              "Tip Elements must be quadratic order for quarter points!" );
      //Pre-emptively assigning Quarterpoint element
      assert( (*eit)->IsQuarterPoint() == false );
      (*eit)->IsQuarterPoint(true) ;

      //getting old node vector
      std::vector<csmp::Node<dim>*> oldNodeConnector = (*eit)->NodeVector();

      //search for element
      if ( tipNode == (*eit)->N(0) ){
        rotated_elements.insert(std::make_pair(*eit,0) );
      } else if (tipNode == (*eit)->N(1) ){
        (*eit)->Assign(0,oldNodeConnector[1]);
        (*eit)->Assign(1,oldNodeConnector[2]);
        (*eit)->Assign(2,oldNodeConnector[0]);
        (*eit)->Assign(3,oldNodeConnector[4]);
        (*eit)->Assign(4,oldNodeConnector[5]);
        (*eit)->Assign(5,oldNodeConnector[3]);
        rotated_elements.insert(std::make_pair(*eit,1) );// One shift to the counter clockwise

      } else if ( tipNode == (*eit)->N(2) ){
        (*eit)->Assign(0,oldNodeConnector[2]);
        (*eit)->Assign(1,oldNodeConnector[0]);
        (*eit)->Assign(2,oldNodeConnector[1]);
        (*eit)->Assign(3,oldNodeConnector[5]);
        (*eit)->Assign(4,oldNodeConnector[3]);
        (*eit)->Assign(5,oldNodeConnector[4]);
        rotated_elements.insert(std::make_pair(*eit,2) ); // two shift to the counterclockwise

      } else
        throw csmp::Exception(ERROR, "Fracture::InitializeQuarterPointElements",
                              "Tip node not found in tip elements!");

      // Move midside nodes to the quarter point position:
      std::vector<Point<dim>> nodePts;
      for ( typename std::vector<Node<dim>*>::iterator nit = (*eit)->NodesBegin(); nit != (*eit)->NodesEnd(); ++nit)
          nodePts.push_back((*nit)->Coordinate());
      Node<dim>* node3 = (*eit)->N(3);
      Node<dim>* node5 = (*eit)->N(5);
      node3->Coordinate(nodePts[0] + (nodePts[1]-nodePts[0])/4.);
      node5->Coordinate(nodePts[0] + (nodePts[2]-nodePts[0])/4.);
      std::cout << "\tMoved nodes " << node3->Idx() << ", " << node5->Idx() << " of Element " << (*eit)->Idx() << " to quarter point" << std::endl;

    }//end of element loop
  } //end of fracture tip loop

  //check we found some elements
  assert(!rotated_elements.empty());

  //Now we rotate the neighbors elmt connector
  for (typename std::map<Element<dim>*,size_t>::iterator it = rotated_elements.begin(); it != rotated_elements.end(); ++it){
    std::vector<Element<dim>*> vc = it->first->NeighborCellVector();
    std::rotate(vc.begin(), vc.begin() + it->second , vc.end() );        //+ vc.size() - it->secondrotate container to the right
    for (size_t i = 0;  i < vc.size(); ++i){
      it->first->Assign(i,vc[i]);
    }
    //for each node change parent element node number
    for (size_t i_n = 0; i_n < it->first->Nodes(); ++i_n){
      it->first->N(i_n)->Unassign(it->first);
      it->first->N(i_n)->ResizeParentStorage( it->first->N(i_n)->Parents() + 1) ;
      it->first->N(i_n)->Assign( i_n, it->first);
    }
  } //end of neighbour connectivity update

  //Updating tip InterFace elements with new faces corresponding to interface side(nodes are the same, just node numbering changed).
  size_t adjusted_interface_sides = 0;
  for (typename std::vector<InterFace<dim>*>::iterator ifit = sb_ref_.PerimeterElementsBegin();
       ifit != sb_ref_.PerimeterElementsEnd(); ++ifit){

    //Search for tip element on INSIDE
    typename std::map<Element<dim>*,size_t>::iterator found_pair = rotated_elements.find((*ifit)->InnerParent());
    if ( found_pair != rotated_elements.end() ){
      //changing face id by rotation
      assert( (*ifit)->InnerParent()->Faces() == 3 );       // code only made for triangular elements
      size_t new_face_id = ( (*ifit)->InnerParentFaceID() + (*ifit)->InnerParent()->Faces() - found_pair->second ) % (*ifit)->InnerParent()->Faces() ;      // Only works for triangles (3 faces)
      (*ifit)->Assign( found_pair->first, new_face_id, INSIDE );
      ++adjusted_interface_sides;

    }

    //searching for tip Element on OUTSIDE
    found_pair = rotated_elements.find( (*ifit)->OuterParent() );
    if (found_pair != rotated_elements.end() ){
      assert( (*ifit)->OuterParent()->Faces() == 3 );       // code only made for triangular elements
      //changing face id by rotating it
      size_t new_face_id = ( (*ifit)->OuterParentFaceID() + (*ifit)->OuterParent()->Faces() - found_pair->second ) % (*ifit)->OuterParent()->Faces() ;     // only works for triangles (3 face)
      (*ifit)->Assign( found_pair->first, new_face_id, OUTSIDE );
      ++adjusted_interface_sides;
    }
  }// end of interface loop

  if (adjusted_interface_sides != 2*fracturetips_.size() )
    throw csmp::Exception(ERROR, "Fracture::InitializeQuarterPointElements()",
                          "Did not adjust the correct amount of interface objects faces");

  //NOTE: Nodes have stayed the same so nothing else should be updated within interface object
  std::cout << "\nFracture::InitializeQuarterPointsElements -> Nodes Coordinates moved to the Quarter points of tip" << std::endl;

}




template<uint32_t dim>
void Fracture<dim>::RestoreQuarterPointElementsToMidPoint(){

  assert(!fracturetips_.empty());

  for (typename std::vector<FractureTip<dim>*>::iterator it = fracturetips_.begin(); it != fracturetips_.end(); ++it){
    (*it)->InitializeFirstSet();
    std::set<Element<dim>*> tipElements = (*it)->TipElements();
    assert(!tipElements.empty());

    for (typename std::set<Element<dim>*>::const_iterator eit = tipElements.begin(); eit != tipElements.end(); ++eit ){
      //Changing FE type
      assert( (*eit)->IsQuarterPoint() == true );
      (*eit)->IsQuarterPoint(false);

      //changing node back to mid point
      assert( (*it)->TipNode()->Coordinate() == (*eit)->N(0)->Coordinate()); //check elm indeed at tip, otherwise initialization is wrong!
      (*eit)->N(3)->Coordinate( 0.5* ( (*eit)->N(0)->Coordinate() + (*eit)->N(1)->Coordinate() ) );
      (*eit)->N(5)->Coordinate( 0.5* ( (*eit)->N(0)->Coordinate() + (*eit)->N(2)->Coordinate() ) );

    } //end of tip element

  } //end of fracture tips

  std::cout << "Fracture -> Quarter points restored to mid point" << std::endl;

} // end of RestoreQuarterPointElementsToMidPoint()



*/


template<uint32_t dim>
std::map<Point<dim>, Node<dim>*> Fracture<dim>::NodeMap(INTERFACE_SIDE side){
  std::map<Point<dim>, Node<dim>*> nodemap;
  for ( typename std::vector<InterFace<dim>*>::const_iterator ifit = sb_ref_.CellsBegin();
        ifit != sb_ref_.CellsEnd(); ++ifit){
    for (uint32_t n= 0u; n < (*ifit)->FE()->Nodes(); ++n)
      nodemap[(*ifit)->N(n,side)->Coordinate()] = (*ifit)->N(n,side) ;
  }

  return nodemap;
}


//orders nodes by x (0), y (1) or z (2) coordinate
template<uint32_t dim>
std::map<double, Node<dim>*> Fracture<dim>::NodeMap(INTERFACE_SIDE side, size_t xyz){
  std::map<double, Node<dim>*> nodemap;
  for ( auto ifit = sb_ref_.CellsBegin();
        ifit != sb_ref_.CellsEnd(); ++ifit){
    for ( uint32_t n=0u; n < (*ifit)->FE()->Nodes(); ++n)
      nodemap[(*ifit)->N(n,side)->Coordinate()[xyz]] = (*ifit)->N(n,side) ;
  }

  return nodemap;
}



template<uint32_t dim>
std::map<Point<dim>, Element<dim>*> Fracture<dim>::ElementMap(INTERFACE_SIDE side){
  std::map<Point<dim>, Element<dim>*> el_map;
  for ( typename std::vector<InterFace<dim>*>::const_iterator ifit = sb_ref_.CellsBegin();
        ifit != sb_ref_.CellsEnd(); ++ifit){
      el_map[(*ifit)->Parent(side)->BaryCenter() ] = (*ifit)->Parent(side) ;
  }

  return el_map;
}


/*
template<uint32_t dim>
std::vector<Node<dim>*> Fracture<dim>::TipNodes(){
  std::vector<Node<dim>*> node_vec;
  for (size_t i = 0; i < fracturetips_.size(); ++i){
    node_vec.push_back(fracturetips_[i]->TipNode());
  }
  return node_vec;
}
*/




//Attention! Use Subdomain flag with caution since not all interfaces that have a tip node will be on the perimeter!
template<uint32_t dim>
std::map<InterFace<dim>*, size_t> Fracture<dim>::FindInterFaceOfNode(Node<dim>* n_ptr, INTERFACE_SIDE side, SUBDOMAIN_PART part  ){

  //setting scope of iterator loop
  typename std::vector<InterFace<dim>*>::const_iterator  f_begin = sb_ref_.CellsBegin();
  typename std::vector<InterFace<dim>*>::const_iterator  f_end   = sb_ref_.CellsEnd();
  switch (part){
  case COMPLETE: break;
  case PERIMETER: f_begin = sb_ref_.PerimeterCellsBegin();
    break;
  case INTERIOR:  f_end   = sb_ref_.PerimeterCellsBegin();
    break;
  }

  std::map<InterFace<dim>*,size_t> InterFace_and_local_node_found;
  for (typename std::vector<InterFace<dim>*>::const_iterator f = f_begin; f != f_end; ++f){
    for ( uint32_t i = 0u; i < (*f)->FE()->Nodes(); ++i){
      if ( (*f)->N(i,side) == n_ptr ){
        InterFace_and_local_node_found.insert( std::make_pair( *f, i) );
      }
    }//end of for loop
  }

  //check we found a node
  if( InterFace_and_local_node_found.size() == 0)
    throw csmp::Exception(ERROR, "Fracture<dim>::FindInterFaceOfNode()",
                          "No InterFace was found which contains for the input Node!" );

  return InterFace_and_local_node_found;

}// end of Find InterFaceOfNode()





/// Calculates length by calling element area
template<uint32_t dim>
double Fracture<dim>::FractureLength(INTERFACE_SIDE side){

  assert(dim == 2);
  double Length = 0.0;
  if (dim == 2){
    for (auto& ifobj : sb_ref_.CellVector())
      Length += ifobj->Area(side);
  } else if (dim == 3){
    throw csmp::Exception (ERROR, "Fracture::FractureLength()", "3D case not implemented yet");
  }

  return Length;
} // end of FractureLength


template<uint32_t dim>
void Fracture<dim>::SetOldFractureLengthToCurrent(double dt){

  double current_length = FractureLength(MIDDLE);
  if (old_fracture_length_ == current_length)
    accumulated_time_ += dt;  //accumulating time spent at current fracture length
  else if (old_fracture_length_ < current_length)
    accumulated_time_ = 0;    //resetting time spent at current fracture length

  //setting the current fracture length as the old one (done before moving to the next time step)
  old_fracture_length_ = current_length;

}




//Iterates over all tips and return the evaluated value measured at the tip
/*
template<uint32_t dim>
std::vector<double> Fracture<dim>::EvaluateTips(){

  if (fracturetips_.size() == 0)
    throw csmp::Exception(ERROR, "Fracture::EvaluateTips()",
                          "No fracture tips are defined in fracture");

  std::vector<double> tip_prop_values;
  bool found_one_tip_to_propgate = false;
  for (size_t i = 0; i < fracturetips_.size(); ++i){
    //0) Check tips are well defined
    assert( fracturetips_[i] != nullptr);

    //1) assesing criteria at tips
    double dt;
    std::pair<bool,bool> continue_and_split_mesh = fracturetips_[i]->EvaluatePropagationCriterion( dt ) ;
    if ( continue_and_split_mesh.second == true ){
      found_one_tip_to_propgate = true;
    }
    //verbose output
    tip_prop_values.push_back(fracturetips_[i]->OutputEvaluation());
  } //end of tip loop

  return tip_prop_values;
}



///General algorithm evaluating tip propagation criteria, propagation lengths, and enacting mesh splitting
template<uint32_t dim>
bool Fracture<dim>::PropagationAlgorithm(Boundary<dim>& b_ref, double& dt, bool propagate_anyway ){

   // 0) Prelimininary checks
   if( b_ref.Cells() == 0)
      throw csmp::Exception(WARNING, "Fracture::PropagateTips(Boundary)",
                            "No face elements for fracture to split!");


   if ( fracturetips_.empty()  )
     throw csmp::Exception(ERROR, "Fracture<dim>::PropagationAlgorithm()",
                           "No tip nodes or fracture tips exist to propagate!");

   size_t                propagated_tips(0);
   std::vector<double> dt_tips (fracturetips_.size(), dt);
   std::deque<bool>      continue_it(fracturetips_.size(), false);

   for ( size_t i = 0; i < fracturetips_.size(); ++i){

     // 0.5) Check tips are configured correctly
     assert( fracturetips_[i] != nullptr);                                //Check fracture tips alsp initialized

     //1) Starting Algorithm:
     //Evaluate propagation criteria of fracture tip
     std::pair<bool,size_t> continue_iteration_and_propagate_mesh  = std::make_pair<bool,size_t>(true,1);
     if (propagate_anyway == false){
     continue_iteration_and_propagate_mesh = fracturetips_[i]->EvaluatePropagationCriterion( dt_tips[i] );
     }
     continue_it[i]      = continue_iteration_and_propagate_mesh.first;
     size_t prop_times   = continue_iteration_and_propagate_mesh.second;
     assert( prop_times >= 0 and prop_times < 20);

     //Propagate mesh based on evaluation criteria result
     if ( prop_times == 0 and propagate_anyway == false){
       std::cout << "\nTip " << i << " is stable" << std::endl;
     } else {
       //check iteration should also continue when we propagate
       assert( continue_iteration_and_propagate_mesh.first == true );


       //Multiple propagation of tip
       for (size_t prop = 0; prop < prop_times; prop++ ){
         std::cout << "\n------------------------------------------------------"
                      "\nPropagation is enacted " << prop_times << " times for Tip " << fracturetips_[i]->TipNode()->Idx()
                   << "\n-----------------------------------------------------" << std::endl;
         //Finding nodes to split on Boundary, and splitting mesh
         propagated_tips++;
         PropagateTip(fracturetips_[i], b_ref);
       }

     } //end of if evaluation
   } //end of fracture tip loop

   std::cout << "\nFracture: Tips propagated " << propagated_tips << std::endl;

   //Minimun time step is chosen as the governing time step
   std::vector<double>::iterator vit = std::min_element(dt_tips.begin(), dt_tips.end());
   size_t controlling_tip = std::distance( dt_tips.begin(), vit);
   dt = *vit;

   // Just fix one tip
   //controlling_tip = 1;
   //dt = dt_tips[controlling_tip];


   assert( dt > 0.0 );
   assert(continue_it.size() == fracturetips_.size());


   //if split occured
   if (propagated_tips != 0 )
     return true;                 //Continue fracture length iteration
   else if ( std::find( continue_it.begin(), continue_it.end(), true ) == continue_it.end() )
     return false;                //All tips have converged and stop iteration (unlikely with multiple tips)
   else if ( continue_it[ controlling_tip ] == true ) {
     std::cout << "\nFracture: Propagation Algorithm\n\tControlling tip " << controlling_tip
               << " Not converged -> Continuing iteration " << std::endl;
     return true;
   } else if ( continue_it[ controlling_tip ] == false ){
     std::cout << "\nFracture: Propagation Algorithm\n\tControlling tip " << controlling_tip
               << " Converged -> Moving to next time step " << std::endl;
     return false;               //Length converged when tip with smallest time step has converged
   } else
     throw csmp::Exception(ERROR, "Fracture::PropagationAlgorithm",
                           "Unforseen propagation result");


} //end of propagation algorithm
*/



/** @author E.Pezzulli
    @date 2020
    @brief General method to split tip nodes along a predefined path specified by
    the boundary object.

*/
/*
template<uint32_t dim>
void Fracture<dim>::PropagateTip(FractureTip<dim>* f_tip, Boundary<dim>& b_ref){

      //getting relevant tipnode
      Node<dim>* tipnode = f_tip->TipNode();


      //1) Create list of fracture faces lying on input boundary and the nodes within them which must be split
      std::map<Face<dim>*, std::set<size_t> > faces_to_split;
      for (typename std::vector<Face<dim>*>::const_iterator fit = b_ref.CellsBegin();  fit != b_ref.CellsEnd(); ++fit){
        assert( (*fit)->FE()->Dim() != 3U);

        //search if top node is in Face object
        for (size_t i = 0; i < (*fit)->Nodes(); ++i){
          if ( (*fit)->N(i) == tipnode ){
            std::set<size_t> nodes_to_split;
            nodes_to_split.insert(i);                         // Add tip node

            //Search for mid side nodes to split
            // Convention - neighbouring mid side nodes are also split
            if ( (*fit)->FE()->MidSideNodes() != 0 ){
              //get id's of all nodes on midside
              std::vector<size_t> midsidenodes;
              (*fit)->FE()->MidSideNodes(midsidenodes);

              //iterate over all segments and search for segment with tip node
              for (size_t s = 0; s < (*fit)->FE()->Segments(); ++s){
                //get segment nodes
                std::vector<size_t> segmentnodes;
                (*fit)->FE()->NodesOfSegment(s, segmentnodes);

                //if we find tip node in segment s
                if (std::find(segmentnodes.begin(), segmentnodes.end(), i) != segmentnodes.end() ){
                  //Loop over all nodes of segment
                  for (size_t n_s = 0; n_s < segmentnodes.size(); ++n_s ){
                    //check if segment node is a mid side node
                    typename std::vector<size_t>::iterator found_mid_side = std::find( midsidenodes.begin(), midsidenodes.end(), segmentnodes[n_s] ) ;
                    if (found_mid_side != midsidenodes.end() ){
                      //if segment node is midside node, then we split it
                      nodes_to_split.insert( *found_mid_side ) ;
                    }
                  } // end of loop over segment nodes
                } // end of if segment has tip node
              } // end of loop over all segments
            } // end of if mid side nodes

            // If also interior nodes exist on face
            if ((*fit)->FE()->InteriorNodes() != 0){
              std::vector<size_t> interiornodes;
              (*fit)->FE()->InteriorNodes(interiornodes);
              for (size_t in = 0; in < interiornodes.size(); ++in){
                nodes_to_split.insert(interiornodes[in]);         //inserting all possible interior nodes
              }
            } // End of if interior nodes

            //Inserting to map: Face, and corresponding set of nodes which need splitting
            faces_to_split.insert( std::make_pair(*fit, nodes_to_split) );
            break;
          } //end of if tip node found in face
        } //end node loop
      }//end face loop



      //2) Splitting nodes and creating new interfaces
      std::vector<InterFace<dim>*> new_interfaces = model_->SplitNodeOnFaces( faces_to_split, sb_ref_.Name(), true);

      //3) Getting new fracture half-length
      double fracture_length = FractureLength(MIDDLE);
      if (two_tips_)
        fracture_length *= 0.5;

      //4) Updating tip with new interfaces
      assert(dim == 2);
      assert(new_interfaces.size() == 1 );
      f_tip->UpdatePropagatedFractureTip(new_interfaces.front(), fracture_length);


} // end of PropagateTips()






*/







/**
  Takes properties defined on the node or element of the outside and inside of interface, and averages
  them to the nodes/elements of the middle region, respectively.

  This method preserves the flag of the side elements and requires the flag Outside and Inside to be the same
 */
template<uint32_t dim>
void Fracture<dim>::AveragePropertyToMiddle( Index property, std::vector<InterFace<dim>*>& interfaces ){

  if (property.type == SCALAR ){
    ScalarVariable sc_in, sc_out ;
      if (property.place == ELEMENT){
        for (size_t i = 0; i < interfaces.size(); ++i){
          //reading elemnts out and in
          interfaces[i]->Parent(INSIDE)->Read(property, sc_in);
          interfaces[i]->Parent(OUTSIDE)->Read(property, sc_out);

          //if flags are not the same then this method is ambiguous!
          assert( sc_in.Flag() == sc_out.Flag() );

          //changing sc_in to have average value but preserve flag
          sc_in.Component(0, 0.5 * ( sc_in() + sc_out() ) );

          //storing scalar variable to mid elmt
          interfaces[i]->Parent(MIDDLE)->Store( property, sc_in);
        }
      } else if (property.place == NODE){
        for (size_t i = 0; i < interfaces.size(); ++i){
          //check same node numbers, otherwise this method should not be used
          assert( interfaces[i]->FE()->Nodes() == interfaces[i]->Parent(MIDDLE)->Nodes() );

          std::vector<ScalarVariable> sc_vec_in, sc_vec_out;
          interfaces[i]->NodePropertyVector(property, sc_vec_in,  INSIDE);
          interfaces[i]->NodePropertyVector(property, sc_vec_out, OUTSIDE);

          //assigning average of nodes on each side to interface
          for ( uint32_t n = 0u; n < interfaces[i]->FE()->Nodes(); ++n){
            assert(sc_vec_in[n].Flag() == sc_vec_out[n].Flag());
            interfaces[i]->N(n,MIDDLE)->Store( property, ScalarVariable(sc_vec_in[n].Flag(), 0.5* (sc_vec_in[n]() + sc_vec_out[n]() ) ));
          }
        }
      } else
        throw csmp::Exception(ERROR, "Fracture<dim>::AveragePropertyToMiddle()",
                              "Property placement not supported!");
  } else if (property.type == VECTOR){

    if (property.place == ELEMENT){
        for (size_t i = 0; i < interfaces.size(); ++i){
          VectorVariable<dim> vc_in, vc_out;
          interfaces[i]->Parent(INSIDE)->Read(property, vc_in);
          interfaces[i]->Parent(OUTSIDE)->Read(property, vc_out);

          //check first and last flag are the same
          assert( vc_in.Flag(0) == vc_out.Flag(0) );
          assert( vc_in.Flag(dim - 1) == vc_out.Flag(dim - 1 ));

          //take average and preserve flag
          vc_in += vc_out;
          vc_in /= 2.0;

          //storing now averaged vector variable
          interfaces[i]->Parent(MIDDLE)->Store(property, vc_in);
        }
      } else if (property.place == NODE){
        for (size_t i = 0; i < interfaces.size(); ++i){
          //check same node numbers, otherwise this method should not be used
          assert( interfaces[i]->FE()->Nodes() == interfaces[i]->Parent(MIDDLE)->Nodes() );

          //reading vector nodal variables
          std::vector<VectorVariable<dim>> vc_vec_in, vc_vec_out;
          interfaces[i]->NodePropertyVector(property, vc_vec_in,  INSIDE);
          interfaces[i]->NodePropertyVector(property, vc_vec_out, OUTSIDE);

          for ( uint32_t n=0u; n < interfaces[i]->FE()->Nodes(); n++){
            VectorVariable<dim> vc = vc_vec_in[n] + vc_vec_out[n] ;
            vc /= 2.0;
            interfaces[i]->N(n,MIDDLE)->Store(property, vc);
          }
        }
      } else
        throw csmp::Exception(ERROR, "Fracture<dim>::AveragePropertyToMiddle()",
                              "Property placement not supported");
    } else
    throw csmp::Exception(ERROR, "Fracture<dim>::AveragePropertyToMiddle()",
                          "Property type not supported");


} //end of AveragePropertyToMiddle


















//Calculated the diffence in displacement between opposite nodes and stores it as the aperture
//This is correct if no existing aperture is already present.
template<uint32_t dim>
void Fracture<dim>::StoreDisplacementDifferenceAsAperture(const std::string displacement, const std::string aperture, INTERFACE_SIDE side ){
  assert( !sb_ref_.CellVector().empty());
  if (side == MIDDLE and midregion_->NodeVector().empty() )
      throw csmp::Exception(ERROR, "Fracture<dim>::StoreDisplacementDifferenceAsAperture",
                            "Middle Region doesn not exist, Therefore cannot store aperture on MIDDLE" );

  csmp::Index d_key = model_->Database().StorageKey(displacement.c_str());
  csmp::Index a_key = model_->Database().StorageKey(aperture.c_str());

  std::vector<Node<dim>*> counted_nodes;
  counted_nodes.reserve( sb_ref_.Cells() );

  for ( InterFace<dim>* IF : sb_ref_.CellVector()  )
      for ( uint32_t i_n = 0u; i_n < IF->FE()->Nodes(); i_n++){
          if ( std::find(counted_nodes.begin(), counted_nodes.end(), IF->MatchingN(i_n,side) ) == counted_nodes.end() ){

              //taking node on either side
              Node<dim>* node_plus (IF->MatchingN(i_n,OUTSIDE));
              Node<dim>* node_minus(IF->MatchingN(i_n,INSIDE));

              //reading displacement key
              VectorVariable<dim> d_plus, d_minus, normal;
              node_plus->Read(  d_key, d_plus);
              node_minus->Read( d_key, d_minus);

              normal = IF->UnitNormal(INSIDE);
              double implicit_aperture = normal.DotProduct( d_plus - d_minus );           //Calculating Implicit aperture without using node coordinates!
              IF->MatchingN(i_n, side )->Store( a_key, ScalarVariable(PLAIN,   implicit_aperture));   // updates the aperture property for latest apertures
              counted_nodes.push_back(IF->N(i_n,side));
        }
      }
} // end of StoreDisplacementDifferenceAsAperture



template<uint32_t dim>
double Fracture<dim>::DistanceFromTip(std::map<Point<dim>, Node<dim>*>& nodemap, Node<dim>* N, BOX_BOUNDARY side){   //distance from rightmost point

  switch (side){
    case RIGHT: {
        double lengthfromtip = 0.0;
        for ( auto i_n = nodemap.find(N->Coordinate()); i_n != --(nodemap.end()); ){
            auto i_n_old = i_n++;                                         //assign old value and then iterate
            assert(i_n->first == i_n->second->Coordinate());
            lengthfromtip += i_n_old->first.DistanceTo(i_n->first);       //add distance between each Node
          }
        return lengthfromtip;
      }
      break;
    case LEFT: {
        double lengthfromtip = 0.0;
        for ( auto i_n = nodemap.find(N->Coordinate()); i_n != ++(nodemap.begin()); ){
            auto i_n_old = i_n--;                                         //assign old value and then iterate
            assert(i_n->first == i_n->second->Coordinate());

            lengthfromtip += i_n_old->first.DistanceTo(i_n->first);       //add distance between each Node
          }
        return lengthfromtip;
      } break;
    default: {
      throw csmp::Exception(ERROR, "Fracture::DistanceFromTip()", "Not working for other box boundaries");
      } break;
    } //end of switch
}


//E.P Iterates over all interface objects and their nodes (or directly over the nodes) of the split boundary, returning the maximum value found of the input property
template<uint32_t dim>
double Fracture<dim>::MaxPropertyValue(const char* property, INTERFACE_SIDE side){

  Index ap_key =  model_->Database().StorageKey(property);

  assert(ap_key.place == NODE); //functionality only for nodes
  assert(ap_key.type == SCALAR ); //can only read  scalar properties

  double mx_prop(0.0);
  if (side != MIDDLE){
    for (typename std::vector<InterFace<dim>*>::const_iterator if_it = sb_ref_.CellsBegin();
         if_it != sb_ref_.CellsEnd(); if_it++){
      for ( uint32_t i = 0u; i < (*if_it)->FE()->Nodes(); i++){
        mx_prop = std::max(mx_prop, (*if_it)->N(i,side)->Read(ap_key)); //taking maximum found
      }
    }
  } else  //quicker process for middle region as no nodes are counted twice
    for (typename std::vector<Node<dim>*>::const_iterator n_it = midregion_->NodesBegin();
         n_it != midregion_->NodesEnd(); n_it++){
      mx_prop = std::max(mx_prop, (*n_it)->Read(ap_key) );
    }

  return mx_prop;
}


//Property aperture must be defined
template<uint32_t dim>
double Fracture<dim>::Volume(const char* aperture){
  double volume = sb_ref_.SurfaceIntegral( model_->Database(), aperture , MIDDLE);
  return volume;
}

template<uint32_t dim>
double Fracture<dim>::SurfaceIntegral( const char* op, INTERFACE_SIDE side){
  return sb_ref_.SurfaceIntegral( model_->Database(), op, side);
}


template<uint32_t dim>
void Fracture<dim>::ConfigureAnalyticalParametersAndSolutions(double constant_flux, double viscosity,
                                                  double critSIF, double youngs_modulus,
                                                  double poisson_ratio, bool plane_strain){

  configured_   = true;
  plane_strain_ = plane_strain;

  //configuring parameters for analytical solutions and fracture touhgnesses
  Q_  = constant_flux;
  mu_ = viscosity;
  Kc_ = critSIF;
  ym_ = youngs_modulus;
  nu_ = poisson_ratio;

}


/**
	Returns the dimensionless toughness using the scaling paramters youngs modulus, poisson ratio, 
	critical stress intensity factor, viscosity and constant flux, for a plane-strain fracture propagating in an infinite medium
*/
template<uint32_t dim>
double Fracture<dim>::DimensionlessToughness(){

  if (configured_ == false)
    throw csmp::Exception(ERROR,"Fracture<dim>::DimensionlessToughness",
                          "Configure Analytical Parameters first!");

  //Calcuulating
  double k_prime = std::sqrt( 32.0/M_PI ) * Kc_ ;
  double mu_prime =  12.0* mu_;
  double ym_prime(0.0);
  if (plane_strain_)
    ym_prime = ym_ / (1.0 - nu_*nu_);
  else
    throw csmp::Exception(ERROR, "Fracture::DimensionlessToughness", "plane stress or 3d sol not implemented");

  return k_prime / (std::pow ( mu_prime * std::pow(ym_prime, 3.0) * Q_ , 0.25 ) );

} //end of Dimensionless toughness



template<uint32_t dim>
double Fracture<dim>::DimensionlessViscosity()
{

  if (configured_ == false)
    throw csmp::Exception(ERROR,"Fracture<dim>::DimensionlessToughness",
                          "Configure Analytical Parameters first!");

  double K = DimensionlessToughness();
  double M = 1.0/(K*K*K*K);       //M = K ^-4

  return M;
}


/** EP
 Analytical half length for a plane strain griffith crack propagating in an infinite elastic medium
 Dahi taleghani thesis pg 130
 "Analysis of hydraulic fracture propagation in fractured reservoirs:
  an improved model for the interaction between induced and natural fractures"

  Toughness Dominated Solution: Garagash 2007 -  Plane Strain propagation of Fracture Durign injection and shut in: Asymptotics of large toughness.
  Small Viscosity Solution: Garagash 2007 -  Plane Strain propagation of Fracture Durign injection and shut in: Asymptotics of large toughness.

  Viscous-dominated solution: Garagash 2005 - Plane-Strain Propagation of a Fluid-Driven Fracture: Small Toughness Solution
*/

template<uint32_t dim>
double Fracture<dim>::AnalyticalAperture(double t, double x){

  if (configured_ == false)
    throw csmp::Exception(ERROR,"Fracture<dim>::AnalyticalAperture()",
                          "Configure Analytical Parameters first!");

  assert(plane_strain_ == true);

  if ( x > 1.001*AnalyticalLength(t) )
    throw csmp::Exception(ERROR, "Fracture<dim>::AnalyticalAperture(double t, double x)",
                          "Along fracture coordinate x is greater then fracture half length!");
  //convenient constants
  double  Ep = ym_ / (1. - nu_*nu_);
  double  mup = 12.0*mu_;
  double  Kp  = 4.0*std::pow( 2.0/M_PI, 0.5) * Kc_;

  //scaled coordinate
  double  xi = x / AnalyticalLength( t );

  /// Analytical solutions based on propagation regime
  //Viscosity Dominated Solution
  if (DimensionlessViscosity() > M_dominant){
    //Garagash 2005
//    double epsi_K = 0.1076 * std::pow( DimensionlessToughness(), 3.16796);
//    double OmegaBar_1_at_well = 0.75802;
    //std::cout << epsi_K*OmegaBar_1_at_well << std::endl;

    double gamma0     = 0.6152;
    double w_star     = std::pow( mup*Q_*Q_*Q_/Ep, 1.0/6.0) * std::pow(t, 1.0/3.0);
    double OmegaBar0  = 1.73205*std::pow( 1.0-xi*xi,2.0/3.0) - 0.15601*std::pow(1.0-xi*xi,5.0/3.0)
                          + 0.13264*(2.0*std::sqrt(1.0-xi*xi) );
    if (xi != 0.0)  //returns nan otherwise
      OmegaBar0        += 0.13264*xi*xi*std::log(std::fabs( (1.0 - std::sqrt(1-xi*xi))/(1.0+std::sqrt(1.0-xi*xi)) )) ;

    double aperture   = w_star * gamma0 * OmegaBar0;
    return aperture;

    //taleghani thesis
    //double Q0 = Q_/2.0;
    //return    2.36 * std::pow( mu_* Q0*Q0*Q0 / Ep , 1.0/6.0) * std::pow(t , 1.0/3.0) * std::pow( 1.0 - xi, 2.0/3.0) ;
  }
  //Small Toughness Solution
  else if (DimensionlessToughness() < K_small  ){
    throw csmp::Exception(ERROR, "NEED TO DO THIS REGIME Still", " ");
  }
  // Toughness Dominated Solution - Zero fluid viscosity and uniform fluid pressure
  else if (DimensionlessToughness() > K_dominant){
    double gamma    = 0.93238815407;
    double w_star   = std::pow( t * Kp*Kp * Q_ / (Ep*Ep) , 1.0/3.0 ) ;
    double OmegaBar_0  = 0.7322959437807616 * std::sqrt( 1.0 - xi*xi);
    double aperture = w_star * gamma * OmegaBar_0;
    return aperture;
  }
  // Small Viscosity Correction
  else if ( DimensionlessViscosity() < M_small ){  //See garagash 2007
    //first order terms
    double gamma_0    = 0.93238815407;
    double w_star     = std::pow( t * Kp*Kp * Q_ / (Ep*Ep) , 1.0/3.0 ) ;
    double OmegaBar_0 = 0.7322959437807616 * std::sqrt( 1.0 - xi*xi);
    //adding correction terms
    double M          = DimensionlessViscosity();
    double delta_M    = M / std::pow( 1.0 + M/0.0333 , 0.5);
    double gamma_1    = -2.722;
    double gamma      = gamma_0 + delta_M * gamma_1;
    double OmegaBar_1 = 1.243184205* (2.0*M_PI - 4.0*xi*std::asin(xi) - (5.0/6.0 - std::log(2.0) ) * std::sqrt(1.0 - xi*xi)
                                        - 1.5*std::log( std::pow( 1.0+std::sqrt(1.0-xi*xi) , 1.0+std::sqrt(1.0-xi*xi) ) / std::pow(1.0-std::sqrt(1.0-xi*xi), 1.0-std::sqrt(1.0-xi*xi))   ) );

    double aperture   = w_star * gamma * (OmegaBar_0 + delta_M * OmegaBar_1);
    return aperture;
  } else
    throw csmp::Exception(ERROR, "Fracture<dim>::AnalyticalAperture(double, double)", "Solution outside of asymptotic regimes!");
} //end of AnalyticalAperture()



/**
 Analytical fluid pressure are the wellbore for a plane strain griffith crack propagating in an infinite elastic medium
*/
template<uint32_t dim>
double Fracture<dim>::AnalyticalPressure(double t, double x)
{
    if ( !configured_ )
        throw csmp::Exception(ERROR,
            "Fracture<dim>::AnalyticalPressure()",
            "Configure Analytical Parameters first!");

    assert(plane_strain_);

    if ( x >= 1.0000001 * AnalyticalLength(t) )
        throw csmp::Exception(ERROR,
            "Fracture<dim>::AnalyticalPressure(double t, double x)",
            "Along fracture coordinate x is greater than fracture half length!");

    // Convenient constants
    const double Ep  = ym_ / (1.0 - nu_*nu_);
    const double mup = 12.0 * mu_;
    const double Kp  = 4.0 * std::sqrt(2.0/M_PI) * Kc_;
    const double xi  = x / AnalyticalLength(t);

    const double Mv = DimensionlessViscosity();
    const double Kt = DimensionlessToughness();

    // -----------------------------------------------------------------------
    // Viscosity dominated (with finite toughness correction)
    // -----------------------------------------------------------------------
    if ( Mv > M_dominant || Kt < K_small )
    {
        const double p_star = std::pow( Ep*Ep * mup / t, 1.0/3.0 );

        constexpr double b01 = 0.475449, b02 = -0.061178, b03 =  0.066322;
        constexpr double b11 = 0.170654, b12 =  0.017132, b13 = -0.039015,
                         b14 = -0.045476;
        // TODO: unused? - constexpr double c11 = 0.36133,  c12 = -1.63867,  c13 = -0.638673;

        const double xi2 = xi * xi;

        // TODO: this code issues warnings because without library, we do not know the return type of the macros 'gsl_sf_hyperg(,,,)' etc.
        const double PI_0 =
              b01 * gsl_sf_hyperg_2F1(-1.0/6.0, 1.0, 0.5,  xi2)
            + b02 * gsl_sf_hyperg_2F1(-7.0/6.0, 1.0, 0.5,  xi2)
            + b03 * (2.0 - M_PI * std::fabs(xi));

        const double PI_1 =
              b11 * gsl_sf_hyperg_2F1( c11,      1.0, 0.5,  xi2)
            + b12 * xi2 * gsl_sf_hyperg_2F1( c12, 1.0, 0.5, xi2)
            + b13 * gsl_sf_hyperg_2F1( c13,      2.0, 1.5,  xi2)
            + b14 * (2.0 - M_PI * std::fabs(xi));

        const double epsi_K = 0.1076 * std::pow(Kt, 3.16796);

        // Verify correction is small relative to leading term
        if ( PI_0 <= epsi_K * PI_1 )
            throw csmp::Exception(ERROR,
                "Fracture<dim>::AnalyticalPressure()",
                "Toughness correction exceeds leading term — outside valid range.");

        return p_star * (PI_0 + epsi_K * PI_1);
    }

    // -----------------------------------------------------------------------
    // Toughness dominated — zero viscosity, uniform pressure
    // -----------------------------------------------------------------------
    if ( Kt > K_dominant )
    {
        const double p_star = std::pow( Kp*Kp*Kp*Kp / (Ep * Q_ * t), 1.0/3.0 );
        constexpr double Pi_0 = 0.1830739859451904;
        return p_star * Pi_0;
    }

    // -----------------------------------------------------------------------
    // Small viscosity correction (Garagash 2007)
    // -----------------------------------------------------------------------
    if ( Mv < M_small )
    {
        const double p_star = std::pow( Kp*Kp*Kp*Kp / (Ep * Q_ * t), 1.0/3.0 );
        constexpr double Pi_0 = 0.1830739859451904;

        const double delta_M = Mv / std::sqrt(1.0 + Mv / 0.0333);
        const double Pi_1    = 1.243184205
                             * (  1.0/24.0
                                + std::log(4.0 * std::sqrt(1.0 - xi*xi))
                                - 3.0*xi * std::acos(xi)
                                  / (4.0 * std::sqrt(1.0 - xi*xi)) );

        return p_star * (Pi_0 + delta_M * Pi_1);
    }

    // -----------------------------------------------------------------------
    // Outside all asymptotic regimes
    // -----------------------------------------------------------------------
    throw csmp::Exception(ERROR,
        "Fracture<dim>::AnalyticalPressure(double, double)",
        "Solution outside of asymptotic regimes!");
}


/**
 Analytical half length for a plane strain griffith crack propagating in an infinite elastic medium
*/
template<uint32_t dim>
double Fracture<dim>::AnalyticalLength(double t){

  if (configured_ == false)
    throw csmp::Exception(ERROR,"Fracture<dim>::AnalyticalLength()",
                          "Configure Analytical Parameters first!");

  assert(plane_strain_ == true);    //Solutions only valid for plane strain propagating fracture

  //convenient constants
  double  Ep = ym_ / (1. - nu_*nu_);
  double  mup = 12.0*mu_;
  double  Kp  = 4.0*std::pow( 2.0/M_PI, 0.5) * Kc_;

  /// Analytical solutions based on propagation regime
  //Viscosity Dominated Solution
  if (DimensionlessViscosity() > M_dominant){
    double gamma0   = 0.61524;
    double gamma1   = -0.01887;
    double epsi_K   = 0.1076 * std::pow(DimensionlessToughness(), 3.16796 ) ;

    double L        = std::pow( Ep * Q_*Q_*Q_ / mup , 1.0/6.0) * std::pow(t, 2.0/3.0);

    double length   = L * ( gamma0 + epsi_K * gamma1);                                                               //garagash solution

    //should only be correction
    assert( epsi_K * gamma1 < gamma0 );

    return length;
    //talegahani thesis solution
    //double Q0 = Q_/2.0;
    //std::cout << "Length tale:\t" << 0.539 * std::pow( Ep * Q0*Q0*Q0 / mu_ , 1.0/6.0) * std::pow(t, 2.0/3.0) ;

  }
  //Small Toughness Solution
  else if (DimensionlessToughness() < K_small ){ //TODO: !!
    throw csmp::Exception(ERROR, "NEED TO DO THIS REGIME Still", " ");
  }
  // Toughness Dominated Solution - Zero fluid viscosity and uniform fluid pressure
  else if (DimensionlessToughness() > K_dominant){
    double gamma_0  = 0.93238815407;
    double L      = std::pow( t * Ep * Q_ / Kp, 2.0/3.0);

    double length = L * gamma_0 ;
    return length;
}
  // Small Viscosity Correction
  else if ( DimensionlessViscosity() < M_small ){             //See garagash 2007
    double gamma_0  = 0.93238815407;
    double L        = std::pow( t * Ep * Q_ / Kp, 2.0/3.0);
    double M        = DimensionlessViscosity();
    double delta_M  = M / std::pow( 1.0 + M/0.0333 , 0.5);
    double gamma_1  = -2.7220;

    double length   = L * ( gamma_0 + delta_M * gamma_1);
    return length;

  } else
    throw csmp::Exception(ERROR, "Fracture<dim>::AnalyticalLength(double)", "Solution outside of asymptotic regimes!");



}//end of analyticallength()





#ifdef DONT_READ

/* E.P TODO: Depracate functionality -- all creation done in Splitboundaryinterface now.

/// InitializeInterFaceObjects()
/// Initializes InterFaceObjects from a set of face pairs. The Faces represent the opposite faces of the interface,
/// which do not share any nodes. I.e the nodes must already be duplicated.
///
/// @brief Implementation is done using the InterFace(Face<dim>, Face<dim>) constructor, which should assign
///   i)   Parent Inner and Outer elements
///   ii)  Face FiniteElement pointer (which should be the same for both faces!)
///   iii) The number of InterFace Neighbors which InterFace has (equal to the number of Neighbours Face<dim> has).
///         note: The InterFace Neighbours are not actually assigned, just the space in interface_connector
///   iv)  outer and inner_parent_face_id, to localise elements face which is an InterFace.
///   v)   size of parent_elements_node_connector_ (just the size is assigned, not the local index of node on face)
///
/// After the construction of InterFace object, method assigns
///   i)   a local interface id,
///   ii)  parents_elements_node_connector (The elements local ID for opposite nodes )
///
/// Minimal construction of InterFace object is complete, and stored in MeshManager.
/// @note No InterFace connectivity is done!! Since SplitBoundary needs to do this!
template<uint32_t dim>
void Fracture<dim>::InitializeInterFaceObjects(){



   const LocalVariables             lvs     = model_->Database().LocalVariablesAt(INTER_FACE);
   const IntegrationPointVariables  lvs_int = model_->Database().IntegrationPointVariablesAt(INTER_FACE);

    // Making InterFace object vector from Faces
    interFacePointVec_.resize(facePairs_.size());
    size_t interface_i = 0;
    for (typename std::map< Face<dim>*,Face<dim>* >::iterator fmap_iter = facePairs_.begin(); fmap_iter!=facePairs_.end(); fmap_iter++){
        size_t face_nodes = fmap_iter->first->Nodes();
        assert( face_nodes == fmap_iter->second->Nodes() );   //total number of nodes on face must be the same

        //Constructing InterFace object from Faces
        InterFace<dim> IF(fmap_iter->first, fmap_iter->second, lvs, lvs_int) ;
        IF.Idx(interface_i);                                          //assigning local indexing

        //Assigning information of local Element node ID's to target nodes on either side of the facture-> Calibrates InterFace::parent_elements_node_connector_
        for (size_t i_node = 0; i_node < face_nodes; i_node++){
            typename std::map<Node<dim>*, Node<dim>*>::iterator nodepair_iterator =  nodePairs_.find( fmap_iter->first->N(i_node)) ;

            //obtaining local ID of Nodes parent element corresponding to inner element of the Face.
            size_t inner_parent_elem_number = -1;
            size_t outer_parent_elem_number = -1;
            for (size_t i_n_parents=0; i_n_parents < nodepair_iterator->first->Parents(); i_n_parents++ ){   //in current mesh, nodes can have up to 4 parents and its not guaranteed opposite nodes share same number of parents.
                if (nodepair_iterator->first->Parent(i_n_parents) == fmap_iter->first->InnerParent())
                  inner_parent_elem_number = i_n_parents; //finding the inner parents index of the node in the node pair
              }
            for (size_t i_n_parents=0; i_n_parents < nodepair_iterator->second->Parents(); i_n_parents++){
                if (nodepair_iterator->second->Parent(i_n_parents) == fmap_iter->second->InnerParent())
                  outer_parent_elem_number = i_n_parents;
              }

            if (inner_parent_elem_number == -1 || outer_parent_elem_number == -1)
              throw csmp::Exception( ERROR,
                                     "Fracture<dim>::Initialize",
                                     "Node Parent element not found to correspond with Face parent element. nodePair and Facepair not properly configured?");

            IF.Assign(i_node, std::make_pair(nodepair_iterator->first->ParentNodeNumber(inner_parent_elem_number) ,     //assigns the elements local node ID to the parent element connector in INterface
                                             nodepair_iterator->second->ParentNodeNumber(outer_parent_elem_number) ) );
        }

        interFacePointVec_[interface_i] = model_->Mesh().UniqueInsertInterFace(IF);     //collection of interface pointers
        assert(nodePairs_.find( IF.N(0,INSIDE))->second  == IF.N(0,OUTSIDE));
        assert(nodePairs_.find( IF.N(1,INSIDE))->second  == IF.N(1,OUTSIDE));
        assert(nodePairs_.find( interFacePointVec_[interface_i]->N(0,INSIDE))->second  == interFacePointVec_[interface_i]->N(0,OUTSIDE));
        assert(nodePairs_.find( interFacePointVec_[interface_i]->N(1,INSIDE))->second == interFacePointVec_[interface_i]->N(1,OUTSIDE));

        interface_i++;
      }
    //End of construction of InterFace object



}

*/


/* E.P TODO:  Use model functionality.
/// InitializeMidPointElements()
/// @brief Further Initalizes InterFace Object with a newly created MidElement which is assigned
///   i)   local element ID
///   ii)  InterFace Finite Element Pointer (so Element FE is lower dimensional)
///   iii) Element At_Boundary
///   iv)  New Nodes, which themselves are created and assigned
///           a) Coordinate
///           b) local index
///           c) At Boundary flag
///           d) Parent Elements (both number and pointer to element are assigned
///
/// MidPointMap_ is then created, which orders the elements in space, (but this is only used for postprocessing).
/// The newly created mid element is inserted into the MeshManager
///
/// @note Neighbour element connectivity must be assigned - Done by Region when calling Fracture<dim>::InitalizeMidPointRegion
template<uint32_t dim>
void Fracture<dim>::InitializeMidPointElements(){

  // Obtaining Interface pointer vector from Initialised split Boundary
  std::vector<InterFace<dim>*> IF_vec = model_->SplitBoundary(split_boundary_name_).CellVector();

  assert(!IF_vec.empty());
  assert(midPointNodes_.empty());

  boxboundaryNodes_.clear();
  //clearing size of Node and element container to be configured

  //making map container for unique nodes
  std::map< Point<dim>, Node<dim>* > uniqueNodes;
  size_t node_idx = 0;                                                                  //id wont be contiguous in space
  for (InterFace<dim>* interface_pointer : IF_vec ){                         //careful -- copying a pointer..?

      //create base element from face type, stored in interface_pointer
      Element<dim> baseElement(interface_pointer->FE());
      baseElement.AtBoundary(NOT);                                          //default to be changed unless node at perimiter found

      //assigning nodes to element
      size_t nodes_per_base_element = interface_pointer->InnerParent()->FE()->NodesPerFace(interface_pointer->InnerParentFaceID());
      for (size_t face_node_i = 0; face_node_i < nodes_per_base_element ; face_node_i++){
          Point<dim> midcoord = interface_pointer->MidpointCoordinateOfNode( face_node_i );              //taking midpoint coordinate of node
          typename std::map<Point<dim>, Node<dim>*>::iterator node_map_it = uniqueNodes.find(midcoord);  //finds the pointer to key with Point with the same coordinate
          if ( node_map_it == uniqueNodes.end() ) {                                                      //if node doesn't already exist with that midpoint coordinate

              Node<dim> midNode;                                                                        //creating node..., filling co-ordinate and adding to the map
              midNode.Coordinate( midcoord );                                                           //assigning coordinate
              midNode.Idx(node_idx);                                                                    //assign idx
              node_idx++;                                                                               //Iterate idx
              uniqueNodes[midcoord] = model_->Mesh().UniqueInsertNode(midNode);                         //Adds new coordinate and Node pointer
              midPointNodes_.push_back(uniqueNodes[midcoord]);                                          //Adds New node pointer to MidPointNodes_
              node_map_it = uniqueNodes.find(midcoord);                                                 //Iterator Points to new node in the unique node map

              assert(node_map_it != uniqueNodes.end());                                                 //check node is found in map

              //Assign Boundary flags to Nodes & Element
              node_map_it->second->AtBoundary(NOT);
              for (size_t node_i = 0; node_i < tipNodes_.size(); node_i++){
                  if ( node_map_it->second->Coordinate() == tipNodes_[node_i]->Coordinate() ){
                    node_map_it->second->AtBoundary(IRREGULAR);       //overwrites previous assignment of if node at boundary
                    baseElement.AtBoundary(IRREGULAR);
                    }
                }

              if ( interface_pointer->AtBoundary(face_node_i, INSIDE) != NOT &&
                   interface_pointer->AtBoundary(face_node_i, INSIDE) != IRREGULAR &&
                   interface_pointer->AtBoundary(face_node_i, INSIDE) != FRACTURE){       // If nodes of inner and outer parents are at a boundary (but not tip nodes)
                  BOX_BOUNDARY b_side = interface_pointer->AtBoundary(face_node_i,INSIDE);
                  assert(b_side == interface_pointer->AtBoundary(face_node_i,OUTSIDE));    // check both sides of interface are at a boundary

                  node_map_it->second->AtBoundary(b_side);                                 // Sets node also to be at Box Boundary of InterFace
                  boxboundaryNodes_.push_back(node_map_it->second);
                  baseElement.AtBoundary(b_side);
                }

            }

          baseElement.Assign(face_node_i, node_map_it->second);             //assigning each node (new or old) to the element

        } //end of nodes in element loop


      //Storing base Element, and assigning pointer to InterFace object
      midPointElementPtr_vec.push_back(model_->Mesh().UniqueInsertElement(baseElement));
      interface_pointer->Assign(midPointElementPtr_vec.back()) ;                // Assigns baseElement pointer to InterFace object

      Element<dim>* newEl = midPointElementPtr_vec.back();
      for (size_t i_n = 0; i_n < newEl->Nodes(); i_n++){
          newEl->N(i_n)->ResizeParentStorage(newEl->N(i_n)->Parents() + 1);     //Adding a parent each time a node is found by an element
          newEl->N(i_n)->Assign(i_n, newEl);                                    //Assigning parent element to each Node (must be done once element is stored in mesh)
        }

    } //End of InterFace obj and Element Creation


  //Calibrating Ordered Node Map
  SortMidRegionNodes();

  std::cout << "\nFracture<dim>::InitializeMidPointElementsMid\nMid point Elements initialized successfully with their nodes. "
               "However element neighbour connectivity yet to be assigned..." << std::endl;

} //End of InitializeMidPointElements



/// InitializeMidPointRegion
/// @brief Finalised integration of Mid Elements into a Region in model.
/// Creates Region from Region::Accumulate(std::vector<Element<dim>*>::iterator ...)
/// Then adds the Region to Model, calling it region, Indecis (indexing of nodes and Elements) are then rerun.
/// Then the new region is added to RegionInterface, at which point model is used to call the new Region, and establish the Element
/// connectivity via Region::EstablishNeighborConnectivity, which then allows the perimiter to be found.

template<uint32_t dim>
void Fracture<dim>::InitializeMidPointRegion(){
  assert(!midPointElementPtr_vec.empty());                            // Need to ensure that InitialiseMidPointElements is run first!

  //Creating Region from vector of element pointers
  Region<dim> fracturefluidflow(midregion_, model_->Database());
  fracturefluidflow.Accumulate(midPointElementPtr_vec.begin(), midPointElementPtr_vec.end());     //assigns elmt_vec and node_vec in ModelSubDomain -- no need for Region::CreateNodePointerVec()

  //Adding Region to "Model" region
  model_->Region("Model").Add(fracturefluidflow);                     // Node Pointer vector and Perimiter automatically called
  model_->UpdateIndices();                                            // renumbers indexing

  //Integrating new Region into Region Interface
  model_->FormRegionFrom(midregion_.c_str(), fracturefluidflow , true);           // This copies the fracturefluidflow region into RegionInterface

  //Calibrating Region in Region Interface
  Region<dim>& region_in_model = model_->Region(midregion_);
  establishNeighborConnectivity( region_in_model.CellVector() );   //Assigns Neighbours to elements by mesh traversal
  region_in_model.IdentifyPerimeter( );                               //Finds perimiter based on neighbours


  //Checking NodeVec ElementVec not equal to zero, and have elements  with nodes and parents
  assert(region_in_model.NodeVector().size() != 0);
  for (Node<dim>* n : region_in_model.NodeVector()){
      assert(n->Parents()!= 0);
    }
  assert(region_in_model.CellVector().size()!=0);
  for (Element<dim>* e : region_in_model.CellVector()){
      assert(e->Nodes()!= 0);
    }

  ///@todo Asserts missing, does ElementptrVec point to the same Elements pointed by in region_in_model ?

  std::cout << "Fracture<dim>::InitializeMidPointRegion \nNew Region '" << midregion_ <<
                "' successfully created and integrated into model." << std::endl;
}



template<uint32_t dim>
void Fracture<dim>::InitializeSplitBoundary(){
 assert(!interFacePointVec_.empty());

 SplitBoundary<dim> splitBoundary( split_boundary_name_, model_->Database());            //Calls split boundary constructor which uses database to calibrate InterFace LocalVariable storage depth
 splitBoundary.Accumulate(interFacePointVec_.begin(), interFacePointVec_.end());       // Initialises elem_vec and node_vec from interface Obj
 model_->FormSplitBoundaryFrom(split_boundary_name_ , splitBoundary);                   // Adds split boundary to SplitBoundaryInterface

 std::cout << "Fracture<dim>::InitializeSplitBoundary() -> Successfully initialised split boundary " << midregion_ << std::endl;
}


*/


/*
template<uint32_t dim>
double Fracture<dim>::SurfaceIntegral(const char* oper){
  model_->SplitBoundary(split_boundary_name_).SurfaceIntegral(model_->Database(), oper, true);
}


*/




//Gets elements on boundary of fluid flow region
//Non Intuititve way to obtain these elements!
//TODO: Implement like BoxBoundaryELements(side) is implemented.
template<uint32_t dim>
std::vector<Element<dim>*> Fracture<dim>::BoxBoundaryElements(){
  assert(model_->Region(midregion_).PerimeterElements() > 0);

  std::vector<Element<dim>*> bound_elements;
  for (Node<dim>* node : boxboundaryNodes_){
      size_t what =  node->Parents();


      if (node->Parents() == 0)
        throw csmp::Exception(ERROR, "Fracture<dim>::BoundaryElements()",
                              "Node doesnt have a parent! " );

      for (size_t parent = 0; parent < node->Parents(); parent++){            //Is this not working!? Make sure nodes on Box Boundary are labelled as such, and not as FRACTURE or other boundaries.
          assert( node->Parent(parent)->AtBoundary() != NOT);
          bound_elements.push_back(node->Parent(parent));
        }

  }

  auto last = std::unique(bound_elements.begin(), bound_elements.end());      // return pointer to the actual size of unique elements and rearrange existing container
  bound_elements.erase(last, bound_elements.end());

  return bound_elements;

  }


//Iterates through perimiter elements of the midRegion and takes the elements which are AtBoundary( side )
template<uint32_t dim>
std::vector<Element<dim>*> Fracture<dim>::BoxBoundaryElements( BOX_BOUNDARY side){

  assert(model_->Region(midregion_).PerimeterElements() > 0);

  std::vector<Element<dim>*> elem_vec;
  for (typename std::vector<Element<dim>*>::const_iterator e_ptr = model_->Region(midregion_).PerimeterElementsBegin() ; e_ptr != model_->Region(midregion_).ElementsEnd(); e_ptr++){
      if ((*e_ptr)->AtBoundary() == side )
        elem_vec.push_back(*e_ptr);
    }

  return elem_vec;
}











/// Setting Boundary Conditions 1D
/// WARNING: THIS ONLY FINDS THE BOUNDARIES BASED ON THEIR CO_ORDINATES, Not on the actual enum "bound_side".
/// TODO: change setting of boundary value based on search for Box_boundary enum! Intuitive approach!
template<uint32_t dim>
void Fracture<dim>::InputPropertyValue(const char* prop, ScalarVariable Sval, BOX_BOUNDARY bound_side){


  Node<dim>* bound_node = nullptr;

  switch (bound_side){
    case LEFT:
      {
        size_t i = 0;
        for (typename std::vector<Node<dim>*>::iterator it = model_->Region(midregion_).NodesBegin(); it != model_->Region(midregion_).NodesEnd(); it++ ){
          if ( (*it)->AtBoundary() != NOT){
            if (i++ == 0)
              bound_node = *it;
            else if ( (*it)->Coordinate() < bound_node->Coordinate() )
              bound_node = *it;                                                     //finding leftmost coordinate
            }
          }

      }
      break;
    case RIGHT:
      {
        size_t i = 0;
        for (typename std::vector<Node<dim>*>::iterator it = model_->Region(midregion_).NodesBegin(); it != model_->Region(midregion_).NodesEnd(); it++ ){
          if ( (*it)->AtBoundary() != NOT){
            if (i++ == 0)
              bound_node = *it;
            else if ( (*it)->Coordinate() > bound_node->Coordinate() )
              bound_node = *it;                                                     //finding rightmost coordinate
            }
          }
      }
      break;
    default: throw csmp::Exception(ERROR, "Fracture<dim>::InputPropertyValue(const char, double, BOX_BOUNDARY)",
                                   "Only use LEFT or RIGHT cases of BOX_BOUNDARY. For a purely vertical fracture the LEFT corresponds to the smallest y directio. "
                                   "Otherwise the leftmost tip is always assigned by LEFT");
      break;
    }

  assert(bound_node != nullptr);

  bound_node->Store(model_->Database().StorageKey(prop), Sval);

  std::cout << "Fracture<dim>::InputPropertyValue(const char*, double, BOX_BOUNDARY\n" << parseBoundary(bound_side) <<
               " Boundary assigned value and Flag" << std::endl;

}


template<uint32_t dim>
void Fracture<dim>::SortMidRegionNodes(){
  midPointMap_.clear();
  assert( !midPointNodes_.empty());
  for (Node<dim>* N : midPointNodes_){
      midPointMap_[N->Coordinate()] = N;
    }

}


/* IS THIS USED?
//TODO: SLOW because doing search of counted nodes for each InterFace!
template<uint32_t dim>
void Fracture<dim>::MoveMidPointCoordinatesBy(const char *displacement){
  assert( !model_->SplitBoundary(split_boundary_name_).CellVector().empty());

  csmp::Index d_key = model_->Database().StorageKey(displacement);
  VectorVariable<dim> d_plus, d_minus, average_displacement;

  std::vector<Node<dim>*> counted_nodes;
  counted_nodes.reserve(model_->Region(midregion_).NodeVector().size());

  for (InterFace<dim>* i_f : model_->SplitBoundary(split_boundary_name_).CellVector()){
      for (size_t node = 0; node < i_f->Nodes() ; node++){
          if ( std::find(counted_nodes.begin(), counted_nodes.end(), i_f->N(node,MIDDLE) ) == counted_nodes.end() ){
            i_f->N(node,OUTSIDE)->Read(d_key,d_plus);
            i_f->N(node,INSIDE)->Read( d_key, d_minus);

            average_displacement = (d_plus + d_minus) / 2.0;

            i_f->N(node, MIDDLE)->Coordinate( i_f->N(node,MIDDLE)->Coordinate() + average_displacement.P() );
            counted_nodes.push_back(i_f->N(node,MIDDLE));
            }
        }
    }
  SortMidRegionNodes(); // UPDATES midPointMap
}
*/

/*
//TODO: SLOW beacuse it does operations twice on nodes which are shared
template<uint32_t dim>
void Fracture<dim>::AverageCoordinatesToMidPoint(){
  assert( !model_->SplitBoundary(split_boundary_name_).CellVector().empty());
  for (InterFace<dim>* i_f : model_->SplitBoundary(split_boundary_name_).CellVector()){
      for (size_t node = 0; node < i_f->Nodes() ; node++){
          i_f->N(node, MIDDLE)->Coordinate(i_f->MidpointCoordinateOfNode(node));
        }
    }
  SortMidRegionNodes(); // UPDATES midPointMap
}
*/


//NEEDS database variable "aperture" and "aperture change" to be defined.
// CURRENTLY SLOW SINCE NEEDS TO DO A SEARCH ALGORITHM FOR EACH NODE!!
// this updates the value of aperture and aperture change given a known old aperture! Not to be used for TIME STEPPING!
template<uint32_t dim>
void Fracture<dim>::IterateApertures(){
  assert(!model_->SplitBoundary(split_boundary_name_).CellVector().empty());
  assert(!model_->Region(midregion_).NodeVector().empty());

  std::vector<Node<dim>*> counted_nodes;
  counted_nodes.reserve(model_->Region(midregion_).NodeVector().size());

  for ( InterFace<dim>* IF : model_->SplitBoundary(split_boundary_name_).CellVector()  ){
      for (size_t i_n = 0; i_n < IF->Nodes(); i_n++){
          if ( std::find(counted_nodes.begin(), counted_nodes.end(), IF->N(i_n,MIDDLE) ) == counted_nodes.end() ){

              double aperture         = IF->ApertureNode(i_n);
              double aperture_change  = IF->ApertureNode(i_n) -  IF->N(i_n,MIDDLE)->Read(model_->Database().StorageKey("aperture old"));

              IF->N(i_n, MIDDLE )->Store(model_->Database().StorageKey("aperture"),        ScalarVariable(PLAIN,  IF->ApertureNode(i_n)));   // updates the aperture property for latest apertures
              IF->N(i_n, MIDDLE )->Store(model_->Database().StorageKey("aperture change"), ScalarVariable(PLAIN,  aperture_change)) ;         // updates the aperture change given a known old aperture

              assert(aperture_change == aperture);  // only for first time step
              counted_nodes.push_back(IF->N(i_n,MIDDLE));
        }

      }
    }
}




//Takes the current displacement solution and applies it to the existing nodal coordinates to then calculate the theoretical aperture.
//Does not move nodes, or need nodes to have moved.
//WARNING make sure displacements coincide with the node co-ordinates from which they come from.
template<uint32_t dim>
void Fracture<dim>::ImplicitlyCalculateAndStoreApertures(const std::string aperture, const std::string displacement ){
  assert(!model_->SplitBoundary(split_boundary_name_).CellVector().empty());
  assert(!model_->Region(midregion_).NodeVector().empty());

  csmp::Index d_key = model_->Database().StorageKey(displacement.c_str());
  csmp::Index a_key = model_->Database().StorageKey(aperture.c_str());

  std::vector<Node<dim>*> counted_nodes;
  counted_nodes.reserve(model_->Region(midregion_).NodeVector().size());

  VectorVariable<dim> d_plus, d_minus, normal;

  for ( InterFace<dim>* IF : model_->SplitBoundary(split_boundary_name_).CellVector()  ){
      for (size_t i_n = 0; i_n < IF->Nodes(); i_n++){
          if ( std::find(counted_nodes.begin(), counted_nodes.end(), IF->N(i_n,MIDDLE) ) == counted_nodes.end() ){

              Node<dim>* node_plus (IF->N(i_n,OUTSIDE));
              Node<dim>* node_minus(IF->N(i_n,INSIDE));

              node_plus->Read(  d_key, d_plus);
              node_minus->Read( d_key, d_minus);

              VectorVariable<dim> Position_plus   = VectorVariable<dim>(node_plus->Coordinate())   +   d_plus;
              VectorVariable<dim> Position_minus  = VectorVariable<dim>(node_minus->Coordinate())  +   d_minus;
              IF->UnitNormal(normal, MIDDLE);

              double implicit_aperture = normal.DotProduct( Position_plus - Position_minus );    //Calculating Implicit aperture without moving nodes!

              IF->N(i_n, MIDDLE )->Store( a_key , ScalarVariable(PLAIN,   implicit_aperture));   // updates the aperture property for latest apertures

              counted_nodes.push_back(IF->N(i_n,MIDDLE));
        }

      }
    }
}


//Calculated the diffence in displacement between opposite nodes and stores it as the aperture
//This is correct if no existing aperture is already present.
template<uint32_t dim>
void Fracture<dim>::StoreDisplacementDifferenceAsAperture(const std::string displacement, const std::string aperture){
  assert(!model_->SplitBoundary(split_boundary_name_).CellVector().empty());
  assert(!model_->Region(midregion_).NodeVector().empty());

  csmp::Index d_key = model_->Database().StorageKey(displacement.c_str());
  csmp::Index a_key = model_->Database().StorageKey(aperture.c_str());

  std::vector<Node<dim>*> counted_nodes;
  counted_nodes.reserve(model_->Region(midregion_).NodeVector().size());

  for ( InterFace<dim>* IF : model_->SplitBoundary(split_boundary_name_).CellVector()  ){
      for (size_t i_n = 0; i_n < IF->Nodes(); i_n++){
          if ( std::find(counted_nodes.begin(), counted_nodes.end(), IF->N(i_n,MIDDLE) ) == counted_nodes.end() ){

              Node<dim>* node_plus (IF->N(i_n,OUTSIDE));
              Node<dim>* node_minus(IF->N(i_n,INSIDE));

              VectorVariable<dim> d_plus, d_minus, normal;
              node_plus->Read(  d_key, d_plus);
              node_minus->Read( d_key, d_minus);

              IF->UnitNormal(normal, MIDDLE);

              double implicit_aperture = normal.DotProduct( d_plus - d_minus );           //Calculating Implicit aperture without using node coordinates!

              IF->N(i_n, MIDDLE )->Store( a_key, ScalarVariable(PLAIN,   implicit_aperture));   // updates the aperture property for latest apertures

              counted_nodes.push_back(IF->N(i_n,MIDDLE));
        }

      }
    }
}




template<uint32_t dim>
void Fracture<dim>::ZeroOutNegativePressures(const std::string pressure, double min_press){
  assert(!model_->SplitBoundary(split_boundary_name_).CellVector().empty());
  assert(!model_->Region(midregion_).NodeVector().empty());

  for (Node<dim>* node : model_->Region(midregion_).NodeVector()){
      if (node->Read(model_->Database().StorageKey(pressure.c_str())) < min_press )
        node->Store(model_->Database().StorageKey(pressure.c_str()), ScalarVariable(PLAIN, min_press));
    }

  std::cout << "Fracture<dim>::ZeroOutNegativePressures() : Successfully eliminated negative pressures in SplitBoundary MidRegion." << std::endl;
}




template<uint32_t dim>
double Fracture<dim>::DistanceFromTip(Node<dim>* N, BOX_BOUNDARY side){   //distance from rightmost point
  assert(!midPointMap_.empty());

  switch (side){
    case RIGHT: {
        double lengthfromtip = 0.0;
        for (typename std::map< Point<dim>,Node<dim>* >::iterator i_n = midPointMap_.find(N->Coordinate()); i_n != --(midPointMap_.end()) ; 0){
            auto i_n_old = i_n++;                                         //assign old value and then iterate
            assert(i_n->first == i_n->second->Coordinate());
            lengthfromtip += i_n_old->first.DistanceTo(i_n->first);       //add distance between each Node
          }
        return lengthfromtip;
      }
      break;
    case LEFT: {
        double lengthfromtip = 0.0;
        for (typename std::map< Point<dim>,Node<dim>* >::iterator i_n = midPointMap_.find(N->Coordinate()); i_n != ++(midPointMap_.begin()) ; 0){
            auto i_n_old = i_n--;                                         //assign old value and then iterate
            assert(i_n->first == i_n->second->Coordinate());

            lengthfromtip += i_n_old->first.DistanceTo(i_n->first);       //add distance between each Node
          }
        return lengthfromtip;
      }
    }


}





template<uint32_t dim>
void Fracture<dim>::CalculateFlux(const char *flux, const char *pressure, const char* aperture, const char* viscosity){

  csmp::Index p_key  = model_->Database().StorageKey(pressure);
  csmp::Index f_key  = model_->Database().StorageKey(flux);
  csmp::Index mu_key = model_->Database().StorageKey(viscosity);
  csmp::Index w_key  = model_->Database().StorageKey(aperture);


  std::map< Point<dim>, Node<dim>*> ordered_nodes;
  for (Element<dim>* e_ptr : model_->Region(midregion_).CellVector()){

      assert(ordered_nodes.size()==0);

      std::vector<size_t> corner_node_ids;
      e_ptr->FE()->CornerNodes(corner_node_ids);                    //Obtains ID's of nodes which lie on the corner of element ( 0 1 returned for linear line)
      assert(corner_node_ids.size() == 2);                          //TODO: Only for 2D. Think of something better for 3d

      for (size_t i = 0; i < corner_node_ids.size() ; i++ ){
        ordered_nodes[ e_ptr->N(corner_node_ids[i])->Coordinate()] =  e_ptr->N(corner_node_ids[i]) ;    //Take corner nodes and insert them ordered by their coordinate
        }

      double dx     =        ordered_nodes.begin()->first.DistanceTo( (--(ordered_nodes.end()))->first);
      double dp     =        (--(ordered_nodes.end()))->second->Read(p_key)   -   ordered_nodes.begin()->second->Read(p_key) ;
      double mu_bar = 0.5 * ((--(ordered_nodes.end()))->second->Read(mu_key)  +   ordered_nodes.begin()->second->Read(mu_key) ) ;
      double w_bar  = 0.5 * ((--(ordered_nodes.end()))->second->Read(w_key)   +   ordered_nodes.begin()->second->Read(w_key)  ) ;

      //Calculalate flux through fracture
      double q      = - w_bar*w_bar*w_bar / (12.0*mu_bar) * dp/dx;

      //Store flux
      assert(f_key.place == ELEMENT);                       //correct evaluation of flux
      e_ptr->Store(f_key, ScalarVariable( PLAIN, q) );

      //Clear ordered nodes
      ordered_nodes.clear();

    }


}












/*
///method returns a map ordered such of Node Coordinates (ordered from the leftmost x co-ordintate to rightmost x coordinate)
///and the element of the map is the Aperture at the node, which is calculated on the fly by the method in the InterFace object
/// This method is useful to plot aperture and to compare to analytical solution.
template<uint32_t dim>
std::map<Point<dim>, double> Fracture<dim>::ApertureMap(){
  std::map<Point<dim>, double> coordApertureMap;
  assert(!model_->SplitBoundary(split_boundary_name_).CellVector().empty());
  for (InterFace<dim>* IF_Obj : model_->SplitBoundary(split_boundary_name_).CellVector()){
      for (size_t local_n = 0; local_n < IF_Obj->Nodes(); local_n++){
          coordApertureMap[IF_Obj->N(local_n, MIDDLE)->Coordinate()] = IF_Obj->ApertureNode(local_n);
          if (coordApertureMap[IF_Obj->N(local_n, MIDDLE)->Coordinate()] < 0.0)
            throw csmp::Exception(ERROR, "Fracture<dim>::ApertureMap()",
                                  "WARNING, Aperture calculated is negative! Check normal directions and Interpenetration!");
        }

    }
  return coordApertureMap;
}
*/



/*
template<uint32_t dim>
std::map<Point<dim>, double> Fracture<dim>::PressureMap(){
  std::map<Point<dim>, double> coordPressureMap;
  assert(!model_->SplitBoundary(split_boundary_name_).CellVector().empty());
  for (InterFace<dim>* IF_Obj : model_->SplitBoundary(split_boundary_name_).CellVector()){
      for (size_t local_n = 0; local_n < IF_Obj->Nodes(); local_n++){
          coordPressureMap[IF_Obj->N(local_n, MIDDLE)->Coordinate()] = IF_Obj->N(local_n, MIDDLE)->Read(model_->Database().StorageKey("fracture fluid pressure"));
        }

    }
  return coordPressureMap;
}
*/

template<uint32_t dim>
std::map<Point<dim>, Node<dim>*> Fracture<dim>::NodeMap(){
  assert(!midPointMap_.empty());
  return midPointMap_;
}






           /**
            This method generates node points which will be used to build a disk filled up with quadratic triangular and line elements.
            This domain will be used to compute domain integrals for fracture growth purposes. Refer to Nejati et al. (IJSS, 2015) for more details;
            */

#ifdef DONT_READ
template<uint32_t dim>
void Fracture<dim>::GeneratePolarStructuredMeshOverADisk (int nr, int nt)
{
    size_t npeT (6), npeL (3);
    DenseMatrix<DM_MIN> XYT_ ((2*nr-1)*nt*npeT,2,0.);  ///< contains the corner nodes' position of disk virtual triangular elements
    DenseMatrix<DM_MIN> XYL_ (nr*npeL,2,0.);          ///< contains the corner nodes' position of disk virtual line elements

    double r (domainRadius_/nr), t (PI/2./nt);

    for (size_t i(0); i < nr; ++i)
    {
        // Generating Surface Elements
        for (int j=0; j < nt; ++j)
        {
            size_t idx = (2*i+1)*j*npeT;

            const Point<2u> p1   (   i   * r *cos(  j     * t),   i   * r * sin(  j     * t));
            const Point<2u> p2   ( (i+1) * r *cos(  j     * t), (i+1) * r * sin(  j     * t));
            const Point<2u> p3   ( (i+1) * r *cos((j+1)   * t), (i+1) * r * sin((j+1)   * t));
            const Point<2u> p4   (   i   * r *cos((j+1)   * t),   i   * r * sin((j+1)   * t));
            const Point<2u> pm   ( (i+1) * r *cos((j+0.5) * t), (i+1) * r * sin((j+0.5) * t)); // Midside node on the perimeter of the disk

            XYT_.AssignRow(idx  , p1       ); // node 1
            XYT_.AssignRow(idx+1, p2       ); // node 2
            XYT_.AssignRow(idx+2, p3       ); // node 3
            XYT_.AssignRow(idx+3,(p1+p2)/2.); // node 4
            XYT_.AssignRow(idx+4,(p2+p3)/2.); // node 5
            XYT_.AssignRow(idx+5,(p3+p1)/2.); // node 6

            if (i>0)
            {
                XYT_.AssignRow(idx+7 , p1       ); // node 1
                XYT_.AssignRow(idx+8 , p2       ); // node 2
                XYT_.AssignRow(idx+9 , p4       ); // node 3
                XYT_.AssignRow(idx+10,(p1+p3)/2.); // node 4
                XYT_.AssignRow(idx+11,(p3+p4)/2.); // node 5
                XYT_.AssignRow(idx+12,(p4+p1)/2.); // node 6
            }

            if (i == 0 && quarterPoint_)
            {
                XYT_.AssignRow(idx+3,(p1*3.+p2)/4.); // node 4
                XYT_.AssignRow(idx+5,(p1*3.+p3)/4.); // node 6
            }

            if (i == nr-1 )
                XYT_.AssignRow(idx+4,pm);           // node 5 on perimeter
        }

        //Generating Line Elements
        //std::vector<Point<2U>> bar (2);
        const Point<2U> p1 ( - i    * r, 0.);
        const Point<2U> p2 ( -(i+1) * r, 0.);

        XYL_.AssignRow(i*npeL  ,  p1       );
        XYL_.AssignRow(i*npeL+1,  p2       );
        XYL_.AssignRow(i*npeL+2, (p1+p2)/2.);

        if (i == 0 && quarterPoint_)
            XYL_.AssignRow(i+2, (p1*3.+p2)/4.);
    }
}


template<uint32_t dim>
void Fracture<dim>::SolOut(){
  for (typename std::map<Point<dim>,Node<dim>*>::iterator it = midPointMap_.begin(); it != midPointMap_.end(); it++){
      std::cout << "\nPos: " << it->first[0] << ", " << it->first[1]
                << "\np        = " << it->second->Read(model_->Database().StorageKey("fracture fluid pressure"))
                << "\np old    = " << it->second->Read(model_->Database().StorageKey("pressure old"))
                << "\nw        = " << it->second->Read(model_->Database().StorageKey("aperture"))
                << "\nw old    = " << it->second->Read(model_->Database().StorageKey("aperture old"))
                << "\nw change = " << it->second->Read(model_->Database().StorageKey("aperture change"))
                << "\nq        = " << it->second->Parent(0)->Read(model_->Database().StorageKey("flux"))
                << " of point x = " << it->second->Parent(0)->N(0)->Coordinate()[0] << std::endl;

    }
}


template<uint32_t dim>
void Fracture<dim>::Out()
{
  for (typename std::map<Point<dim>,Node<dim>*>::iterator it = midPointMap_.begin(); it != midPointMap_.end(); it++){
      std::cout << "\nPos: " << it->first[0] << ", " << it->first[1] << " :  p = " << it->second->Read(model_->Database().StorageKey("fracture fluid pressure"))
                << "\nw        = " << it->second->Read(model_->Database().StorageKey("aperture"))
                << "\nw change = " << it->second->Read(model_->Database().StorageKey("aperture change")) << std::endl;
    }

  /*
std::cout << "\n\nFracture FACES: " << faces_.size();
for (typename std::vector<Face<dim>*>::const_iterator iter = faces_.begin(); iter != faces_.end(); iter++)
    (*iter)->Out();
std::cout << "\nFracture NODES: " << nodes_.size();
for (typename std::vector<Node<dim>*>::const_iterator iter = nodes_.begin(); iter != nodes_.end(); iter++)
    (*iter)->Out();
    */
}

#endif

#endif

#ifdef CSMP_WITH_SAMG_SOLVER

template<uint32_t dim>
void Fracture<dim>::SetSolverSettings(SAMG_Settings& settings)
 {


  std::cout << "\nsetSolverSettings";
  settings.Set_isym(2);		/*	1	A is symmetric    //Was at 2 !
                      2	A is not symmetric*/
  settings.Set_itypu(0);		/*	itypu	0 Actual content of u is chosen as first approximation
                      1	First approximation u==0
                      2	First approximation u==1
                      3	First approximation is a random function */
  settings.Set_eps(1.e-12);		/* =0.0	Stopping criteria based in eps is de-activated
                        >0.0	Iteration stops if res <= eps.res0 (res0 = frist residual)
                        <0.0	Iteration stops if res <= |eps| */
  settings.Set_rel_eps(1.e-11); /*Stopping criterion when first guess is set (itypu == 0).
                         As stopping criterion rel_eps * ||rhs|| will be used. Make sure that rel_eps
                         has a negative value, so that the absolute value of rel_eps is chosen as
                         convergence criterion. (See explanation in Set_eps())*/
  settings.Set_napproach(2);	/*	1	Scalar approach (regardless of nsys)
                      2	Unknown-based (if used in scalar system napproach will be reset to 1)
                      3-5	Point32 based approaches - selects type of interpolation to use
                      3: interp. is separate for each unknown
                      4: interp. is same " " "
                      5: interp. is point- (block-) wise*/
  settings.Set_nxtyp(0);      /*	0	Gauss-Seidel relaxation
                      1	ILU(0) (substantial increase in required memory)
                      2	ILUT
                      3	Special box relaxation
                      5	Gauss-Seidel blockwise
                      ILU - best for coupled problems? - but more expensive*/
  settings.Set_igam(1);		/*	igam		1 V-cycle (standard)
                      2	F-cycle
                      3	W-cycle
                      4	WW-cycle (very expensive)*/
  settings.Set_ncgrad(1);		/*	0	default accelerator
                      1	Preconditioner for CG (standard)
                      2	Precon. for BI-CGSTAB
                      3	Precon. for GMRES*/
  settings.Set_nkdim(0);		/*	Subswitch for ncyc.
                      0	Select default dimension
                      1-8	Dimension = nkdim+1
                      9	Dimension = 20*/

  settings.Set_iswit(5);		/*	iswit	Controls re-use of SAMG decompositions during repeated calls.
                      5 	Complete SAMG run. [Upon return, memory is released
                      4	Same as 5 except memory not released
                      3	partial setup: Re-use coarser grids and interp. but update
                      Galerkin operators. Memory not released
                      2	No setup: Re-use coarser grids, interp. and Galerkin from prev. run
                      1	Same as 2 except SAMG assumes matrix A to be the same as prev. run*/
  settings.Set_iextent(1);	/*		iextent	Memory extension switch. Selects beahaviour when limits of initial
                      dimensioning have been reached
                      0	SAMG returns with error code
                      1	SAMG allocates ext. memory and continues (if no core space,
                      writes prev. allocated data to disk
                      2	SAMG allocates ext. memory and continues (if no core space,
                      SAMG terminates)
                      3	SAMG allocates ext. memory and continues (prev. allocated data
                      is written to disk)  */
  settings.Set_norm_typ(0);	/*	norm_typ	Selects type of norm to be used in computing residuals
                      0	L2-norm
                      1	L1-norm
                      2	Maximum norm*/
  settings.Set_idmp(1);		/*  idmp   0 Coarsening history
                      1 Standard print output (coarsening history).*/
  settings.Set_igdp(1);		/* igdp  >1 is only relevant for coupled systems. Otherwise: ignored.
                      0 No particular output.
                      1 Display table on grids (full problem).
                      2 Same for all submatrices (only if nsys>1).*/
  settings.Set_iadp(1);		/*  >1 is only relevant for coupled systems.
                      0 No particular output.
                      1 Display table on coarse-level matrices (full problem).
                      2 In addition: same info for all submatrices.
                      3 In addition: connectivity info between unknowns.*/
  settings.Set_iwdp(1);		/*  iwdp   >1 is only relevant for coupled systems.
                      0 No particular output.
                      1 Display table on interpolation matrices (full problem).
                      2 In addition: same info for all submatrices.
                      3 In addition: connectivity info between unknowns.*/
  settings.Set_iout1(2);		/*	1 Table of input data and work statistics.
                      2 Standard history of cycling process.
                      3 Extended history: including all levels.
                      4 Extended history: including all levels and even partial smoothing*/
  settings.Set_iout2(0);		/*  0 No action.
                      1 Display all of SAMG’s hidden parameters.*/
  settings.Set_ncgtyp(1);		/*	1 Standard process. This is supposed to be used if A has
                      mostly negative off-diagonals. Positive off-diagonal
                      elements (if any) should be small. Variables with only
                      positive couplings will become C-variables.
                      2 Standard process except that variables which have only
                      positive couplings are treated by absolute value.
                      3 Standard process except that, for mixed-sign rows, all
                      "large" positive entries (threshold parameter ewt2) are
                      eliminated before a decision on strong connectivity is
                      made. If, for some variable i, this does not lead to a clear
                      decision, i will become a C-variable.
                      4 Same as 3 except that, if the elimination of positive
                      couplings of variable i fails to give a clear picture, this
                      option temporarily switches to the standard process 1.
                      5 Same as 4 except that variables which have only
                      positive couplings are treated by absolute value.*/
  settings.Set_nred(1);		/*	0 Standard coarsening.
                      1-4 Aggressive coarsening. 1 is most, 4 is least
                      aggressive;
                      2-3 are in between. Recommendation:
                      use 1 for anisotropic and 2 for isotropic problems.
                      5 Cluster coarsening & piecewise constant interpolation.
                      6 Cluster coarsening & multi-pass interpolation.*/
  settings.Set_ncycle(20000);

} // end of SetSettings()

#endif




template<uint32_t dim>
void Fracture<dim>::SolOut(INTERFACE_SIDE side, const char* aperture, const char* pressure)
{
  Index a_key = model_->Database().StorageKey(aperture);
  Index p_key = model_->Database().StorageKey(pressure);

  std::map<Point<dim>, Node<dim>*> nodemap = this->NodeMap(side);
  for (typename std::map<Point<dim>,Node<dim>*>::iterator it = nodemap.begin(); it != nodemap.end(); it++){
      std::cout << "\nPos: " << it->first[0] << ", " << it->first[1]
                << "\n\tp = " << it->second->Read(p_key)
                << "\n\tw = " << it->second->Read( a_key) << std::endl;
    }
}

/**
	E.P Tests made specific to a KGD fracture, which is propagating from the left boundary inwards. 
	It tests the split is done correctly by taking a copy of the interfaces before they are split, and afterwards and veryfying that
		-) Parents of interface objects have stayed the same (including middle)
		-) splitboundary grew by one
		-) Neighbors before are the same as the neibhbors after, apart from the old tip which was a nullptr and is now connected
		-) Perimiter elements correctly configured by ModelSubdomain
		-) Nodes before split are the same as after (apart from the old tip node) ( including middle elements)

*//*
template<uint32_t dim>
void Fracture<dim>::TestKGDMeshPropagation2D( Boundary<dim>& b_path, bool propagate_anyway){

  std::vector<InterFace<dim>*> Interfaces_before = sb_ref_.CellVector();

  if (dim == 2){
    assert( sb_ref_.PerimeterCells() == 2) ;
  }

  this->TestSplitNodeAssignment(true);

  //creating map of ordered interfaces before propagation
  std::map<Point<dim>, InterFace<dim>> ordered_if;
  for (typename std::vector<InterFace<dim>*>::iterator ifit = Interfaces_before.begin(); ifit != Interfaces_before.end(); ++ifit){
    ordered_if.insert( std::make_pair( (*ifit)->BaryCenter(), *(*ifit)) );
  }


  //finding tip node, and local node number
  size_t n_id = 0;
  Node<dim>* tipNode = nullptr ;
  InterFace<dim> tip_if = std::prev(ordered_if.end())->second;
  for (size_t n = 0; n < tip_if.FE()->Nodes(); ++n){
    if ( tip_if.N(n,INSIDE) == tip_if.N(n,OUTSIDE)){
      n_id =  n;
      tipNode = tip_if.N(n,INSIDE);
      break;
    }
  }
  if(tipNode == nullptr)
    throw csmp::Exception(ERROR, "Fracture::TestKGDMeshPropagation2D()",
                          "No tip interface found... (test meant for fracture tips located on right most part of fracture");


  assert(tip_if.ConnectedNeighbors() != tip_if.Neighbors());
  //checking left boundary correct
  assert(ordered_if.begin()->second.ConnectedNeighbors() != ordered_if.begin()->second.Neighbors());



  //WE PROPAGATE MESH-------------------------
  //------------------------------------------
  double dt_factor;
  this->PropagationAlgorithm(b_path, dt_factor, propagate_anyway);
  //------------------------------------------



  std::vector<InterFace<dim>*> Interfaces_after = sb_ref_.CellVector();
  assert( sb_ref_.PerimeterCells() == 2 );

  //Getting ordered IF post splitting
  std::map<Point<dim>, InterFace<dim>*> ordered_post;
  for (typename std::vector<InterFace<dim>*>::iterator ifit = Interfaces_after.begin(); ifit != Interfaces_after.end(); ++ifit){
    ordered_post.insert( std::make_pair( (*ifit)->BaryCenter(), (*ifit)) );
  }

  //checking sides
  assert(ordered_post.size() == 1 + ordered_if.size());

  typename std::map<Point<dim>,InterFace<dim>*>::iterator ifit_p = ordered_post.begin();
  for (typename std::map<Point<dim>,InterFace<dim>>::iterator ifit = ordered_if.begin();
       ifit != ordered_if.end(); ++ifit, ++ifit_p){

    //check that interfaces are the same!
    assert( ifit->second == *ifit_p->second );

    //check elements are the same
    assert( ifit->second.Parent(INSIDE)   == ifit_p->second->Parent(INSIDE));
    assert( ifit->second.Parent(OUTSIDE)  == ifit_p->second->Parent(OUTSIDE));
    assert( ifit->second.Parent(MIDDLE)   == ifit_p->second->Parent(MIDDLE));

    //check neighbors have stayed the same
    if ( !(ifit->second == tip_if) ){
      for (size_t ifs = 0; ifs < ifit->second.Neighbors(); ++ifs ){
        assert(ifit->second.Neighbor(ifs) == ifit_p->second->Neighbor(ifs));
      }
    } else {
      assert( ifit->second.Neighbor(n_id) == nullptr);          //check interface was originally correct
      assert( ifit_p->second->Neighbor(n_id) != nullptr);       //check same interface now has a neighbor
    }

    //check nodes are the same or duplicated
    for (size_t n = 0; n < ifit->second.FE()->Nodes(); ++n){
      if (ifit->second.N(n,INSIDE) != tipNode ){ //if not tip node
        assert( ifit->second.N(n,INSIDE)  == ifit_p->second->N(n,INSIDE) );
        assert( ifit->second.N(n,OUTSIDE) == ifit_p->second->N(n,OUTSIDE) );
        assert( ifit->second.N(n,MIDDLE)  == ifit_p->second->N(n,MIDDLE) );
      } else {
        assert( ifit_p->second->N(n,INSIDE) != ifit_p->second->N(n,OUTSIDE) );
        assert( ifit_p->second->N(n,INSIDE) == ifit->second.N(n,INSIDE));
        assert( ifit_p->second->N(n,INSIDE) == tipNode);
      }
    }

  }// end of ordered post interface loop




}

*/






/** Adapted from SplitboundaryInterface_Test::TestSplitNodeAssignment
 Test which checks the following of a splitboundary
  1) Nodes Outside != nodes Inside (for non perimiter elements)
  2) Consistent unit normals (inside and outside are opposite and inside middle the same
  3) Coordinates match on Inside, Outside and Middle
  4) Checks Boundary flags of perimeter elements & that boundary nodes are split
  5) Checks Inner Elements have no nodes on outside and vica versa
 */
/*
template<uint32_t dim>
void Fracture<dim>::TestSplitNodeAssignment( bool check_middle_elements )
{



  // loop over SplitBoundaries
  for ( typename Model<dim>::splitBoundaryIterator spbit = model_->SplitBoundariesBegin(); spbit != model_->SplitBoundariesEnd(); ++spbit )
  {

    std::set<Element<dim>*> outerParents;
    std::set<Element<dim>*> innerParents;
    std::set<Element<dim>*> middleParents;

    std::set<size_t>     outerParentsNodes;
    std::set<size_t>     innerParentsNodes;
    std::set<size_t>     middleParentNodes;
    std::set<Node<dim>*> not_duplicated_nodes;
    std::set<InterFace<dim>*> Perimeter_InterFaces;

    // loop over InterFaces
    for ( typename std::vector<InterFace<dim>* >::const_iterator ifit = (*spbit).second.ElementsBegin(); ifit != spbit->second.ElementsEnd(); ++ifit )
    {


      //0) Check that unit normals are consistently pointing
      //std::cout << "Inner Parent ID: " << (*ifit)->InnerParent()->Idx() << " and FE object id " << (*ifit)->InnerParent()->FE()->CurrentID() << std::endl;
      //Getting unit normal and coord matrix
      VectorVariable<dim> nrml_in, nrml_out, nrml_mid;
      (*ifit)->UnitNormal( nrml_in, INSIDE);

      //Comparing coordinate matrixes of element and of finiteElement
      DenseMatrix<DM_MIN> XY = (*ifit)->InnerParent()->FE()->XY;
      (*ifit)->InnerParent()->NodeCoordinateMatrix(XY);
      assert( XY(0,0) == (*ifit)->InnerParent()->FE()->XYZ(0,0) and XY(0,1) == (*ifit)->InnerParent()->FE()->XYZ(0,1) );
      assert( XY(1,0) == (*ifit)->InnerParent()->FE()->XYZ(1,0) and XY(1,1) == (*ifit)->InnerParent()->FE()->XYZ(1,1) );
      assert( XY(2,0) == (*ifit)->InnerParent()->FE()->XYZ(2,0) and XY(2,1) == (*ifit)->InnerParent()->FE()->XYZ(2,1) );

      //repeating for second unit normal
      (*ifit)->UnitNormal( nrml_out, OUTSIDE);
      //std::cout << "Fracture: Checking unit normal\t" << nrml_in.DotProduct(nrml_out) << std::endl;
      //Actual check of unit normal
      assert(nrml_in.DotProduct(nrml_out) < 0.0);

      if ( (*ifit)->ConnectedNeighbors() == (*ifit)->Neighbors() )  //if not on perimeter
      {
        //making list of in out parents
        innerParents.insert( (*ifit)->InnerParent());
        outerParents.insert( (*ifit)->OuterParent());
        middleParents.insert( (*ifit)->InterveningElement());

        //0) Check that unit normals are consistently pointing
        VectorVariable<dim> nrml_in, nrml_out, nrml_mid;
        (*ifit)->UnitNormal( nrml_in, INSIDE);
        (*ifit)->UnitNormal( nrml_out, OUTSIDE);
        assert(nrml_in.DotProduct(nrml_out) < 0.0);
        if (check_middle_elements){
          //0.1) Unit normals consistency check
          (*ifit)->UnitNormal( nrml_mid, MIDDLE);
          assert(nrml_mid.DotProduct(nrml_in) > 0.0);
          assert(nrml_mid.DotProduct(nrml_out) < 0.0);
        }

        //Check Interface nodes are indeed nodes stored on higher dim elements
        for (size_t i(0); i < (*ifit)->FE()->Nodes(); ++i){

          // 1) Check nodes are different on either side
          assert( (*ifit)->N( i, INSIDE )->Idx() != (*ifit)->N( i, OUTSIDE )->Idx() );
          assert( (*ifit)->N( i, INSIDE ) != (*ifit)->N( i, OUTSIDE ) );

          // 2) Check nodes match on both side (coordinate wise) //if nodes have been displaced this would break!
          assert( std::fabs( ((*ifit)->N( i, INSIDE )->Coordinate() -  (*ifit)->N( i, OUTSIDE )->Coordinate()).Length() ) < 0.0001 );

          // 3) Check that node is present in corresponding parent element
          assert( (*ifit)->N(i,INSIDE) == (*ifit)->InnerParent()->N( (*ifit)->ParentNodeNumber(i, INSIDE) ));
          assert( (*ifit)->N(i,OUTSIDE) == (*ifit)->OuterParent()->N( (*ifit)->ParentNodeNumber(i, OUTSIDE) ));

          // 4) Check that all nodes are on perimeter or Model
          assert( model_->Region("Model").IsPerimeterNode( (*ifit)->N(i,INSIDE) ) == true);
          assert( model_->Region("Model").IsPerimeterNode( (*ifit)->N(i,OUTSIDE) ) == true);

          if (check_middle_elements){
            //1) Assignment
            assert( (*ifit)->N(i,INSIDE)->Idx()   != (*ifit)->N(i,MIDDLE)->Idx()  );
            assert( (*ifit)->N(i,OUTSIDE)->Idx()  != (*ifit)->N(i,MIDDLE)->Idx()  );
            assert( (*ifit)->N(i,INSIDE)  != (*ifit)->N(i,MIDDLE) );
            assert( (*ifit)->N(i,OUTSIDE)  != (*ifit)->N(i,MIDDLE) );

            //2) Coordinates Match
            assert( std::fabs( ((*ifit)->N(i,MIDDLE)->Coordinate() - (*ifit)->N(i,INSIDE)->Coordinate() ).Length() ) < 0.0001);

            //making list of middle parent nodes
            middleParentNodes.insert((*ifit)->N(i,MIDDLE)->Idx());

            //4) Checking perimter of middle elements - must be false unless at boundary
            //assert( model_->Region("Model").IsPerimeterNode( (*ifit)->N(i,MIDDLE) ) == false);
          }

          //making list of in out nodes for later test
          innerParentsNodes.insert( (*ifit)->N(i,INSIDE)->Idx() );
          outerParentsNodes.insert( (*ifit)->N(i,OUTSIDE)->Idx());
        }
     } else {
        Perimeter_InterFaces.insert( *ifit );
        // Identifying
        bool atboundary = true;
        //seeing if nodes are split at perimeter
        for ( size_t i = 0; i<(*ifit)->FE()->Nodes(); ++i ){

          // Check that all nodes are on perimeter or Model
          assert( model_->Region("Model").IsPerimeterNode( (*ifit)->N(i,INSIDE) ) == true);
          assert( model_->Region("Model").IsPerimeterNode( (*ifit)->N(i,OUTSIDE) ) == true);


          if ( (*ifit)->N( i, INSIDE ) == (*ifit)->N( i, OUTSIDE ) ){
            atboundary = false;
            //check nodes are the same
            assert( (*ifit)->N( i, INSIDE ) == (*ifit)->N( i, OUTSIDE ) ); //test also pointers are the same (not just idx)
            //add to list of not duplicated nodes
            not_duplicated_nodes.insert( (*ifit)->N( i, INSIDE ) );
            //check middle nodes still new
            if (check_middle_elements){
              assert( (*ifit)->N(i,MIDDLE)->Idx() != (*ifit)->N(i,INSIDE)->Idx());  //the middle node should still be duplicate
              assert( (*ifit)->N(i,MIDDLE) != (*ifit)->N(i,INSIDE) );

              //check middle node is at perimiter
              assert( model_->Region("Model").IsPerimeterNode( (*ifit)->N(i,MIDDLE) ) == true);
            }
          }
        }
        if (atboundary){ // check that we are indeed at boudnary
          bool found_bound = false;
          for ( size_t i=0; i < (*ifit)->FE()->Nodes(); ++i){
            //if on boundary
            if ( (*ifit)->N(i,INSIDE)->AtBoundary() != INTERNAL and (*ifit)->N(i,INSIDE)->AtBoundary() != NOT){
              found_bound = true;
              //checkong nodes in and out are split at boundary
              assert( (*ifit)->N(i,INSIDE)                != (*ifit)->N(i,OUTSIDE));
              assert( (*ifit)->N(i,INSIDE)->Idx()         != (*ifit)->N(i,OUTSIDE)->Idx());
              assert( (*ifit)->N(i,OUTSIDE)->AtBoundary() == (*ifit)->N(i,INSIDE)->AtBoundary()) ;

              //checking nodes still on perimiter
              assert( model_->Region("Model").IsPerimeterNode( (*ifit)->N(i,INSIDE) ) == true);
              assert( model_->Region("Model").IsPerimeterNode( (*ifit)->N(i,OUTSIDE) ) == true);


              //checking middle element nodes
              if (check_middle_elements){
                assert( (*ifit)->N(i,MIDDLE)->AtBoundary() == (*ifit)->N(i,INSIDE)->AtBoundary());

                //check middle node is original
                assert( (*ifit)->N(i,MIDDLE) != (*ifit)->N(i,INSIDE) );
                assert( (*ifit)->N(i,MIDDLE) != (*ifit)->N(i,OUTSIDE) );
                assert( (*ifit)->N(i,MIDDLE)->Idx() != (*ifit)->N(i,INSIDE)->Idx() );
                assert( (*ifit)->N(i,MIDDLE)->Idx() != (*ifit)->N(i,OUTSIDE)->Idx() );

                //check middle element parent is at boundary
                assert( (*ifit)->N(i,MIDDLE)->Parent(0)->AtBoundary() == (*ifit)->N(i,MIDDLE)->AtBoundary());

                //check middle region node at boundary is on perimeter
                assert( model_->Region("Model").IsPerimeterNode( (*ifit)->N(i,MIDDLE) ) == true);

              }

            } // end of if found boundary node
          }// end of nodes loop

          assert(found_bound == true);
        }
      } //end of non-perimiter if

    } // end of interface loop



    //Global check of all nodes on either side - no OUTSIDE are the same as INSIDE
    //(similar to check 1, but does not rely on access operator N(i,side) )
    for ( std::set<size_t>::const_iterator it( outerParentsNodes.begin() ); it != outerParentsNodes.end(); ++it )
          assert( innerParentsNodes.find( *it ) == innerParentsNodes.end() );

    //check no nodes were found in any of the nodes in the outside elements
    if (check_middle_elements){
      for (std::set<size_t>::const_iterator it( outerParentsNodes.begin() ); it != outerParentsNodes.end(); ++it ){
        assert( middleParentNodes.find( *it ) == middleParentNodes.end());
      }
      for (std::set<size_t>::const_iterator it( innerParentsNodes.begin() ); it != innerParentsNodes.end(); ++it ){
        assert( middleParentNodes.find(*it) == middleParentNodes.end());
      }

      assert( middleParents.size() == innerParents.size());
    }


    ///Checking Element parent correct calibration
    //trivial check of same amount on both sides
    assert( innerParents.size() == outerParents.size());

    //Similar global checks
    for (typename std::set<Element<dim>*>::const_iterator it( outerParents.begin()); it != outerParents.end(); ++it){
      assert( innerParents.find( *it) == innerParents.end());                         //test no parents on both sides
      for ( size_t n = 0; n < (*it)->Nodes(); ++n){
         assert( innerParentsNodes.find( (*it)->N(n)->Idx() ) == innerParentsNodes.end() ) ; //test no nodes on inside are found on outer parent
      }
    }

    for (typename std::set<Element<dim>*>::const_iterator it(innerParents.begin()); it != innerParents.end(); ++it){
         for ( size_t n = 0; n < (*it)->Nodes(); ++n){
            assert( outerParentsNodes.find( (*it)->N(n)->Idx() ) == outerParentsNodes.end() ) ; //test no nodes on inside are found on outer parent
         }
    }


    //Checking Perimeter Nodes are so in model subdomain
    assert( Perimeter_InterFaces.size() == sb_ref_.PerimeterElements());
    for( typename std::set<Node<dim>*>::const_iterator ndit = not_duplicated_nodes.begin(); ndit != not_duplicated_nodes.end(); ndit++ ){
      assert( spbit->second.IsPerimeterNode( *ndit ) ) ;
      assert( (*ndit)->AtBoundary() == INTERNAL ); //non split node should only be internal node (cannot be bound node, that should be split by convention)
    }
    for (typename std::set<InterFace<dim>*>::const_iterator ifit = Perimeter_InterFaces.begin(); ifit != Perimeter_InterFaces.end(); ++ifit ){
      assert( spbit->second.IsPerimeterElement( *ifit ) );     //check interface is on subdomain perimeter
      for (size_t i = 0; i < (*ifit)->FE()->Nodes() ; ++i){
        //if the node is on a box boundary
        if ( (*ifit)->N(i,INSIDE)->AtBoundary() != NOT and (*ifit)->N(i,INSIDE)->AtBoundary() != INTERNAL  ){
          assert( (*ifit)->N(i,INSIDE) != (*ifit)->N(i,OUTSIDE) );               //check that it is split
          assert( (*ifit)->N(i,INSIDE)->Idx() != (*ifit)->N(i,OUTSIDE)->Idx() );
        }
      }
    }

    //Testing parent element node indexex
    for (auto& node : (*spbit).second.NodeVector() ){
      for (size_t i = 0 ; i < node->Parents(); ++i){
        //testing parent_node_indexes and parent_element_pointers is configured correctly!
        assert(node == node->Parent(i)->N( node->ParentNodeNumber(i)));
      }
    }

  } // end of splitboundaries loop

} // TestSplitNodeAssignment

*/

/*
//How can we test all initialisations?
//Code gives memory error at destruction
template<uint32_t dim>
void Fracture<dim>::TestPropertyInitialization(){

  std::list<Index> indexes;
  PropertyDatabase<dim> pref = model_->Database();
  pref.ListKeys(indexes);

  double vmin = 0 ,vmax = 0;
  for (std::list<Index>::iterator k_it = indexes.begin(); k_it != indexes.end(); ++k_it ){
    pref.RangeOf( pref.Name(*k_it), vmin,vmax);
    std::cout << "\nChecking Range of: " << pref.Name(*k_it) ;
    if ( k_it->place != ELEMENT and k_it->place != NODE){
      std::cout << "\n\tplace not on element or node" << std::endl;
      continue;
    }
    for (typename std::vector<InterFace<dim>*>::iterator ifit = sb_ref_.ElementsBegin(); ifit != sb_ref_.ElementsEnd(); ++ifit){
        if (k_it->place == ELEMENT){
          ScalarVariable sc;
          //checking values
          (*ifit)->InnerParent()->IsWithinRange(*k_it, vmin, vmax );
          (*ifit)->OuterParent()->IsWithinRange(*k_it, vmin, vmax );
          (*ifit)->InterveningElement()->IsWithinRange(*k_it, vmin, vmax);
          std::cout << "\n\tElement Value : " ;
          std::cout << (*ifit)->InnerParent()->Read(*k_it)  << ", " ;
          std::cout << (*ifit)->OuterParent()->Read(*k_it)  << ", " ;
          std::cout << (*ifit)->InterveningElement()->Read(*k_it) ;
          //checking flags
          std::cout << "\tElement Status: " ;
          std::cout << parseStatus( (*ifit)->InnerParent()->Status(*k_it) ) << ", " ;
          std::cout << parseStatus( (*ifit)->OuterParent()->Status(*k_it) ) << ", " ;
          std::cout << parseStatus( (*ifit)->InterveningElement()->Status(*k_it) ) ;

        } else if (k_it->place == NODE){
          for (size_t n = 0; n < (*ifit)->FE()->Nodes(); n++){
            (*ifit)->N(n,INSIDE)->IsWithinRange(*k_it, vmin, vmax );
            (*ifit)->N(n,OUTSIDE)->IsWithinRange(*k_it, vmin, vmax );
            (*ifit)->N(n,MIDDLE)->IsWithinRange(*k_it, vmin, vmax );

            if (k_it->type == SCALAR){
              std::cout << "\tNode Value : " ;
              std::cout << (*ifit)->N(n,INSIDE)->Read(*k_it)  << ", " ;
              std::cout << (*ifit)->N(n,OUTSIDE)->Read(*k_it)  << ", " ;
              std::cout << (*ifit)->N(n,MIDDLE)->Read(*k_it)  << std::endl;
              //checking flags
              std::cout << "\tNode Status :" ;
              std::cout << parseStatus( (*ifit)->N(n,INSIDE)->Status(*k_it) ) << ", " ;
              std::cout << parseStatus( (*ifit)->N(n,OUTSIDE)->Status(*k_it) ) << ", " ;
              std::cout << parseStatus( (*ifit)->N(n,MIDDLE)->Status(*k_it) ) << std::endl;
            } else if (k_it->type == VECTOR){
                VectorVariable<dim> vc_i, vc_o,vc_m;
                (*ifit)->N(n,INSIDE)->Read(*k_it, vc_i);
                (*ifit)->N(n,OUTSIDE)->Read(*k_it, vc_o);
                (*ifit)->N(n,MIDDLE)->Read(*k_it, vc_m);
                std::cout << "\n\tNode Value I: " ;
                for (size_t d = 0; d < dim;  ++d){
                  std::cout << vc_i(d) << ", " ;
                }
                std::cout << "\tNode Status I: " ;
                for (size_t d = 0; d < dim;  ++d){
                  std::cout << parseStatus(vc_i.Flag(d)) << ", " ;
                }
                std::cout << "\n\tNode Value O: " ;
                for (size_t d = 0; d < dim;  ++d){
                  std::cout << vc_o(d) << ", " ;
                }
                std::cout << "\tNode Status O: " ;
                for (size_t d = 0; d < dim;  ++d){
                  std::cout << parseStatus(vc_o.Flag(d)) << ", " ;
                }
                std::cout << "\n\tNode Value M: " ;
                for (size_t d = 0; d < dim;  ++d){
                  std::cout << vc_m(d) << ", " ;
                }
                std::cout << "\tNode Status M: " ;
                for (size_t d = 0; d < dim;  ++d){
                  std::cout << parseStatus(vc_m.Flag(d)) << ", " ;
                }

              }
          }// end of node loop

        } else {
            std::cout << "\n\tNot on element or node" << std::endl;
          }

    } // end index iterator
  } //end interface loop

  indexes.clear();

}
*/
template class Fracture<1U>;
template class Fracture<2U>;
template class Fracture<3U>;



} //csmp
