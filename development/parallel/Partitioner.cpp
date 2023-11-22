#include "Partitioner.h"

extern "C" {
//#include "metis.h"
}

extern "C" void METIS_PartMeshDual(int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern "C" void METIS_PartMeshNodal(int *, int *, int *, int *, int *, int *, int *, int *, int *);


using namespace std;

namespace csp {

extern csp::ErrorHandler           skm_err;

template<typename fT, stl_index dim>
Partitioner<fT,dim>::Partitioner( const VSet<fT,dim>& vs, int n )
 : vset(vs),
   n_parts(n)
 {
 }


// Write a new function which uses the node-partition as given by METIS,
// create node_vectors for each partition
// loop over elements of vset, loop over nodes of element
// element belongs to partition to which node belongs
// use a set to store elements per partition (no multiple entries)
template<typename fT, stl_index dim>
void Partitioner<fT,dim>::FastPartition( std::vector<std::vector<stl_index> >& element_vector, bool three_D )
 {
    cout <<"\nPartitioner::FastPartition: Start partitioning for "<<n_parts<<" processors."<<endl;

    // Initialise arguments for METIS function
    ne = vset.Elements();
    nn = vset.Vertices();
    int etype = 1;
    if(three_D) etype = 2;
    int numflag = 1;
    int edgecut = 0;

    // Write plist to single vector
    std::vector<int> elmnts;
    elmnts.resize(0);
    epart.resize(ne);
    npart.resize(nn);

    for (std::deque<std::vector<stl_index> >::const_iterator it=vset.PlistBegin(); it!=vset.PlistEnd(); it++)
      for (unsigned int i=0;i!=(*it).size();i++)
	elmnts.push_back((*it)[i]);

    for (unsigned int i=0;i!=100;i++)
      cout <<elmnts[i]<< endl;

    cout <<"\nPartitioner<fT,dim>::FastPartition: Calling METIS library." << endl;

//    METIS_PartMeshDual(&ne,&nn,&elmnts[0],&etype,&numflag,&n_parts,&edgecut,&epart[0],&npart[0]);
    METIS_PartMeshNodal(&ne,&nn,&elmnts[0],&etype,&numflag,&n_parts,&edgecut,&epart[0],&npart[0]);

    // create vector (sized n_parts) continaing empty sets
    std::vector<std::set<stl_index> > element_set;
    element_set.resize(n_parts);
    for(unsigned int i=0;i!=n_parts;i++)
      element_set[i].clear();

    // Fill each set with element ID's belonging to that partition
    for(stl_index eid=1;eid!=(vset.Elements()+1U);eid++)
      for(std::vector<stl_index>::const_iterator nit=vset.PlistBegin(eid); nit!=vset.PlistEnd(eid); nit++)
        element_set[npart[(*nit)-1]-1].insert(eid);

    // create empty vectors
    element_vector.resize(n_parts);
    for(unsigned int i=0;i!=n_parts;i++)
      element_vector[i].resize(0);

    cout <<"\nPartitioner<fT,dim>::FastPartition: Creating overlapping elements." << endl;

    // create element vector from element set
    for(unsigned int i=0;i!=n_parts;i++)
      for(std::set<stl_index>::const_iterator it=element_set[i].begin(); it!=element_set[i].end(); it++)
        element_vector[i].push_back((*it));

    // create vector (sized n_parts) continaing empty sets
    std::vector<std::set<stl_index> > node_set;
    node_set.resize(n_parts);
    for(unsigned int i=0;i!=n_parts;i++)
      node_set[i].clear();

    // Fill each set with node ID's belonging to that partition
    for(unsigned int i=0;i!=n_parts;i++)
      for(std::set<stl_index>::const_iterator it=element_set[i].begin(); it!=element_set[i].end(); it++)
        for(std::vector<stl_index>::const_iterator nit=vset.PlistBegin(*it); nit!=vset.PlistEnd(*it); nit++)
          node_set[i].insert(*nit);

    cout <<"\nPartitioner<fT,dim>::FastPartition: Determine outerhalo node ID's." << endl;

    // clear outerhalo nodes set
    outerhalo.resize(n_parts);
    for(unsigned int i=0;i!=n_parts;i++)
      outerhalo[i].clear();

    // Fill outerhalo node set
    for(unsigned int i=0;i!=n_parts;i++)
      for(std::set<stl_index>::const_iterator nit=node_set[i].begin(); nit!=node_set[i].end(); nit++)
        if((npart[(*nit)-1]-1)!=i)
          // this node is outerhalo node of partition i
          // it belongs to partition npart[(*nit)-1]-1
          outerhalo[i].insert(make_pair((*nit),(npart[(*nit)-1]-1)));

    cout <<"\nPartitioner<fT,dim>::FastPartition: Successfully created "<<n_parts<<" partitions" << endl;

 }





template<typename fT, stl_index dim>
void Partitioner<fT,dim>::MetisPartition( std::vector<std::vector<stl_index> >& element_vector, bool element_partition )
 {
 
    // Initialise arguments for METIS function
    ne = vset.Elements();
    nn = vset.Vertices();
    int n_(n_parts);
    int etype = 1;
    int numflag = 1;
    int edgecut = 0;    

    // Write plist to single vector
    std::vector<int> elmnts;
    element_vector.resize(n_parts);
    std::deque<std::vector<stl_index> >::const_iterator it=vset.PlistBegin();
    elmnts.resize((*it).size()*ne);
    epart.resize(ne);
    npart.resize(nn);
    
    unsigned int k(0);
    for (std::deque<std::vector<stl_index> >::const_iterator it=vset.PlistBegin(); it!=vset.PlistEnd(); it++)
      for (unsigned int i=0;i!=(*it).size();i++)
        {
	      elmnts[k]=(*it)[i];
	      k++;
	}

    cout <<"\nPartitioner::MetisPartition: Start partitioning for "<<n_parts<<" processors."<<endl;

    cout <<"\nPartitioner<fT,dim>::MetisPartition: Calling METIS library." << endl;

    if (element_partition)
      METIS_PartMeshDual(&ne,&nn,&elmnts[0],&etype,&numflag,&n_,&edgecut,&epart[0],&npart[0]);
    else
      METIS_PartMeshNodal(&ne,&nn,&elmnts[0],&etype,&numflag,&n_,&edgecut,&epart[0],&npart[0]);

    cout <<"\n\tMetis Returned, total edgecut = "<<edgecut<<endl;
    
    // Set all sizes in vectors of element_vector to 0
    for (unsigned int i=0; i!=element_vector.size(); i++)
      element_vector[i].resize(0);

    unsigned int part(1);

    // Fill element_vector with elements numbers
    for (std::vector<std::vector<stl_index> >::iterator it=element_vector.begin(); it!=element_vector.end(); it++, part++)
      for (unsigned int j=0; j!=epart.size(); j++)
        // If this element is part of partition i
 	    if(epart[j]==part)
 	      (*it).push_back(j+1);

    CreateNodeVector(element_vector);

    SharedNodes();
    
    cout <<"\nPartitioner<fT,dim>::MetisPartition: Creating overlapping elements." << endl;

    CreateHaloElements(element_vector);  // do only once!!
    AddHaloElements(element_vector);
    CreateHaloNodes();
    ErroneousElements(element_vector);
    RemoveDoubleEntries(element_vector);

    cout <<"\nPartitioner<fT,dim>::MetisPartition: Dsetermine outerhalo node ID's." << endl;

    RemoveDoubleOuterhaloEntries();
    CreateNodeVector(element_vector);

    for(unsigned int i=0;i!=element_vector.size();i++) assert (!element_vector[i].empty());

    cout <<"\nPartitioner<fT,dim>::MetisPartition: Successfully created "<<n_parts<<" partitions" << endl;
    PartitionInfoToScreen(element_vector);
 }

template<typename fT, stl_index dim>
void Partitioner<fT,dim>::PartitionInfoToScreen(std::vector<std::vector<stl_index> >& element_vector)
 {
    int part=0;
    for (std::vector<std::vector<stl_index> >::iterator it=element_vector.begin(); it!=element_vector.end(); it++, part++)
     {
      cout<<"\nPartitioner<fT,dim>::PartitionInfoToScreen: Partition "<<part<<":\tTotal number of elements = "<<(*it).size()<< endl;
      cout<<"\t\t\t\t\tNumber of halo elements   = "<<halo_elements[part].size()<< endl;
      cout<<"\t\t\t\t\tTotal number of nodes     = "<<node_vector[part].size()<< endl;
      cout<<"\t\t\t\t\tNumber of innerhalo nodes = "<<innerhalo[part].size()<< endl;
      cout<<"\t\t\t\t\tNumber of outerhalo nodes = "<<outerhalo[part].size()<< endl;
     }
 }

template<typename fT, stl_index dim>
void Partitioner<fT,dim>::CreateNodeVector(std::vector<std::vector<stl_index> >& element_vector)
 {
   // THIS FUNCTION MUST BE CALLED WHEN THE HALO-ELEMENTS ARE NOT YET ADDED TO ELEMENT-VECTOR
   
   // resize vector, clear set
   node_vector.resize(n_parts);
   for (unsigned int p=0;p!=n_parts;p++) node_vector[p].clear();

   // loop over partitions
   for (unsigned int p=0;p!=n_parts;p++)
     // Loop over elements of partition
     for (std::vector<stl_index>::const_iterator it=element_vector[p].begin();it!=element_vector[p].end();it++)
       // loop over nodes element
       for (std::vector<stl_index>::iterator i=vset.PlistBegin(*it);i!=vset.PlistEnd(*it);i++)
         // Insert node in set
	 node_vector[p].insert(*i);

 }



template<typename fT, stl_index dim>
void Partitioner<fT,dim>::ErroneousElements(std::vector<std::vector<stl_index> >& element_vector)
 {
    // Some tripple elements need to be removed
    // If all nodes of one element are in the outer halo of a partition
    // --> remove this element from halo_element vector
    
    std::vector<stl_index> temp;
    std::vector<std::set<std::pair<stl_index,int> > >  erroneous_halo_elements;
    std::vector<std::set<std::pair<stl_index,int> > >  switched_halo_elements;
    // Resize:
    erroneous_halo_elements.resize(n_parts);
    for(unsigned int p=0;p!=n_parts;p++) erroneous_halo_elements[p].clear();
    switched_halo_elements.resize(n_parts);
    for(unsigned int p=0;p!=n_parts;p++) switched_halo_elements[p].clear();

    // Remove old halo_element vector from element_vector
    RemoveHaloElements(element_vector);

    // Loop over partitions
    for(int p=0;p!=n_parts;p++)
     {
       // pump outerhalo in temp vector
       unsigned int i(0);
       temp.resize(outerhalo[p].size());
       for(std::set<std::pair<stl_index,int> >::const_iterator it=outerhalo[p].begin();it!=outerhalo[p].end();it++,i++)
         temp[i]=(*it).first;
       // Loop over halo element vector
       for (std::set<std::pair<stl_index,int> >::iterator eid=halo_elements[p].begin();eid!=halo_elements[p].end();eid++)
        {
          unsigned int number(0);
	  // Loop over nodes of halo element
	  for (vector<stl_index>::const_iterator nit=vset.PlistBegin((*eid).first);nit!=vset.PlistEnd((*eid).first);nit++)
	   {
	     // find node id in outerhalo
             std::vector<stl_index>::iterator f=find(temp.begin(),temp.end(),*nit);
	     if(f!=temp.end()) number++;
	   }
          if(number==vset.PlistSize((*eid).first)) // All nodes are outerhalo nodes
	    {
	      // This element needs to be deleted from element_vector of partition p
              std::vector<stl_index>::iterator f=find(element_vector[p].begin(),element_vector[p].end(),(*eid).first);
	      if(f!=element_vector[p].end())
	       {
		 // This means an "original" METIS element has been removed from partition p
		 // Its nodes are all outerhalo nodes
		 // This generally means that this element is at a junction of several partitions
		 // Determine halo's to which element belongs, add entries in halo vector WITH CORRECT 2nd entry

		 // loop over halo_element[p] -> it
		 // if it.first == eid.first, push back it.second (other partition) in vector
		 // if vector.empty : error
		 // if vector.size > n_parts: error
		 // some kind of loop over vector
		 // halo_elements[vector[i]].insert(*eid.first,vector[i-1])...

 	         switched_halo_elements[p].insert(*eid);
		 cout << "\nPartitioner::ErroneousElements: Switched element "<<(*eid).first<<" from partition "<<p<<" to partition "<<(*eid).second<< endl;

	       }
	      else
	       {
	      //   cout<<"ERROR: Inconsistency..."<< endl;
	      // This entry (eid) also needs to be deleted from halo_element vector of partition p (shared with any other processor)
	      // And from halo_element vector to which the element belongs, shared with partition p
                 erroneous_halo_elements[p].insert(*eid);
		 cout<<"\nPartitioner::ErroneousElements: Removed element "<<(*eid).first<<" from partition "<<p<< endl;
               }

	    }
	}
      }

    // Remove erroneous entries from halo_elements vector
    for(int p=0;p!=n_parts;p++)
      for(std::set<std::pair<stl_index,int> >::iterator it=erroneous_halo_elements[p].begin();it!=erroneous_halo_elements[p].end();it++)
        {
	  // first remove entry from halo_elements[p]
          std::set<std::pair<stl_index,int> >::iterator f=find(halo_elements[p].begin(),halo_elements[p].end(),(*it));
	  if(f!=halo_elements[p].end())
	    halo_elements[p].erase(*f);
          else
	    skm_err.notice( CSP_ERROR, "Partitioner<fT,dim>::ErroneousElements",
                                       "\nError removing erroneous element from partition.");

	  // Now delete halo_elements entry (*it).first of partition (*it).second which it shares with p
	  f=find(halo_elements[(*it).second].begin(),halo_elements[(*it).second].end(),make_pair((*it).first,p));
	  if(f!=halo_elements[(*it).second].end())
	    halo_elements[(*it).second].erase(*f);
	  else
	    skm_err.notice( CSP_ERROR, "Partitioner<fT,dim>::ErroneousElements",
                                       "\nError removing erroneous element from partition.");
	}
     
     // Add switched entries to halo elements vector
     for(int p=0;p!=n_parts;p++)
      {
       // Remove entry from halo_element vector
       for(std::set<std::pair<stl_index,int> >::iterator it=switched_halo_elements[p].begin();it!=switched_halo_elements[p].end();it++)
        {
	  // first remove entry from halo_elements[p]
          std::set<std::pair<stl_index,int> >::iterator f=find(halo_elements[p].begin(),halo_elements[p].end(),(*it));
	  if(f!=halo_elements[p].end())
	    halo_elements[p].erase(*f);
          else
	    skm_err.notice( CSP_ERROR, "Partitioner<fT,dim>::ErroneousElements",
                                       "\nError removing erroneous element from partition.");
        }
       // Remove entry from element_vector
       for(std::set<std::pair<stl_index,int> >::iterator it=switched_halo_elements[p].begin();it!=switched_halo_elements[p].end();it++)
        {
          std::vector<stl_index>::iterator f=find(element_vector[p].begin(),element_vector[p].end(),(*it).first);
	  if(f!=element_vector[p].end())
	    element_vector[p].erase(f);
	}
       // Deal with changing halo's in other partitions (!=p)
       for(std::set<std::pair<stl_index,int> >::const_iterator it=switched_halo_elements[p].begin();it!=switched_halo_elements[p].end();it++)
         for(std::set<std::pair<stl_index,int> >::const_iterator it2=switched_halo_elements[p].begin();it2!=switched_halo_elements[p].end();it2++)
           if((*it).first==(*it2).first && (*it).second!=(*it2).second)
	     {
	       // find corresponding halo elements, first in
               std::set<std::pair<stl_index,int> >::iterator f=find(halo_elements[(*it).second].begin(),halo_elements[(*it).second].end(),make_pair((*it).first,p));
	       if(f!=halo_elements[(*it).second].end())
	        {
                 //(*f).second = (*it2).second;
                 halo_elements[(*it).second].erase(f);
                 halo_elements[(*it).second].insert(make_pair((*f).first,(*it2).second));
                 cout<<"\nPartitioner::ErroneousElements: Switched halo element "<<(*f).first<<" of partition "<<(*it).second<<" \
from neighbor = "<<p<<" to neighbor = "<<(*f).second<< endl;
	        }
	       // and vice versa...
               f=find(halo_elements[(*it2).second].begin(),halo_elements[(*it2).second].end(),make_pair((*it2).first,p));
	       if(f!=halo_elements[(*it2).second].end())
	        {
	        //(*f).second = (*it).second;
                 halo_elements[(*it2).second].erase(f);
                 halo_elements[(*it2).second].insert(make_pair((*f).first,(*it).second));
                 cout<<"\nPartitioner::ErroneousElements: Switched halo element "<<(*f).first<<" of partition "<<(*it2).second<<" \
from neighbor = "<<p<<" to neighbor = "<<(*f).second<< endl;
	        }

	     }

     }



     // Add new halo_element vector
     AddHaloElements(element_vector);

     // Recreate Halo_nodes
     CreateHaloNodes();

 }

template<typename fT, stl_index dim>
void Partitioner<fT,dim>::AddHaloElements(std::vector<std::vector<stl_index> >& element_vector)
 {
     // Add Halo Elements to element_vector
    for (unsigned int p=0;p!=element_vector.size();p++)
      for (std::set<std::pair<stl_index,int> > ::const_iterator it_halo=halo_elements[p].begin();it_halo!=halo_elements[p].end();it_halo++)
        element_vector[p].push_back((*it_halo).first);
 }

template<typename fT, stl_index dim>
void Partitioner<fT,dim>::RemoveDoubleEntries(std::vector<std::vector<stl_index> >& element_vector)
 {
    // Add Halo Elements to element_vector, use temporary set such that there are NO double entries in element_vector
    std::vector<std::set<stl_index> > temp;
    //Resize vector and clear sets
    temp.resize(n_parts);
    for(unsigned int i=0;i!=temp.size();i++)temp[i].clear();

    // 1st write core elements into set
    for (unsigned int p=0;p!=n_parts;p++)
      for (std::vector<stl_index>::const_iterator it=element_vector[p].begin();it!=element_vector[p].end();it++)
        temp[p].insert(*it);
    
    // 2nd write halo elements into set
    for (unsigned int p=0;p!=n_parts;p++)
      for (std::set<std::pair<stl_index,int> > ::const_iterator it_halo=halo_elements[p].begin();it_halo!=halo_elements[p].end();it_halo++)
        temp[p].insert((*it_halo).first);

    // Clear element vector
    for (unsigned int p=0;p!=n_parts;p++)
      element_vector[p].clear();

    // Write set into element vector
    for (unsigned int p=0;p!=n_parts;p++)
      element_vector[p].resize(temp[p].size());
    for (unsigned int p=0;p!=n_parts;p++)
      copy(temp[p].begin(),temp[p].end(),element_vector[p].begin());
 }

 
template<typename fT, stl_index dim>
void Partitioner<fT,dim>::RemoveHaloElements(std::vector<std::vector<stl_index> >& element_vector)
 {
     // Remove Halo Elements from element_vector
    for (unsigned int p=0;p!=element_vector.size();p++)
      for (std::set<std::pair<stl_index,int> >::const_iterator it_halo=halo_elements[p].begin();it_halo!=halo_elements[p].end();it_halo++)
       {
         std::vector<stl_index>::iterator f=find(element_vector[p].begin(),element_vector[p].end(),(*it_halo).first);
	 if(f!=element_vector[p].end())
	   element_vector[p].erase(f);
       }

 }


template<typename fT, stl_index dim>
void Partitioner<fT,dim>::SharedNodes()
 {
    // resize vector shared_nodes
    shared_nodes.resize(n_parts);
    // vector of sets of pairs. Size vector is # of partitions.
    std::pair<stl_index,int> node_p;

    unsigned int p(1);
    // Loop over shared_nodes
    for (std::vector<std::set<pair<stl_index,int> > >::iterator it=shared_nodes.begin();it!=shared_nodes.end();it++, p++)
      // Loop over epart
      for (unsigned int el_id=1; el_id!=(epart.size()+1);el_id++)
        // Loop over nodes of this element
        for (unsigned int n_loc=0; n_loc!=(*vset.PlistBegin()).size(); n_loc++)
            if (epart[el_id-1]==p && epart[el_id-1]!=npart[vset.Plist(el_id,n_loc)-1])
              // node is not part of the same partition as element
              // store node number + partition to which it belongs in shared_nodes
              {
                node_p.first = vset.Plist(el_id,n_loc);
                node_p.second= npart[vset.Plist(el_id,n_loc)-1];
                (*it).insert(node_p);
              }

    p = 1;
    unsigned int p2;
    // Rewrite such that every set contains all nodes which bound with other partition
    for (std::vector<std::set<pair<stl_index,int> > >::iterator it=shared_nodes.begin();it!=shared_nodes.end();it++, p++)
      {
        p2=1;
        for (std::vector<std::set<pair<stl_index,int> > >::iterator it2=shared_nodes.begin();it2!=shared_nodes.end();it2++,p2++)
          if((*it)!=(*it2))
            for (std::set<pair<stl_index,int> >::const_iterator pit=(*it2).begin();pit!=(*it2).end();pit++)
              if((*pit).second==p)
                {
			      node_p.first = (*pit).first;
			      node_p.second= p2;
			      (*it).insert(node_p);
                 }
      }

    // Triple points need to be added to list as well
    p2=1;
    for (std::vector<std::set<pair<stl_index,int> > >::iterator it=shared_nodes.begin();it!=shared_nodes.end();it++,p2++)
      for (std::set<pair<stl_index,int> >::const_iterator pit=(*it).begin();pit!=(*it).end();pit++)
        {
        p = 1;
        for (std::vector<std::set<pair<stl_index,int> > >::const_iterator it2=shared_nodes.begin();it2!=shared_nodes.end();it2++, p++)
          if((*it)!=(*it2))
            for (std::set<pair<stl_index,int> >::const_iterator pit2=(*it2).begin();pit2!=(*it2).end();pit2++)
              if((*pit2).first==(*pit).first && (*pit2).second!=p2)
			    {
			      // Node is triple point, add to list.
			      node_p.first = (*pit2).first;
			      node_p.second= p;
			      (*it).insert(node_p);
			    }
        }


    // Renumber 0...n_parts-1
    for (std::vector<std::set<pair<stl_index,int> > >::iterator it=shared_nodes.begin();it!=shared_nodes.end();it++)
      for (std::set<pair<stl_index,int> >::iterator pit=(*it).begin();pit!=(*it).end();pit++)// (*pit).second--;
       {        
        (*it).erase(pit);
        (*it).insert(make_pair((*pit).first,((*pit).second - 1)));
       }

 }


template<typename fT, stl_index dim>
void Partitioner<fT,dim>::GetSharedNodes( unsigned int proc, std::vector<std::pair<int,int> >& shared_n )
 {
    //
    unsigned int p(0);
    shared_n.resize(0);
    for (std::vector<std::set<pair<stl_index,int> > >::const_iterator it=shared_nodes.begin();it!=shared_nodes.end();it++, p++)
      if(p==proc)
		for (std::set<pair<stl_index,int> >::const_iterator pit=(*it).begin();pit!=(*it).end();pit++)
		  shared_n.push_back(*pit); 	      

 }

template<typename fT, stl_index dim>
void Partitioner<fT,dim>::GetHaloNodes(std::vector<std::set<std::pair<stl_index,int> > >&  inhalo, std::vector<std::set<std::pair<stl_index,int> > >&  outhalo )
 {
    inhalo = innerhalo;
    outhalo = outerhalo;
 }

template<typename fT, stl_index dim>
void Partitioner<fT,dim>::GetOuterhaloNodes(std::vector<std::set<std::pair<stl_index,int> > >&  outhalo )
 {
    outhalo = outerhalo;
 }


template<typename fT, stl_index dim>
void Partitioner<fT,dim>::CreateHaloElements(const std::vector<std::vector<stl_index> >& element_vector)
 {

   // partitions share elements with each other these need to be added to element_vector
   // Partition 0 takes elements from Partition 1
   // Partition 1 takes elements from Partition 2
   // etc..


   std::deque<std::vector<stl_index> >::const_iterator plist_it;
   
   // 1. First create halo_element vector (vector<set<pair<stl_index,int> > >)
   
   // Resize vector
   halo_elements.resize(n_parts);
   int eid(1);

   // loop over partitions
   for (unsigned int p=0;p!=n_parts;p++)
     // Loop over elements of partition
     for (std::vector<stl_index>::const_iterator it=element_vector[p].begin();it!=element_vector[p].end();it++)
       // loop over nodes element
       for (std::vector<stl_index>::iterator i=vset.PlistBegin(*it);i!=vset.PlistEnd(*it);i++)
         // loop over shared_nodes[p]
	 for (std::set<std::pair<stl_index,int> >::const_iterator nit=shared_nodes[p].begin();nit!=shared_nodes[p].end();nit++)
	   if((*i)==(*nit).first)
	    {
	     // (*i) is a shared node on boundary
	     // between partitions p and (*nit).second

	     if (p > (*nit).second)
	       // eid belong to p and is a halo element of both p and (*nit).second
              {
	        halo_elements[p].insert(make_pair((*it),(*nit).second));
	        halo_elements[(*nit).second].insert(make_pair((*it),p));
              }
	    }
 }



template<typename fT, stl_index dim>
void Partitioner<fT,dim>::CreateHaloNodes()
 {

   // 2. Create Inner and outerhalo vectors (vector<set<pair<stl_index,int> > >)

   // Resize vector / Clear set
   innerhalo.resize(n_parts);
   outerhalo.resize(n_parts);
   for(int p=0;p!=n_parts;p++) innerhalo[p].clear();
   for(int p=0;p!=n_parts;p++) outerhalo[p].clear();

   // First create outerhalo
   // Loop over partitions
   for (int p=0;p!=n_parts;p++)
    {
      outerhalo[p].clear();
      // Loop over halo elements
      for(std::set<std::pair<stl_index,int> >::const_iterator it=halo_elements[p].begin();it!=halo_elements[p].end();it++)
        // loop over nodes it
	for (vector<stl_index>::const_iterator nit=vset.PlistBegin((*it).first);nit!=vset.PlistEnd((*it).first);nit++)
         {
	   std::set<stl_index>::const_iterator f=find(node_vector[p].begin(),node_vector[p].end(),(*nit));
	   if(f==node_vector[p].end()) // node is outerhalo node
	     outerhalo[p].insert(make_pair(*nit,(*it).second));
	   else if(p>(*it).second)
	    {
	      std::set<std::pair<stl_index,int> >::const_iterator f2=find(shared_nodes[p].begin(),shared_nodes[p].end(),make_pair(*nit,(*it).second));
              if(f2!=shared_nodes[p].end()) // node is outerhalo node
	        outerhalo[p].insert(make_pair(*nit,(*it).second));
	    }
          }
     }
   std::vector<stl_index> temp;
   // Now create innerhalo nodes
   for (int p=0;p!=n_parts;p++)
    {
      innerhalo[p].clear();
      // pump outerhalo in temp vector
      temp.resize(outerhalo[p].size());
      unsigned int i(0);
      for(std::set<std::pair<stl_index,int> >::const_iterator it=outerhalo[p].begin();it!=outerhalo[p].end();it++,i++)
        temp[i]=(*it).first;

      // Loop over halo elements
      for(std::set<std::pair<stl_index,int> >::const_iterator it=halo_elements[p].begin();it!=halo_elements[p].end();it++)
        // loop over nodes it
	for (vector<stl_index>::const_iterator nit=vset.PlistBegin((*it).first);nit!=vset.PlistEnd((*it).first);nit++)
         {
	   std::vector<stl_index>::const_iterator f=find(temp.begin(),temp.end(),(*nit));
	   if (f==temp.end()) // nit is not an outerhalo, ie it is an innerhalo
	     innerhalo[p].insert(make_pair(*nit,(*it).second));
	 }
    }

 }
 

template<typename fT, stl_index dim>
void Partitioner<fT,dim>::RemoveDoubleOuterhaloEntries()
 {
    
   std::vector<std::set<std::pair<stl_index,int> > >  double_entry;
   double_entry.resize(n_parts);
   for(int p=0;p!=n_parts;p++) double_entry[p].clear();

   // Loop over partitions
   for(int p=0;p!=n_parts;p++)
     // 2 Loops over outerhalo[p]
     for(std::set<std::pair<stl_index,int> >::iterator it=outerhalo[p].begin();it!=outerhalo[p].end();it++)
       for(std::set<std::pair<stl_index,int> >::iterator it2=outerhalo[p].begin();it2!=outerhalo[p].end();it2++)
         if((*it).first==(*it2).first && (*it).second!=(*it2).second)
	  {
	    // Double entry:
	    // Only that one should be saved which is innerhalo node of neighbour partition
	    std::set<std::pair<stl_index,int> >::const_iterator f=find(innerhalo[(*it).second].begin(),innerhalo[(*it).second].end(),make_pair((*it).first,p));
	    if(f==innerhalo[(*it).second].end())
	     {
	      // this ID belongs NOT to innerhalo of partition (*it).second
	      // and therefor needs to be removed from outerhalo[p]
              double_entry[p].insert(*it);
	     }
          }

   // Loop over partitions
   for(int p=0;p!=n_parts;p++)
    // Loop over doubles
    for(std::set<std::pair<stl_index,int> >::const_iterator it=double_entry[p].begin();it!=double_entry[p].end();it++)
     {
       std::set<std::pair<stl_index,int> >::const_iterator f=find(outerhalo[p].begin(),outerhalo[p].end(),(*it));
       if (f!=outerhalo[p].end())
         outerhalo[p].erase(*f);
       else
         skm_err.notice( CSP_ERROR, "Partitioner<fT,dim>::RemoveDoubleOuterhaloEntries",
                                    "\nError while removing erroneous element from partition.");
     }

 }

template class Partitioner<csp_float,1U>;
template class Partitioner<csp_float,2U>;
template class Partitioner<csp_float,3U>;

} // end namespace csp
