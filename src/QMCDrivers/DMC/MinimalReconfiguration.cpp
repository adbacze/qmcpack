//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Andrew D. Baczewski, adbacze@sandia.gov, Sandia National Laboratories
//
// File created by: Andrew D. Baczewski, adbacze@sandia.gov, Sandia National Laboratories               
//////////////////////////////////////////////////////////////////////////////////////
    


#include "QMCDrivers/DMC/MinimalReconfiguration.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

/** default constructor
 *
 * set SwapMode
 */
MinimalReconfiguration::MinimalReconfiguration(Communicate* c) :WalkerControlBase(c)
{
  SwapMode=1;
}

int MinimalReconfiguration::reconfigureWalkers(MCWalkerConfiguration& W)
{
  int nw(W.getActiveWalkers()); //get the number of active walkers and store it in nw
  if(wConfTilde.size()!=nw)    //check to see if instance has set size of scaled weight vector appropriately 
  {
    wConfTilde.resize(nw);     //if not (because 1st call from instance) then resize appropriately
  }
  //process walker contents 
  RealType esum=0.0,e2sum=0.0,wtot=0.0,ecum=0.0;
  MCWalkerConfiguration::iterator it(W.begin());
  RealType r2_accepted=0.0,r2_proposed=0.0;
  for(int iw=0; iw<nw; iw++)
  {
    r2_accepted+=(*it)->Properties(R2ACCEPTED);
    r2_proposed+=(*it)->Properties(R2PROPOSED);
    RealType wgt((*it)->Weight);
    RealType e((*it)->Properties(LOCALENERGY));
    esum += wgt*e;
    e2sum += wgt*e*e;
    ecum += e;
    wtot += wConfTilde[iw]=wgt; //note that we will be scaling each entry of wConfTilde  later
    ++it;
  }
  curData[ENERGY_INDEX]=esum;
  curData[ENERGY_SQ_INDEX]=e2sum;
  curData[WALKERSIZE_INDEX]=nw;
  curData[WEIGHT_INDEX]=wtot;
  curData[EREF_INDEX]=ecum;
  curData[R2ACCEPTED_INDEX]=r2_accepted;
  curData[R2PROPOSED_INDEX]=r2_proposed;

  //divide the sum of weights by the number of walkers to get the global weight in Eqn. 49 of ACK 2000
  wGlobal=wtot/static_cast<RealType>(nw);
  std::vector<int> plus,minus;
  RealType NreconfPlus=0.0, NreconfMinus=0.0;
  for(int iw=0; iw<nw; iw++)
  {
    wConfTilde[iw]=wConfTilde[iw]/wGlobal;
    if(wConfTilde[iw]>1.0){
       plus.push_back(iw);
       NreconfPlus += std::abs(wConfTilde[iw]-1.0);
    }
    else{
       minus.push_back(iw);	   
       NreconfMinus += std::abs(wConfTilde[iw]-1.0); 
    }
  }
 
  app_log() << "   NreconfPlus: " << NreconfPlus << std::endl;
  Nreconf = static_cast<int>(std::round(NreconfPlus + Random()));

  curData[FNSIZE_INDEX]=nw-Nreconf;
  curData[RNONESIZE_INDEX]=Nreconf;

  int copyIndex, killIndex, copyWalker, killWalker;
  if(Nreconf){
     for(int r=0; r<Nreconf; r++){
	
	// check to see if there is a preferred way to generate random indices
	copyIndex = static_cast<int>(Random()*plus.size());
	killIndex = static_cast<int>(Random()*minus.size());
	copyWalker = plus[copyIndex];
	killWalker = minus[killIndex];

        W[killWalker]->makeCopy(*(W[copyWalker]));
        W[killWalker]->ParentID=W[copyWalker]->ID;
        W[killWalker]->ID=(++NumWalkersCreated)*NumContexts+MyContext;

	// remove walkers that were copied/overwritten so we don't try twice
	plus.erase(plus.begin()+copyIndex);
	minus.erase(minus.begin()+killIndex);
     }
  }

  return nw;
}

int MinimalReconfiguration::branch(int iter, MCWalkerConfiguration& W, RealType trigger)
{
  int nwkept = reconfigureWalkers(W);
  app_log() << "   # of reconfigurations in current generation: " << Nreconf << std::endl;
  //update EnsembleProperty
  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end)
  {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }
  return nwkept;
}


