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
    


#include "QMCDrivers/DMC/MinimalReconfigurationMPI.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

/** default constructor
 *
 * set SwapMode
 */
MinimalReconfigurationMPI::MinimalReconfigurationMPI(Communicate* c) :WalkerControlBase(c)
{
  SwapMode=1;
}

int MinimalReconfigurationMPI::reconfigureWalkers(MCWalkerConfiguration& W)
{
  int nw(W.getActiveWalkers()); // get the number of active walkers and store it in nw
  if(TotalWalkers !=nw*NumContexts)   	   
  {
    FirstWalker=nw*MyContext;
    LastWalker=FirstWalker+nw;
    TotalWalkers=nw*NumContexts;
    
    wConfTilde.resize(nw);   
    wOffset.resize(NumContexts+1);
    NreconfPlus.resize(NumContexts);
    NreconfMinus.resize(NumContexts);     
    dN.resize(NumContexts+4);
  }

  //accumulate the energies
  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  int iw=0;
  RealType esum=0.0,e2sum=0.0,wtot=0.0,ecum=0.0;
  RealType r2_accepted=0.0,r2_proposed=0.0;
  while(it != it_end)
  {
    r2_accepted+=(*it)->Properties(R2ACCEPTED);
    r2_proposed+=(*it)->Properties(R2PROPOSED);
    RealType wgt((*it)->Weight);
    RealType e((*it)->Properties(LOCALENERGY));
    esum += wgt*e;
    e2sum += wgt*e*e;
    wtot += wgt;
    ecum += e;
    wConfTilde[iw++]=wgt;
    ++it;
  }
  curData[ENERGY_INDEX]=esum;
  curData[ENERGY_SQ_INDEX]=e2sum;
  curData[WALKERSIZE_INDEX]=nw;
  curData[WEIGHT_INDEX]=wtot;
  curData[EREF_INDEX]=ecum;
  curData[R2ACCEPTED_INDEX]=r2_accepted;
  curData[R2PROPOSED_INDEX]=r2_proposed;
  std::fill(curData.begin()+LE_MAX,curData.end(),0.0);
  curData[LE_MAX+MyContext]=wtot;
  //collect everything
  myComm->allreduce(curData);
  //update EnsembleProperty
  W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];
  wOffset[0]=0.0;
  for(int ip=0,jp=LE_MAX; ip<NumContexts; ip++,jp++)
    wOffset[ip+1] = wOffset[ip] + curData[jp];
  
  wGlobal=wOffset[NumContexts]/static_cast<RealType>(TotalWalkers); 

  std::vector<int> plus,minus;
  std::fill(NreconfPlus.begin(),NreconfPlus.end(),0.0);
  std::fill(NreconfMinus.begin(),NreconfMinus.end(),0.0);
  for(int iw=0; iw<nw; iw++)
  {
    wConfTilde[iw]=wConfTilde[iw]/wGlobal;
    if(wConfTilde[iw]>1.0){
       plus.push_back(iw);
       NreconfPlus[MyContext] += std::abs(wConfTilde[iw]-1.0);
    }
    else{
       minus.push_back(iw);	   
       NreconfMinus[MyContext] += std::abs(wConfTilde[iw]-1.0); 
    }
  }
  myComm->allreduce(NreconfPlus);
  myComm->allreduce(NreconfMinus);
 
  RealType CumlPlus=0.0, CumlMinus=0.0; 
  for(int i=0; i<NumContexts; i++)
  {
	  CumlPlus += NreconfPlus[i];
	  CumlMinus += NreconfMinus[i];
  }	  	 
  
  app_log() << "cumulative: " << CumlPlus << " " << CumlMinus << std::endl;

  std::fill(dN.begin(),dN.end(),0);
  dN[MyContext] = plus.size();
  dN[NumContexts+1] = static_cast<int>(std::round(CumlPlus+Random()));

  myComm->bcast(dN);
  app_log() << " # reconfigurations: " << dN[NumContexts+1] << std::endl;
  
  // app_log() << " size(plus), size(minus) " << plus.size() << " " << minus.size() << std::endl;
  // app_log() << " nreconf " << NreconfPlus << " " << NreconfMinus << std::endl;  
  
  //int Nreconf = static_cast<int>(NreconfPlus + Random());
  ////app_log() << " # of reconfigurations: " << Nreconf << std::endl;

  //std::vector<int> copyList(Nreconf,0);
  //std::vector<int> killList(Nreconf,0);
  //for(int r; r<Nreconf; r++){
  //   
  //   copyList[r] = plus[static_cast<int>(Random()*plus.size())-1];
  //   killList[r] = minus[static_cast<int>(Random()*minus.size())-1];
  //   
  //   app_log() << " copy, kill " << copyList[r] << " " << killList[r] << std::endl;

  //   W[killList[r]]->makeCopy(*(W[copyList[r]]));
  //   W[killList[r]]->ParentID=W[copyList[r]]->ID;
  //   W[killList[r]]->ID=(++NumWalkersCreated)*NumContexts+MyContext;
  //}

  return nw;
}

void MinimalReconfigurationMPI::sendWalkers(MCWalkerConfiguration& W,
    const std::vector<IndexType>& NAME)
{


}

void MinimalReconfigurationMPI::recvWalkers(MCWalkerConfiguration& W,
    const std::vector<IndexType>& NAME)
{


}

int MinimalReconfigurationMPI::branch(int iter, MCWalkerConfiguration& W, RealType trigger)
{
  //apply minimal stochastic reconfiguration to walker population
  int nwkept = reconfigureWalkers(W);

  //update EnsembleProperty
  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;

  //all weights and multiplicities are 1 in minimal stochastic reconfiguration
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end)
  {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }
  return nwkept;
}


