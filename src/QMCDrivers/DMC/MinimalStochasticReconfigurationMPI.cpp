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
    


#include "QMCDrivers/DMC/MinimalStochasticReconfigurationMPI.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

/** default constructor
 *
 * set SwapMode
 */
MinimalStochasticReconfigurationMPI::MinimalStochasticReconfigurationMPI(Communicate* c) :WalkerControlBase(c)
{
  SwapMode=1;
}

int MinimalStochasticReconfigurationMPI::reconfigureWalkers(MCWalkerConfiguration& W)
{
  int nw(W.getActiveWalkers()); // get the number of active walkers and store it in nw
  if(wConfScaled.size()!=nw)    // check to see if instance has set size of scaled weight vector appropriately 
  {
    wConfScaled.resize(nw);     // if not (because 1st call from instance) then resize appropriately
  }

  //accumulate the energies
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
    wtot += wConfScaled[iw]=wgt; // note that we will be scaling each entry of wConfScaled later
    ++it;
  }
  curData[ENERGY_INDEX]=esum;
  curData[ENERGY_SQ_INDEX]=e2sum;
  curData[WALKERSIZE_INDEX]=nw;
  curData[WEIGHT_INDEX]=wtot;
  curData[EREF_INDEX]=ecum;
  curData[R2ACCEPTED_INDEX]=r2_accepted;
  curData[R2PROPOSED_INDEX]=r2_proposed;

  // divide the sum of weights by the number of walkers to get the global weight
  wGlobal=wtot/static_cast<RealType>(nw);
  std::vector<int> plus,minus;
  RealType NreconfPlus=0.0, NreconfMinus=0.0;
  for(int iw=0; iw<nw; iw++)
  {
    wConfScaled[iw]=wConfScaled[iw]/wGlobal;
    if(wConfScaled[iw]>1.0){
       plus.push_back(iw);
       NreconfPlus += std::abs(wConfScaled[iw]-1.0);
    }
    else{
       minus.push_back(iw);	   
       NreconfMinus += std::abs(wConfScaled[iw]-1.0); 
    }
  }
  
  // app_log() << " size(plus), size(minus) " << plus.size() << " " << minus.size() << std::endl;
  // app_log() << " nreconf " << NreconfPlus << " " << NreconfMinus << std::endl;  
  
  int Nreconf = static_cast<int>(NreconfPlus + Random());
  //app_log() << " # of reconfigurations: " << Nreconf << std::endl;

  std::vector<int> copyList(Nreconf,0);
  std::vector<int> killList(Nreconf,0);
  for(int r; r<Nreconf; r++){
     
     copyList[r] = plus[static_cast<int>(Random()*plus.size())-1];
     killList[r] = minus[static_cast<int>(Random()*minus.size())-1];
     
     app_log() << " copy, kill " << copyList[r] << " " << killList[r] << std::endl;

     W[killList[r]]->makeCopy(*(W[copyList[r]]));
     W[killList[r]]->ParentID=W[copyList[r]]->ID;
     W[killList[r]]->ID=(++NumWalkersCreated)*NumContexts+MyContext;
  }

/* 
  curData[FNSIZE_INDEX]=nw-minus.size();
  curData[RNONESIZE_INDEX]=minus.size();
  for(int i=0; i<plus.size(); i++)
  {	  
    int im=minus[i],ip=plus[i];
    W[im]->makeCopy(*(W[ip]));
    W[im]->ParentID=W[ip]->ID;
    W[im]->ID=(++NumWalkersCreated)*NumContexts+MyContext;
  }
  //int killed = shuffleIndex(nw);
  //fout << "# Total weight " << wtot << " " << killed <<  std::endl;
  //cout << "<<<< CopyIndex " << std::endl;
  //std::copy(IndexCopy.begin(), IndexCopy.end(), std::ostream_iterator<int>(std::cout, " "));
  //cout << std::endl << "<<<<<<" << std::endl;
  //for(int iw=0; iw<nw; iw++) {
  //  if(IndexCopy[iw] != iw) {
  //    W[iw]->assign(*(W[IndexCopy[iw]]));
  //  }
  //}
  return icdiff;
*/
  return nw;
}

void MinimalStochasticReconfigurationMPI::sendWalkers(MCWalkerConfiguration& W,
    const std::vector<IndexType>& NAME)
{


}

void MinimalStochasticReconfigurationMPI::recvWalkers(MCWalkerConfiguration& W,
    const std::vector<IndexType>& NAME)
{


}

int MinimalStochasticReconfigurationMPI::branch(int iter, MCWalkerConfiguration& W, RealType trigger)
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


