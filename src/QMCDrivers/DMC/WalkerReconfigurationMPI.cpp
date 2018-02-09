//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Andrew D. Baczewski, adbacze@sandia.gov, Sandia National Laboratories
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/DMC/WalkerReconfigurationMPI.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

/** default constructor
 *
 * set SwapMode
 */
WalkerReconfigurationMPI::WalkerReconfigurationMPI(Communicate* c):
  WalkerControlBase(c), TotalWalkers(0)
{
  SwapMode=1;
}

int
WalkerReconfigurationMPI::branch(int iter, MCWalkerConfiguration& W, RealType trigger)
{
  int nwkept = swapWalkers(W);
  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;
  //RealType wgtInv(1.0/curData[WEIGHT_INDEX]);
  //accumData[ENERGY_INDEX]     += curData[ENERGY_INDEX]*wgtInv;
  //accumData[ENERGY_SQ_INDEX]  += curData[ENERGY_SQ_INDEX]*wgtInv;
  //accumData[WALKERSIZE_INDEX] += nwkept;
  ////accumData[WALKERSIZE_INDEX] += curData[WALKERSIZE_INDEX];
  //accumData[WEIGHT_INDEX]     += curData[WEIGHT_INDEX];
  //set Weight and Multiplicity to default values
  return nwkept;
}

int WalkerReconfigurationMPI::swapWalkers(MCWalkerConfiguration& W)
{
  //ostringstream o;
  //o << "check." << MyContext << ".dat";
  //ofstream fout(o.str().c_str(),std::ios::app);
  //
  int nw=W.getActiveWalkers(); // get the number of walkers in this context
  if(TotalWalkers !=nw*NumContexts) // each context should have the same number of walkers
  {
    FirstWalker=nw*MyContext; // get the global index of the first walker
    LastWalker=FirstWalker+nw; // get the global index of the last walker (rather that index+1, for delimiting a loop)
    TotalWalkers = nw*NumContexts;
    nwInv = 1.0/static_cast<RealType>(TotalWalkers); // this will go into the reconfigured weights (each walker will 'weight' total weight/total population)
    ncopy_w.resize(nw); // for storing the number of copies of each local walker upon reconfiguration
    wConf.resize(nw); // for storing the weights of the local walkers
    //wSum.resize(NumContexts);
    wOffset.resize(NumContexts+1); // for storing partial sums of the accumulated weight on each context
    dN.resize(NumContexts+4); // for storing the total change in population on each context (i.e., subtracting local copy/overwrite pairs)
  }

  UnitZeta=Random(); // the random offset for spacing the 'teeth' of the comb
  myComm->bcast(UnitZeta); // all contexts need the same random offset to have consistent 'teeth'
  DeltaStep=UnitZeta*nwInv;

  //std::fill(wSum.begin(),wSum.end(),0.0);
  
  // sum up the total weight of the local population of walkers
  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  int iw=0; 
  RealType esum=0.0,e2sum=0.0,wtot=0.0,ecum=0.0;
  RealType r2_accepted=0.0,r2_proposed=0.0;
  while(it != it_end)
  {
    r2_accepted+=(*it)->Properties(R2ACCEPTED);
    r2_proposed+=(*it)->Properties(R2PROPOSED);
    RealType wgt((*it)->Weight);
    wtot += wConf[iw++]=wgt;
    RealType e((*it)->Properties(LOCALENERGY));
    esum += wgt*e;
    e2sum += wgt*e*e;
    ecum += e;
    ++it;
  }

  //x //wSum[MyContext]=wtot;
  curData[ENERGY_INDEX]=esum;
  curData[ENERGY_SQ_INDEX]=e2sum;
  curData[WALKERSIZE_INDEX]=nw;
  curData[WEIGHT_INDEX]=wtot;
  curData[EREF_INDEX]=ecum;
  curData[R2ACCEPTED_INDEX]=r2_accepted;
  curData[R2PROPOSED_INDEX]=r2_proposed;
  std::fill(curData.begin()+LE_MAX,curData.end(),0.0);
  curData[LE_MAX+MyContext]=wtot;
  
  //sum up across contexts 
  myComm->allreduce(curData); // note: this implicitly reduces some quantities that we can't reduce yet (but this is a drop in the ocean)
  
  //update EnsembleProperty
  W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];
  
  //wOffset[ip] is the partial sum update to ip
  wOffset[0]=0;
  //for(int ip=0; ip<NumContexts; ip++) wOffset[ip+1]=wOffset[ip]+wSum[ip];
  
  // loop over contexts, get the accumulated weight at the boundary of each context
  for(int ip=0,jp=LE_MAX; ip<NumContexts; ip++,jp++)
    wOffset[ip+1]=wOffset[ip]+curData[jp];
  wtot=wOffset[NumContexts]; // wtot was the total weight of the local population, now it is the total weight of the entire population

  //find the lower and upper bound of index
  int minIndex=static_cast<int>((wOffset[MyContext]/wtot-DeltaStep)*static_cast<RealType>(TotalWalkers))-1;
  int maxIndex=static_cast<int>((wOffset[MyContext+1]/wtot-DeltaStep)*static_cast<RealType>(TotalWalkers))+1;
  int nb=maxIndex-minIndex+1;
  std::vector<RealType> Zeta(nb);
  for(int i=minIndex, ii=0; i<maxIndex; i++,ii++)
  {
    Zeta[ii]= wtot*(DeltaStep+static_cast<RealType>(i)*nwInv); // place the 'teeth' of the comb
  }
  RealType wCur=wOffset[MyContext];
  int ind=0;
  while(Zeta[ind]<wCur)
  {
    ind++;
  }
  //surviving walkers
  int icdiff=0;
  for(iw=0; iw<nw; iw++)
  {
    RealType tryp=wCur+std::abs(wConf[iw]);
    int ni=0;
    while(Zeta[ind]<tryp && Zeta[ind] >= wCur)
    {
      ind++;
      ni++;
    }
    wCur+=std::abs(wConf[iw]);
    if(ni)
    {
      icdiff++;
    }
    ncopy_w[iw]=ni;
  }
  //plus: a list of walkers to be duplicated
  //minus: a list of walkers to be removed
  std::vector<int> plus, minus;
  for(iw=0; iw<nw; iw++)
  {
    int m=ncopy_w[iw];
    if(m>1)
      // add the index of this walker to plus, duplicate m-1 times
    {
      plus.insert(plus.end(),m-1,iw);
    }
    else if(m==0)
      // add the walker index to be killed/overwritten
    {
      minus.push_back(iw);
    }
  }
  //save the number of local walkers to be removed. This will be collected later.
  int nw_removed=minus.size();
  //copy within the local node
  int lower=std::min(plus.size(),minus.size());
  while(lower>0)
  {
    --lower;
    int im=minus[lower]; //walker index to be replaced
    int ip=plus[lower]; //walker index to be duplicated
    W[im]->makeCopy(*(W[ip])); //copy the walker
    W[im]->ParentID=W[ip]->ID;
    W[im]->ID=(++NumWalkersCreated)*NumContexts+MyContext;
    minus.pop_back();//remove it
    plus.pop_back();//remove it
  }
  std::fill(dN.begin(),dN.end(),0);
  //dN[ip] extra/missing walkers
  if(plus.size())
  {
    dN[MyContext]=plus.size();
  }
  if(minus.size())
  {
    dN[MyContext]=-minus.size();
  }
  //dN[NumContexts] contains the number of surviving walkers
  dN[NumContexts]=icdiff;
  //other data to compute survival rate and how many walkers are swapped
  dN[NumContexts+1]=nw_removed;
  dN[NumContexts+2]=plus.size();
  dN[NumContexts+3]=minus.size();
  //collect the data
  myComm->allreduce(dN);
  //Each task will send or recv not both.
  if(plus.size())
    sendWalkers(W,plus);
  if(minus.size())
    recvWalkers(W,minus);
  //app_log() << "RECONFIGURATION  plus= " << dN[NumContexts+2] << " minus= " << dN[NumContexts+3] << std::endl;
  //app_log() << "RECONFIGURATION ";
  //for(int i=0; i<NumContexts; ++i) app_log() << dN[i] << " ";
  //app_log() << std::endl;
  //record the number of walkers created/destroyed

  // loop over all local walkers in the reconfigured population, accumulate statistics with the new weights
  it=W.begin();
  it_end=W.end();
  iw=0;
  wtot=0.0;
  esum=0.0;
  e2sum=0.0;
  ecum=0.0;
  r2_accepted=0.0;
  r2_proposed=0.0;
  while(it != it_end)
  {
    r2_accepted+=(*it)->Properties(R2ACCEPTED);
    r2_proposed+=(*it)->Properties(R2PROPOSED);
    (*it)->Weight=curData[WEIGHT_INDEX]/curData[WALKERSIZE_INDEX];
    RealType wgt((*it)->Weight);
    RealType e((*it)->Properties(LOCALENERGY));
    esum += wgt*e;
    e2sum += wgt*e*e;
    wtot += wgt;
    ecum += e;
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
  myComm->allreduce(curData);
  W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];

  curData[RNONESIZE_INDEX]=dN[NumContexts+1];
  curData[FNSIZE_INDEX]=curData[WEIGHT_INDEX]-curData[RNONESIZE_INDEX];
  curData[SENTWALKERS_INDEX]=dN[NumContexts+2];
  //collect surviving walkers
  return dN[NumContexts];
}

void WalkerReconfigurationMPI::sendWalkers(MCWalkerConfiguration& W,
    const std::vector<IndexType>& plus)
{
  std::vector<int> minusN, plusN;
  for(int ip=0; ip<NumContexts; ip++)
  {
    if(dN[ip]>0)
    {
      plusN.insert(plusN.end(),dN[ip],ip);
    }
    else
      if(dN[ip]<0)
      {
        minusN.insert(minusN.end(),-dN[ip],ip);
      }
  }
  int wbuffer_size=W[0]->byteSize();
  int nswap=plusN.size();
  int last = std::abs(dN[MyContext])-1;
  int ic=0;
  while(ic<nswap && last>=0)
  {
    if(plusN[ic]==MyContext)
    {
      //OOMPI_Packed sendBuffer(wbuffer_size,OOMPI_COMM_WORLD);
      OOMPI_Packed sendBuffer(wbuffer_size,myComm->getComm());
      W[plus[last]]->putMessage(sendBuffer);
      //OOMPI_COMM_WORLD[minusN[ic]].Send(sendBuffer);
      myComm->getComm()[minusN[ic]].Send(sendBuffer);
      --last;
    }
    ++ic;
  }
}

void WalkerReconfigurationMPI::recvWalkers(MCWalkerConfiguration& W,
    const std::vector<IndexType>& minus)
{
  std::vector<IndexType> minusN, plusN;
  for(int ip=0; ip<NumContexts; ip++)
  {
    if(dN[ip]>0)
    {
      plusN.insert(plusN.end(),dN[ip],ip);
    }
    else
      if(dN[ip]<0)
      {
        minusN.insert(minusN.end(),-dN[ip],ip);
      }
  }
  int wbuffer_size=W[0]->byteSize();
  int nswap=plusN.size();
  int last = std::abs(dN[MyContext])-1;
  int ic=0;
  while(ic<nswap && last>=0)
  {
    if(minusN[ic]==MyContext)
    {
      //OOMPI_Packed recvBuffer(wbuffer_size,OOMPI_COMM_WORLD);
      //OOMPI_COMM_WORLD[plusN[ic]].Recv(recvBuffer);
      OOMPI_Packed recvBuffer(wbuffer_size,myComm->getComm());
      myComm->getComm()[plusN[ic]].Recv(recvBuffer);
      int im=minus[last];
      W[im]->getMessage(recvBuffer);
      W[im]->ParentID=W[im]->ID;
      W[im]->ID=(++NumWalkersCreated)*NumContexts+MyContext;
      --last;
    }
    ++ic;
  }
}


