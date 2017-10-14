//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Andrew D. Baczewski, adbacze@sandia.gov, Sandia National Laboratories
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "OhmmsData/ParameterSet.h"
#include "QMCDrivers/DMC/WalkerControlFactory.h"
#include "QMCDrivers/DMC/CombReconfiguration.h"
#include "QMCDrivers/DMC/MinimalReconfiguration.h"
#include "QMCDrivers/DMC/WalkerPureDMC.h"
#if defined(HAVE_MPI)
#include "QMCDrivers/DMC/WalkerControlMPI.h"
#include "QMCDrivers/DMC/CombReconfigurationMPI.h"
#include "QMCDrivers/DMC/MinimalReconfigurationMPI.h"
#endif

namespace qmcplusplus
{

WalkerControlBase* createWalkerController(int nwtot, Communicate* comm, xmlNodePtr cur,
    bool reconfig)
{
  app_log() << "  Creating WalkerController: target  number of walkers = " << nwtot << std::endl;
  ///set of parameters
  int nmax=0;
  std::string reconfigopt("no");
  ParameterSet m_param;
  m_param.add(nwtot,"targetWalkers","int");
  m_param.add(nwtot,"targetwalkers","int");
  m_param.add(nmax,"max_walkers","int");
  m_param.add(reconfigopt,"reconfiguration","string");
  m_param.put(cur);
  //if(nmax<0) nmax=2*nideal;
  //if(nmin<0) nmin=nideal/2;
  WalkerControlBase* wc=0;
  int ncontexts = comm->size();
  bool fixw= (reconfig || reconfigopt == "comb" || reconfigopt == "min" || reconfigopt == "pure" );
  if(fixw)
  {
    int nwloc=std::max(omp_get_max_threads(),nwtot/ncontexts);
    nwtot=nwloc*ncontexts; 
  }
#if defined(HAVE_MPI)
  if(ncontexts>1)
  {
    if(reconfigopt=="comb")
    {
      app_log() << "  Using CombReconfigurationMPI for population control." << std::endl;
      wc = new CombReconfigurationMPI(comm);
    }
    else if(reconfigopt=="min")
    {
      app_log() << "  Using MinimalReconfigurationMPI for population control." << std::endl;	    
      wc = new MinimalReconfigurationMPI(comm);
    }
    else if(reconfigopt=="pure")
    {
      app_log() << "  Using PureDMCMPI for population control." << std::endl;
    }   	   
    else
    {
      app_log() << "  Using WalkerControlMPI for dynamic population control." << std::endl;
      wc = new WalkerControlMPI(comm);
    }
  }
  else
#endif
  {
    if(reconfigopt=="comb")
    {
      app_log() << "  Using CombReconfiguration for population control." << std::endl;
      wc = new CombReconfiguration(comm);
    }
    else if(reconfigopt=="min")
    {
      app_log() << "  Using MinimalReconfiguration for population control." << std::endl;
      wc = new MinimalReconfiguration(comm);
    }
    else if(reconfigopt=="pure")
    {
      app_log() << "  Using PureDMC for population control." << std::endl;	    
    }   
    else
    {
      app_log() << "  Using WalkerControlBase for dynamic population control." << std::endl;
      wc = new WalkerControlBase(comm);
    }
  }
  wc->MyMethod=fixw;
  wc->setMinMax(nwtot,nmax);
  return wc;
}


}

