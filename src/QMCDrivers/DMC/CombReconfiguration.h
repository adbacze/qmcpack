//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Andrew D. Baczewski, adbacze@sandia.gov, Sandia National Laboratories
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_COMB_RECONFIGURATION_WALKER_CONTROL_H
#define QMCPLUSPLUS_COMB_RECONFIGURATION_WALKER_CONTROL_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{

/** Class to handle walker controls with simple global sum
 *
 * Base class to handle serial mode with branching only
 */
struct CombReconfiguration: public WalkerControlBase
{

  //random number [0,1)
  RealType UnitZeta;

  std::vector<int>      IndexCopy;
  //weight per walker
  std::vector<RealType> wConf;
  //comb
  std::vector<RealType> Zeta;
  //index 1 = total # of walkers, index 2 = size of plus, index 3 = size of minus
  std::vector<IndexType> dN;

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  CombReconfiguration(Communicate* c);

  /** perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

  /** return 0.0 to disable feedback method */
  RealType getFeedBackParameter(int ngen, RealType tau)
  {
    return 0.0;
  }

  /** return the surviving Walkers
   */
  int getIndexPermutation(MCWalkerConfiguration& W);
  int shuffleIndex(int nw);
};
}
#endif

