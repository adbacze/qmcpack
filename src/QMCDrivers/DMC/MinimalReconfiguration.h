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
    
    

#ifndef QMCPLUSPLUS_MIN_RECONFIGURATION_WALKER_CONTROL_H
#define QMCPLUSPLUS_MIN_RECONFIGURATION_WALKER_CONTROL_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{

/** Class to implement minimal stochastic reconfiguration as described in Assaraf, Caffarel, and Khelif PRE (2000)
 *
 * Serial implementation doesn't have to worry about swapping walkers between nodes
 */
struct MinimalReconfiguration: public WalkerControlBase
{

  // sum of individual weights divided by the total number of walkers
  RealType wGlobal;

  // weight per walker scaled by population average weight 
  std::vector<RealType> wConfScaled;

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  MinimalReconfiguration(Communicate* c);

  /** perform branch, copy positive walkers into the space occupied by negative walkers */
  int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

  /** return 0.0 to disable feedback method */
  RealType getFeedBackParameter(int ngen, RealType tau)
  {
    return 0.0;
  }

  int reconfigureWalkers(MCWalkerConfiguration& W);
};
}
#endif

