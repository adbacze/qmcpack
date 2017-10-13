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
    
    

#ifndef QMCPLUSPLUS_MINSTOCHRECONF_WALKER_CONTROLMPI_H
#define QMCPLUSPLUS_MINSTOCHRECONF_WALKER_CONTROLMPI_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{

/** Class to implement minimal stochastic reconfiguration as described in Assaraf, Caffarel, and Khelif PRE (2000)
 *
 * Serial implementation doesn't have to worry about swapping walkers between nodes
 */
struct MinimalStochasticReconfiguration: public WalkerControlBase
{
  ///total number of walkers
  int TotalWalkers;
  /// first index of the local walkers
  int FirstWalkers;
  /// last index of the local walkers
  int LastWalkers;
  


  // sum of individual weights divided by the total number of walkers
  RealType wGlobal;

  // weight per walker scaled by population average weight 
  std::vector<RealType> wConfScaled;

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  MinimalStochasticReconfiguration(Communicate* c);

  /** perform branch, copy positive walkers into the space occupied by negative walkers */
  int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

  /** return 0.0 to disable feedback method */
  RealType getFeedBackParameter(int ngen, RealType tau)
  {
    return 0.0;
  }
  /** return the number of surviving walkers
  */
  int reconfigureWalkers(MCWalkerConfiguration& W);

  /** send excess walkers to another node
   * @param plus local indices of the walkers to be sent
   */
  void sendWalkers(MCWalkerConfiguration& W, const std::vector<IndexType>& plus);

  /** receive excess walkers from another node, overwrite dying walkers
   * @param minus local indices of the walkers to be overwritten
   */
  void recvWalkers(MCWalkerConfiguration& W, const std::vector<IndexType>& minus);
};
}
#endif

