//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Andrew D. Baczewski, adbacze@sandia.gov, Sandia National Laboratories
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_PUREDMC_WALKER_CONTROL_H
#define QMCPLUSPLUS_PUREDMC_WALKER_CONTROL_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{

/** Class to handle walker controls with simple global sum
 *
 * Base class to handle serial mode with branching only
 */
struct PureDMC: public WalkerControlBase
{

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  PureDMC(Communicate* c):WalkerControlBase(c) {};

  /** perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, RealType trigger)
  {
    return doNotBranch(iter, W);
  };

  /** return 0.0 to disable feedback method */
  RealType getFeedBackParameter(int ngen, RealType tau)
  {
    return 0.0;
  }
};
}
#endif

