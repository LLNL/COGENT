#ifndef _GYROAVGOPFAC_H_
#define _GYROAVGOPFAC_H_

#include "GyroaverageOperator.H"
#include "KineticSpecies.H"

#include "NamespaceHeader.H"

namespace GyroaverageOperatorFactory {
  void createOps(std::map<std::string,GyroaverageOperator*>&, const KineticSpeciesPtrVect&);
  void deleteOps(std::map<std::string,GyroaverageOperator*>&);
}

#include "NamespaceFooter.H"

#endif
