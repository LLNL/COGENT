#include "GyroaverageOperatorFactory.H"

#include "NamespaceHeader.H"

void GyroaverageOperatorFactory::createOps(std::map<std::string, GyroaverageOperator*>&  a_ops,
                                           const KineticSpeciesPtrVect&                  a_species_vec )
{
  /* clear any previous operators */
  deleteOps(a_ops);

  /* number of species */
  int nspecies = a_species_vec.size();

  for (int i=0; i<nspecies; i++) {
    const KineticSpecies& species = *(a_species_vec[i]);

    if (species.isGyrokinetic()) {

      /* get species characteristics */
      const std::string&  name = species.name();
      const Real          mass = species.mass(),
                          charge = species.charge();
  
      /* create the gyroaveraging operator for this species
       * if it doesn't already exist */
      std::map<std::string, GyroaverageOperator*>::iterator it;
      it = a_ops.find(name);
      if (it == a_ops.end()) {
        GyroaverageOperator* op = new GyroaverageOperator();
        op->define(species.phaseSpaceGeometry(), name, mass, charge);
        a_ops[name] = op;
      }
    }
  }

  return;
}

void GyroaverageOperatorFactory::deleteOps(std::map<std::string, GyroaverageOperator*>&  a_ops )
{
  for ( std::map<std::string, GyroaverageOperator*>::iterator it = a_ops.begin();
        it != a_ops.end(); ++it) {
    delete it->second;
  }
  a_ops.clear();
}

#include "NamespaceFooter.H"
