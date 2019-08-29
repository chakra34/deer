#ifndef NEIGHBORANDVALUES_H
#define NEIGHBORANDVALUES_H

// MOOSE includes
#include "PropertyElement.h"
#include "MooseMesh.h"
#include "Material.h"
#include "neml_interface.h"

// Forward Declarations
class NeighborAndValues;

template <>
InputParameters validParams<NeighborAndValues>();

class NeighborAndValues : public PropertyElement
{
public:
  NeighborAndValues(const InputParameters & parameters);

  virtual ~NeighborAndValues();

  virtual void initialSetup() override;

  virtual void initialize() override {};

  virtual void execute() override;

  virtual void threadJoin(const UserObject & /*y*/) override {return; };

  /**
   * Called _once_ after execute has been called all all "objects".
   */
  virtual void finalize() override {return; };

  Real getMisorientationFromPair( std::tuple<unsigned int, unsigned int>) const;
  Real Misori(std::vector<Real> q1, std::vector<Real> q2) const;
  // Real Misori_neml(std::vector<Real> q1, std::vector<Real> q2) const;

protected:
  // Pointer to the mesh
//  MooseMesh *_mesh_ptr;
   unsigned int _m = 1;
   unsigned int to_do;
   std::map<std::tuple<unsigned int, unsigned int>, unsigned int> _to_do_misorientation;
// third variable to check same block id for the tuple
   std::map<std::tuple<unsigned int, unsigned int>, Real> _misorientation_map;
   std::set<std::tuple<unsigned int,unsigned int>> _unique_set_of_tuples;


};

#endif //NEIGHBORANDVALUES_H
