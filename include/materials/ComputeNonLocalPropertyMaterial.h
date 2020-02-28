#pragma once
// MOOSE includes
#include "Material.h"
#include "NeighborAndValues.h"
// Forward Declarations
class ComputeNonLocalPropertyMaterial : public Material
{
public:
  static InputParameters validParams();
  ComputeNonLocalPropertyMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

private:
   MaterialProperty<Real> & _misorientation;
   MaterialProperty<Real> & _gbenergy;

   const NeighborAndValues & _gets_misorientation;
   Real _ref_gb_energy;

};
