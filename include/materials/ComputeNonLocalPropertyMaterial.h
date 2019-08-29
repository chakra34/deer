#ifndef COMPUTENONLOCALPROPERTYMATERIAL_H
#define COMPUTENONLOCALPROPERTYMATERIAL_H

// MOOSE includes
#include "Material.h"
// #include "PropertyElement.h"
// #include "NeighborElement.h"
#include "NeighborAndValues.h"
// Forward Declarations
class ComputeNonLocalPropertyMaterial;

template <>
InputParameters validParams<ComputeNonLocalPropertyMaterial>();

class ComputeNonLocalPropertyMaterial : public Material
{
public:
  ComputeNonLocalPropertyMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

private:
   MaterialProperty<Real> & _misorientation;

   const NeighborAndValues & _gets_misorientation;

};

#endif //COMPUTENONLOCALPROPERTYMATERIAL_H
