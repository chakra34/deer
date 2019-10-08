#ifndef MAKE_SCALAR_FROMAUX_H
#define MAKE_SCALAR_FROMAUX_H

#include "Material.h"

class Make_scalar_fromAux;

template <>
InputParameters validParams<Make_scalar_fromAux>();

/**
 * Declares a constant material property of type RankTwoTensor.
 */
class Make_scalar_fromAux : public Material
{
public:
  Make_scalar_fromAux(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  const VariableValue & _val;
  MaterialProperty<Real> & _prop;
};

#endif // MAKE_SCALAR_FROMAUX_H
