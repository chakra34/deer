#pragma once

#include "Material.h"

/**
 * Declares a constant material property of type RankTwoTensor.
 */
class Make_scalar_fromAux : public Material
{
public:
  static InputParameters validParams();
  Make_scalar_fromAux(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  const VariableValue & _val;
  MaterialProperty<Real> & _prop;
};
