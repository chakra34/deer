#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "RankFourTensor.h"
#include "GuaranteeProvider.h"

class ComputeDislocationDensityBase : public DerivativeMaterialInterface<Material>,
                                    public GuaranteeProvider
{
public:
  static InputParameters validParams();
  ComputeDislocationDensityBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpDislocationDensity() = 0;

  const std::string _base_name;
  std::string _dislocation_density_name;

  MaterialProperty<Real> & _op_dislocation_density;

};
