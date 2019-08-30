//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "RankFourTensor.h"
#include "GuaranteeProvider.h"

class ComputeDislocationDensityBase;

template <>
InputParameters validParams<ComputeDislocationDensityBase>();

/**
 * ComputeElasticityTensorBase the base class for computing elasticity tensors
 */
class ComputeDislocationDensityBase : public DerivativeMaterialInterface<Material>,
                                    public GuaranteeProvider
{
public:
  ComputeDislocationDensityBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpDislocationDensity() = 0;

  const std::string _base_name;
  std::string _dislocation_density_name;

  MaterialProperty<Real> & _op_dislocation_density;

};