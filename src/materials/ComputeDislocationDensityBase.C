//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDislocationDensityBase.h"
#include "Function.h"

template <>
InputParameters
validParams<ComputeDislocationDensityBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  return params;
}

ComputeDislocationDensityBase::ComputeDislocationDensityBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    GuaranteeProvider(this),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _dislocation_density_name(_base_name + "op_dislocation_density"),
    _op_dislocation_density(declareProperty<Real>(_dislocation_density_name))
{
}

void
ComputeDislocationDensityBase::computeQpProperties()
{
  computeQpDislocationDensity();

}
