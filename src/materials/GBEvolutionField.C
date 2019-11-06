//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBEvolutionField.h"

registerMooseObject("DeerApp", GBEvolutionField);

template <>
InputParameters
validParams<GBEvolutionField>()
{
  InputParameters params = validParams<GBEvolutionBase>();
  return params;
}

GBEvolutionField::GBEvolutionField(const InputParameters & parameters)
  : GBEvolutionBase(parameters),
  _GBEnergy(getMaterialPropertyByName<Real>("GB_Energy")) // in J/m^2
{
}

void
GBEvolutionField::computeQpProperties()
{
  // eV/nm^2
  Real val;
  val = _GBEnergy[_qp];
  if(_GBEnergy[_qp] == 0.0){
    val = 0.608;
  }
  _sigma[_qp] = val * _JtoeV * (_length_scale * _length_scale);

  GBEvolutionBase::computeQpProperties();
}
