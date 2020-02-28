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

InputParameters GBEvolutionField::validParams() {
  InputParameters params = GBEvolutionBase::validParams();
  params.addParam<Real>("ref_gb_energy",0.608," Reference GB energy in J/m^2");
  return params;
}

GBEvolutionField::GBEvolutionField(const InputParameters & parameters)
  : GBEvolutionBase(parameters),
  _GBEnergy(getMaterialPropertyByName<Real>("GB_Energy")), // in J/m^2
  _ref_gb_energy(getParam<Real>("ref_gb_energy"))
{
}

void
GBEvolutionField::computeQpProperties()
{
  // eV/nm^2
  Real val;
  val = _GBEnergy[_qp];
  if(_GBEnergy[_qp] == 0.0){
    val = _ref_gb_energy;
  }
  _sigma[_qp] = val * _JtoeV * (_length_scale * _length_scale);

  GBEvolutionBase::computeQpProperties();
}
