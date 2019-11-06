//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GBEvolutionBase.h"

// Forward Declarations
class GBEvolutionField;

template <>
InputParameters validParams<GBEvolutionField>();

/**
 * Grain boundary energy parameters for isotropic uniform grain boundary energies
 */
class GBEvolutionField : public GBEvolutionBase
{
public:
  GBEvolutionField(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const MaterialProperty<Real> & _GBEnergy;
};