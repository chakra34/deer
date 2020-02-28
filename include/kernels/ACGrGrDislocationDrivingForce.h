//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ACBulk.h"

// Forward Declarations

/**
 * Calculates the porton of the Allen-Cahn equation that results from the deformation energy.
 * Must access the elastic_strain stored as a material property
 * Requires the name of the elastic tensor derivative as an input.
 */
class ACGrGrDislocationDrivingForce : public ACBulk<Real>
{
public:
  static InputParameters validParams();
  ACGrGrDislocationDrivingForce(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);

private:
  const MaterialProperty<Real> & _D_dislocation_density;
};
