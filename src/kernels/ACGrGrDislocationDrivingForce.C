//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACGrGrDislocationDrivingForce.h"

#include "Material.h"

registerMooseObject("DeerApp", ACGrGrDislocationDrivingForce);

template <>
InputParameters
validParams<ACGrGrDislocationDrivingForce>()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Adds dislocation energy contribution to the Allen-Cahn equation");
  params.addRequiredParam<MaterialPropertyName>(
      "D_dislocation_density_name", "The dislocation density derivative for the specific order parameter");
  return params;
}

ACGrGrDislocationDrivingForce::ACGrGrDislocationDrivingForce(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _D_dislocation_density(getMaterialProperty<Real>("D_dislocation_density_name"))
{
}

Real
ACGrGrDislocationDrivingForce::computeDFDOP(PFFunctionType type)
{
  // shear modulus and burgers vector
  Real mu = 80000;  // 80 GPa shear modulus
  Real b = 2.86*pow(10,-10);  // Burgers vector 2.86 A

  switch (type)
  {
    case Residual:
      return 0.5 *
             _D_dislocation_density[_qp] * mu * pow(b,2); // Compute the deformation energy driving force

    case Jacobian:
      return 0.0;
  }

  mooseError("Invalid type passed in");
}
