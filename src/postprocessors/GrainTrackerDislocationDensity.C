//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrainTrackerDislocationDensity.h"
#include "EulerAngleProvider.h"

registerMooseObject("DeerApp", GrainTrackerDislocationDensity);

template <>
InputParameters
validParams<GrainTrackerDislocationDensity>()
{
  InputParameters params = validParams<GrainTracker>();
  params.addParam<bool>("random_rotations",
                        true,
                        "Generate random rotations when the Euler Angle "
                        "provider runs out of data (otherwise error "
                        "out)");
  params.addRequiredParam<Real>("dislocation", "dislocation density in grain tracker");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  return params;
}

GrainTrackerDislocationDensity::GrainTrackerDislocationDensity(const InputParameters & parameters)
  : GrainDataTracker<Real>(parameters),
    _random_rotations(getParam<bool>("random_rotations")),
    _dislocation(getParam<Real>("dislocation")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider"))
{
}

Real
GrainTrackerDislocationDensity::newGrain(unsigned int new_grain_id)
{
  EulerAngles angles;

  if (new_grain_id < _euler.getGrainNum())
    angles = _euler.getEulerAngles(new_grain_id);
  else
  {
    if (_random_rotations)
      angles.random();
    else
      mooseError("GrainTrackerDislocationDensity has run out of grain rotation data.");
  }

  Real dislocation = _dislocation * new_grain_id;

  return dislocation;
}
