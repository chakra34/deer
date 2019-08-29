//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GrainDataTracker.h"

class GrainTrackerDislocationDensity;
class EulerAngleProvider;

template <>
InputParameters validParams<GrainTrackerDislocationDensity>();

/**
 * Manage a list of elasticity tensors for the grains
 */
class GrainTrackerDislocationDensity : public GrainDataTracker<Real>
{
public:
  GrainTrackerDislocationDensity(const InputParameters & parameters);

protected:
  Real newGrain(unsigned int new_grain_id);

  /// generate random rotations when the Euler Angle provider runs out of data (otherwise error out)
  const bool _random_rotations;

  /// unrotated elasticity tensor
  Real _dislocation;
  /// object providing the Euler angles
  const EulerAngleProvider & _euler;

};
