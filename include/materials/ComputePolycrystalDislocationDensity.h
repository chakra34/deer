//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeDislocationDensityBase.h"
#include "GrainDataTracker.h"

// Forward Declarations
class ComputePolycrystalDislocationDensity;

template <>
InputParameters validParams<ComputePolycrystalDislocationDensity>();

/**
 * Compute an evolving elasticity tensor coupled to a grain growth phase field model.
 */
class ComputePolycrystalDislocationDensity : public ComputeDislocationDensityBase
{
public:
  ComputePolycrystalDislocationDensity(const InputParameters & parameters);

protected:
  virtual void computeQpDislocationDensity();

  Real _length_scale;
  Real _pressure_scale;

  /// Grain tracker object
  const GrainDataTracker<Real> & _grain_tracker;

  /// Number of order parameters
  const unsigned int _op_num;

  /// Order parameters
  std::vector<const VariableValue *> _vals;

  /// actual Material property dislocation density
  const MaterialProperty<Real> & _dislocation_density;

  /// vector of dislocation density material properties
  std::vector<MaterialProperty<Real>* > _D_dislocation_density;

  /// Conversion factor from J to eV
  const Real _JtoeV;

};
