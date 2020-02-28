//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputePolycrystalDislocationDensity.h"

registerMooseObject("DeerApp", ComputePolycrystalDislocationDensity);

InputParameters ComputePolycrystalDislocationDensity::validParams() {
  InputParameters params = ComputeDislocationDensityBase::validParams();
  params.addClassDescription(
      "Compute an evolving dislocation_density coupled to a grain growth phase field model.");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that keeps track of the active grains");
  params.addParam<Real>("length_scale", 1.0e-9, "Lengthscale of the problem, in meters");
  params.addParam<Real>("pressure_scale", 1.0e6, "Pressure scale of the problem, in pa");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

ComputePolycrystalDislocationDensity::ComputePolycrystalDislocationDensity(
    const InputParameters & parameters)
  : ComputeDislocationDensityBase(parameters),
    _length_scale(getParam<Real>("length_scale")),
    _pressure_scale(getParam<Real>("pressure_scale")),
    _grain_tracker(getUserObject<GrainDataTracker<Real>>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(_op_num),
    _dislocation_density(getMaterialPropertyByName<Real>("dislocation_density")),
    _D_dislocation_density(_op_num),
    _JtoeV(6.24150974e18)
{
  // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    // Initialize variables
    _vals[op_index] = &coupledValue("v", op_index);

//    declare dislocation density derivative properties
    _D_dislocation_density[op_index] = &declarePropertyDerivative<Real>(
        _dislocation_density_name, getVar("v", op_index)->name());
  }
}

void
ComputePolycrystalDislocationDensity::computeQpDislocationDensity()
{
  // Moose::out << " INSIDE ComputePolycrystalDislocationDensity first block\n";
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  // Moose::out << " after getting feature from grainTracker\n";

  // Calculate elasticity tensor
  _op_dislocation_density[_qp] = 0.0;
  Real sum_h = 0.0;
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Interpolation factor for dislocation_density
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Sum all dislocation_density
    _op_dislocation_density[_qp] += _dislocation_density[_qp] * h;
    sum_h += h;
  }

  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _op_dislocation_density[_qp] /= sum_h;

  // Calculate dislocation_density derivative: rho_deriv = dhdopi/sum_h * (rho_op - _rho)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    (*_D_dislocation_density[op_index])[_qp] = 0.0;

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;
    Real & rho_deriv = (*_D_dislocation_density[op_index])[_qp];

    rho_deriv = (_op_dislocation_density[_qp] -_dislocation_density[_qp]) * dhdopi / sum_h;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
    rho_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;
  }
}
