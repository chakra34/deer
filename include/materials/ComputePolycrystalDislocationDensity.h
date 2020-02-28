#pragma once

#include "ComputeDislocationDensityBase.h"
#include "GrainDataTracker.h"

class ComputePolycrystalDislocationDensity : public ComputeDislocationDensityBase
{
public:
  static InputParameters validParams();
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
