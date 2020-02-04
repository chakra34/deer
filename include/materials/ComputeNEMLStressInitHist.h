#pragma once

#include "ComputeNEMLStressInitHist.h"

class ComputeNEMLStressInitHist;

template <>
InputParameters validParams<ComputeNEMLStressInitHist>();

class ComputeNEMLStressInitHist: public ComputeNEMLStress
{
 public:
  ComputeNEMLStressInitHist(const InputParameters & parameters);

 protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;
  MaterialProperty<Real> & _history;
 private:
  Real _init_hist;

};
