#pragma once

#include "ComputeNEMLStress.h"

class ComputeNEMLStressInitHist: public ComputeNEMLStress
{
 public:
   static InputParameters validParams();
  ComputeNEMLStressInitHist(const InputParameters & parameters);

 protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;
  MaterialProperty<Real> & _history;
 private:
  Real _init_hist;

};
