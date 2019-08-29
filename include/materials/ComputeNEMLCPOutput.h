#ifndef COMPUTENEMLCPOUTPUT_H
#define COMPUTENEMLCPOUTPUT_H

#include "ComputeNEMLStressUpdate.h"

class ComputeNEMLCPOutput;

template <>
InputParameters validParams<ComputeNEMLCPOutput>();

class ComputeNEMLCPOutput: public ComputeNEMLStressUpdate
{
 public:
  ComputeNEMLCPOutput(const InputParameters & parameters);
  virtual ~ComputeNEMLCPOutput() {};

 protected:
   MaterialProperty<std::vector<Real>> & _orientation_q;
   MaterialProperty<Real> & _dislocation_density;
   const bool _cpOut;

   void getCPOutput(
          const double * const e_np1, const double * const e_n,
          const double * const w_np1, const double * const w_n,
          double T_np1, double T_n, double t_np1, double t_n,
          double * const s_np1, const double * const s_n,
          double * const h_np1, const double * const h_n,
          double * const A_np1, double * const B_np1,
          double & u_np1, double u_n, double & p_np1, double p_n) override;

};


#endif // COMPUTENEMLCPOUTPUT_H
