#ifndef ComputeNEMLCPGrainGrowthOutput_H
#define ComputeNEMLCPGrainGrowthOutput_H

#include "ComputeNEMLStressUpdate.h"
#include "Material.h"

class ComputeNEMLCPGrainGrowthOutput;

template <>
InputParameters validParams<ComputeNEMLCPGrainGrowthOutput>();

class ComputeNEMLCPGrainGrowthOutput: public ComputeNEMLStressUpdate
{
 public:
  ComputeNEMLCPGrainGrowthOutput(const InputParameters & parameters);
  virtual ~ComputeNEMLCPGrainGrowthOutput() {};

 protected:
   neml::SingleCrystalModel *_cpmodel = nullptr;

   MaterialProperty<std::vector<Real>> & _orientation_q;
   MaterialProperty<Real> & _dislocation_density;
   /// object providing the Euler angles
   const EulerAngleProvider * _euler;
   /// grain id
  const VariableValue & _grain_id;
   // const MaterialProperty<Real> & _grain_id;

   MaterialProperty<Real> & _changed_grain;
   MaterialProperty<Real> & _history;
   unsigned int _grain;

    virtual void initQpStatefulProperties();

   virtual void stressUpdate(
       const double * const e_np1, const double * const e_n,
       const double * const w_np1, const double * const w_n,
       double T_np1, double T_n, double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1, double * const B_np1,
       double & u_np1, double u_n, double & p_np1, double p_n);

   void getCPOutput(double *const h_np1);
};


#endif // ComputeNEMLCPGrainGrowthOutput_H
