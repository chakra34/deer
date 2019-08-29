#include "ComputeNEMLCPOutput.h"

registerMooseObject("DeerApp", ComputeNEMLCPOutput);

template <>
InputParameters
validParams<ComputeNEMLCPOutput>()
{
  InputParameters params = validParams<ComputeNEMLStressUpdate>();
  params.addParam<bool>("add_cp_output", false, "Give CP output");
  return params;
}

ComputeNEMLCPOutput::ComputeNEMLCPOutput(const InputParameters & parameters)
   : ComputeNEMLStressUpdate(parameters),
   _orientation_q(declareProperty<std::vector<Real>>("orientation_q")),
   _dislocation_density(declareProperty<Real>("dislocation_density")),
    _cpOut(getParam<bool>("add_cp_output"))
{

}

void
ComputeNEMLCPOutput::getCPOutput(
       const double * const e_np1, const double * const e_n,
       const double * const w_np1, const double * const w_n,
       double T_np1, double T_n, double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1, double * const B_np1,
       double & u_np1, double u_n, double & p_np1, double p_n)
{
  // First do some declaration
  // Adding CP properties:
  // 1) orientaion in quaternion space

  // Actually call the update
  _orientation_q[_qp].resize(4);
  if (_cpOut) {
    neml::Orientation Q = _cpmodel->get_active_orientation(h_np1);
       for (unsigned int i = 0; i < 4; i++){
         _orientation_q[_qp][i] = Q.quat()[i];   // assigning quaternion
 // based on Kocks-Mecking strength = b*ksi*G*sqrt(_rho); Here as an example took G of steel 80GPa, ksi = 0.9, b=2.86A
      double strength = (double)_cpmodel->get_current_strength(h_np1);
     _dislocation_density[_qp] = _t*pow(10,6); //strength * strength / (2.86*pow(10,-10) * 0.9 * 80000.0);
       }
  }
  else {
    for (unsigned int i = 0 ; i < 4; i++){
      _orientation_q[_qp][i] = 0.0;
    }
  }

}
