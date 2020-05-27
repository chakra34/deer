#include "ComputeNEMLCPGrainGrowthOutput.h"

registerMooseObject("DeerApp", ComputeNEMLCPGrainGrowthOutput);

InputParameters ComputeNEMLCPGrainGrowthOutput::validParams() {
  InputParameters params = ComputeNEMLStressUpdate::validParams();
  params.addParam<UserObjectName>("euler_angle_provider","dummy"
                                        "Name of Euelr angle provider user object");
  params.addRequiredCoupledVar("unique_grains", "unique grains");
  params.addParam<UserObjectName>(
      "get_avg_orientation",
      "The name of the UserObject that is give the "
      "current average grain orientation");

  return params;
}

ComputeNEMLCPGrainGrowthOutput::ComputeNEMLCPGrainGrowthOutput(const InputParameters & parameters)
   : ComputeNEMLStressUpdate(parameters),
    _orientation_q(declareProperty<std::vector<Real>>("orientation_q")),
    _dislocation_density(declareProperty<Real>("dislocation_density")),
    _euler(parameters.isParamSetByUser("euler_angle_provider") ? &getUserObject<EulerAngleProvider>("euler_angle_provider") : nullptr), // May be making it a required parameter
    _grain_id(coupledValue("unique_grains")),
    _grain_id_old(coupledValueOld("unique_grains")),
    _gets_avg_ori(getUserObject<BlockAverageValue>("get_avg_orientation")),
    _changed_grain(declareProperty<Real>("changed_grain")),
    _history(declareProperty<Real>("history"))
{
  _cpmodel = static_cast<neml::SingleCrystalModel *>(_model.get());
 }

// Assigning Euler angles from file
void
ComputeNEMLCPGrainGrowthOutput::initQpStatefulProperties()
{
  ComputeNEMLStressBase::initQpStatefulProperties();
  _history[_qp] = 0.0;
  _orientation_q[_qp].resize(4);

  if (_euler != nullptr) {
      EulerAngles angles;
      // auto grains = _euler.getGrainNum(); // total grains
       _grain = std::max(0,_current_elem->subdomain_id() - 1);
       _changed_grain[_qp] = 0.0;
        angles = _euler->getEulerAngles(_grain); // current orientation
        neml:: Orientation e = neml::Orientation::createEulerAngles(angles.phi1, angles.Phi, angles.phi2,"degrees");
       _cpmodel->set_active_orientation(&_hist[_qp].front(),e);
       _orientation_q[_qp][0] = e.quat()[0];
       _orientation_q[_qp][1] = e.quat()[1];
       _orientation_q[_qp][2] = e.quat()[2];
       _orientation_q[_qp][3] = e.quat()[3];
  }
   else{
     _orientation_q[_qp][0] = 0.0;
     _orientation_q[_qp][1] = 0.0;
     _orientation_q[_qp][2] = 0.0;
     _orientation_q[_qp][3] = 1.0;
   }
  _dislocation_density[_qp] = 0.0;
}

void
ComputeNEMLCPGrainGrowthOutput::stressUpdate(
      const double * const e_np1, const double * const e_n,
      const double * const w_np1, const double * const w_n,
      double T_np1, double T_n, double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1, double * const B_np1,
      double & u_np1, double u_n, double & p_np1, double p_n)
{

  ComputeNEMLStressUpdate::stressUpdate(e_np1, e_n, w_np1, w_n, T_np1, T_n, t_np1, t_n,
               s_np1, s_n, h_np1, h_n, A_np1, B_np1, u_np1, u_n,
               p_np1, p_n);

     if (_t > 0.0){
       EulerAngles angles;
       // Real block_id;
       // block_id = std::max(0,_current_elem->subdomain_id() - 1);
       std::vector<Real> current_ori = {0.0, 0.0, 0.0, 0.0};
       if ( std::round(_grain_id[_qp]) != std::round(_grain_id_old[_qp]) ) {
           unsigned int id = (unsigned int) std::round(_grain_id[_qp]) + 1;
           current_ori = _gets_avg_ori.getBlockAvgValue(id); // change the map format
           neml:: Orientation e = neml::Orientation::Quaternion(current_ori);

           _cpmodel->set_active_orientation(&_hist[_qp].front(), e);
           _changed_grain[_qp] = 1.0;
           _hist[_qp].back()   = 0.0;
           _history[_qp]       = _hist[_qp].back();
         }
         else{
            _changed_grain[_qp] = 0.0;
            _history[_qp] = _hist[_qp].back();
          }
          getCPOutput(h_np1); // passing the history for outputs
         }
}

// Method to store CP output as material parameters
void
ComputeNEMLCPGrainGrowthOutput::getCPOutput(double * const h_np1){

  _orientation_q[_qp].resize(4);
  neml::Orientation Q = _cpmodel->get_active_orientation(h_np1);
     for (unsigned int i = 0; i < 4; i++){
       _orientation_q[_qp][i] = Q.quat()[i];   // assigning quaternion
     }
// based on Kocks-Mecking strength = b*ksi*G*sqrt(_rho); Here as an example took G of steel 80GPa, ksi = 0.9, b=2.86A
    double strength = _cpmodel->strength(h_np1,300.0);
    // Moose::out<<"strength "<<_cpmodel->strength(h_np1,c,300.0);
    // unsigned int id = _current_elem->subdomain_id();
   _dislocation_density[_qp] =  pow(strength  / (2.86*pow(10,-7) * 0.5 * 80000.0),2.0); //should be square based on formula _t*id*pow(10,6);  //
}
