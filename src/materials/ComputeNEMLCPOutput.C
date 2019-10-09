#include "ComputeNEMLCPOutput.h"

registerMooseObject("DeerApp", ComputeNEMLCPOutput);

template <>
InputParameters
validParams<ComputeNEMLCPOutput>()
{
  InputParameters params = validParams<ComputeNEMLStressUpdate>();
  params.addParam<bool>("add_cp_output", false, "Give CP output");
  params.addParam<bool>("crystal_plasticity", false, "Use crystal plasticity material models");
  params.addParam<UserObjectName>("euler_angle_provider","dummy"
                                        "Name of Euelr angle provider user object");
  params.addParam<unsigned int>("grain_id", 0,"ID of the grain for this material");
  return params;
}

ComputeNEMLCPOutput::ComputeNEMLCPOutput(const InputParameters & parameters)
   : ComputeNEMLStressUpdate(parameters),
    _orientation_q(declareProperty<std::vector<Real>>("orientation_q")),
    _dislocation_density(declareProperty<Real>("dislocation_density")),
    _euler(parameters.isParamSetByUser("euler_angle_provider") ? &getUserObject<EulerAngleProvider>("euler_angle_provider") : nullptr), // May be making it a required parameter
    _grain(getParam<unsigned int>("grain_id"))
{
  _cpmodel = static_cast<neml::SingleCrystalModel *>(_model.get());
  if(!parameters.isParamSetByUser("grain_id")){
   mooseWarning("grain id's not provided, block id will be used for the cp");
   _given = 0;
  }

 if(_euler == nullptr ) {
   mooseWarning("no euler angle file is given for a single default orientation will be used !!!!");
  }

 }

// Assigning Euler angles from file
void
ComputeNEMLCPOutput::initQpStatefulProperties()
{
  ComputeNEMLStressBase::initQpStatefulProperties();
  _orientation_q[_qp].resize(4);

  if (_euler != nullptr) {
      EulerAngles angles;
      // auto grains = _euler.getGrainNum(); // total grains
      if (_given == 0){
        unsigned int grain = std::max(0,_current_elem->subdomain_id() - 1);
        angles = _euler->getEulerAngles(grain);
      }
      else{
        angles = _euler->getEulerAngles(_grain); // current orientation
      }
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
ComputeNEMLCPOutput::stressUpdate(
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

  getCPOutput(h_np1); // passing the history for outputs
}

// Method to store CP output as material parameters
void
ComputeNEMLCPOutput::getCPOutput(double * const h_np1){

  _orientation_q[_qp].resize(4);
  neml::Orientation Q = _cpmodel->get_active_orientation(h_np1);
     for (unsigned int i = 0; i < 4; i++){
       _orientation_q[_qp][i] = Q.quat()[i];   // assigning quaternion
     }
// based on Kocks-Mecking strength = b*ksi*G*sqrt(_rho); Here as an example took G of steel 80GPa, ksi = 0.9, b=2.86A
    // double strength = 0.0 ; //(double)_cpmodel->get_current_strength(h_np1);
    unsigned int id = _current_elem->subdomain_id();
   _dislocation_density[_qp] = _t*id*pow(10,6); //strength * strength / (2.86*pow(10,-10) * 0.9 * 80000.0);
}
