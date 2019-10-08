#include "ComputeNEMLCPGrainGrowthOutput.h"

registerMooseObject("DeerApp", ComputeNEMLCPGrainGrowthOutput);

template <>
InputParameters
validParams<ComputeNEMLCPGrainGrowthOutput>()
{
  InputParameters params = validParams<ComputeNEMLStressUpdate>();
  params.addParam<UserObjectName>("euler_angle_provider","dummy"
                                        "Name of Euelr angle provider user object");
  params.addRequiredCoupledVar("unique_grains", "unique grains");
  return params;
}

ComputeNEMLCPGrainGrowthOutput::ComputeNEMLCPGrainGrowthOutput(const InputParameters & parameters)
   : ComputeNEMLStressUpdate(parameters),
    _orientation_q(declareProperty<std::vector<Real>>("orientation_q")),
    _dislocation_density(declareProperty<Real>("dislocation_density")),
    _euler(parameters.isParamSetByUser("euler_angle_provider") ? &getUserObject<EulerAngleProvider>("euler_angle_provider") : nullptr), // May be making it a required parameter
    // _grain_id(getMaterialPropertyByName<Real>("unique_grains")),
    _grain_id(coupledValue("unique_grains")),
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
  if (_euler != nullptr) {
      EulerAngles angles;
      // auto grains = _euler.getGrainNum(); // total grains
       _grain = std::max(0,_current_elem->subdomain_id() - 1);
       _changed_grain[_qp] = 0.0;
        angles = _euler->getEulerAngles(_grain); // current orientation
        neml:: Orientation e = neml::Orientation::createEulerAngles(angles.phi1, angles.Phi, angles.phi2,"degrees");
       _cpmodel->set_active_orientation(&_hist[_qp].front(),e);

       _orientation_q[_qp].resize(4);
         for (unsigned int i = 0; i < 4; i++){
           _orientation_q[_qp][i] = e.quat()[i];   // assigning quaternion
         }
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
       unsigned int block_id = std::max(0,_current_elem->subdomain_id() - 1);
       if (_grain_id[_qp] != block_id){
           angles = _euler->getEulerAngles(_grain_id[_qp]); // changed orientation
           neml:: Orientation e = neml::Orientation::createEulerAngles(angles.phi1, angles.Phi, angles.phi2,"degrees");
           _cpmodel->set_active_orientation(&_hist[_qp].front(), e);
           _changed_grain[_qp] = 1.0;
           _hist[_qp].front()  = 0.0;
           _history[_qp]       = _hist[_qp].front();
           // Moose::out<<"Current unique_grain id "<<_grain_id[_qp]<<"\n";
           // Moose::out<<"Orientation "<<angles.phi1<<" "<<angles.Phi<<" "<<angles.phi2<<"\n";
         }
         else{
            _changed_grain[_qp] = 0.0;
            _history[_qp] = _hist[_qp].front();
          }
          getCPOutput(h_np1); // passing the history for outputs
         }
   // Moose::out<<"Comparing histories "<<*h_np1<<" hist qp"<<_history[_qp]<<" \n";
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
    double strength = (double)_cpmodel->get_current_strength(h_np1);
    unsigned int id = _current_elem->subdomain_id();
   _dislocation_density[_qp] = _t*id*pow(10,6); //strength * strength / (2.86*pow(10,-10) * 0.9 * 80000.0);
}
