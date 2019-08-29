#include "ComputeNEMLStressBase.h"


template <>
InputParameters
validParams<ComputeNEMLStressBase>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<FileName>("database", "Path to NEML XML database.");
  params.addRequiredParam<std::string>("model", "Model name in NEML database.");
  params.addCoupledVar("temperature", 0.0, "Coupled temperature");

  params.addParam<bool>("large_kinematics", false, "Use large displacement kinematics");
  params.addParam<bool>("crystal_plasticity", false, "Use crystal plasticity material models");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<UserObjectName>("euler_angle_provider","dummy"
                                        "Name of Euelr angle provider user object");
  params.addParam<unsigned int>("grain_id", 0,"ID of the grain for this material");
  return params;
}

ComputeNEMLStressBase::ComputeNEMLStressBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _fname(getParam<FileName>("database")),
    _mname(getParam<std::string>("model")),
    _mechanical_strain_inc(getMaterialPropertyByName<RankTwoTensor>("mechanical_strain_inc")),
    _vorticity_inc(getMaterialPropertyByName<RankTwoTensor>("vorticity_inc")),
    _temperature(coupledValue("temperature")),
    _temperature_old(coupledValueOld("temperature")),
    _mechanical_strain(declareProperty<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>("mechanical_strain")),
    _linear_rot(declareProperty<RankTwoTensor>("linear_rot")),
    _linear_rot_old(getMaterialPropertyOld<RankTwoTensor>("linear_rot")),
    _stress(declareProperty<RankTwoTensor>("stress")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>("stress")),
    _material_jacobian(declareProperty<RankFourTensor>("material_jacobian")),
    _hist(declareProperty<std::vector<Real>>("hist")),
    _hist_old(getMaterialPropertyOld<std::vector<Real>>("hist")),
    _energy(declareProperty<Real>("energy")),
    _energy_old(getMaterialPropertyOld<Real>("energy")),
    _dissipation(declareProperty<Real>("dissipation")),
    _dissipation_old(getMaterialPropertyOld<Real>("dissipation")),
    _elastic_strain(declareProperty<RankTwoTensor>("elastic_strain")),
    _inelastic_strain(declareProperty<RankTwoTensor>("inelastic_strain")),
    _ld(getParam<bool>("large_kinematics")),
    _cp(getParam<bool>("crystal_plasticity")),
    _euler(parameters.isParamSetByUser("euler_angle_provider") ? &getUserObject<EulerAngleProvider>("euler_angle_provider") : nullptr),
    _grain(getParam<unsigned int>("grain_id"))
{
  // I strongly hesitate to put this here, may change later
  _model = neml::parse_xml_unique(_fname, _mname);
  _cpmodel = static_cast<neml::SingleCrystalModel *>(_model.get());

  if((_cp == true) && (!parameters.isParamSetByUser("grain_id"))){
    mooseWarning("grain id's not provided, block id will be used for the cp");
    _given = 0;
  }
  if((_cp == true) && (_euler == nullptr) ) {
    mooseWarning("no euler angle file is given for a CP model default orientation will be used");
  }
}

void ComputeNEMLStressBase::computeQpProperties()
{
  // Get the strains added together
  updateStrain();

  // First do some declaration and translation
  double s_np1[6];
  double s_n[6];
  tensor_neml(_stress_old[_qp], s_n);

  double e_np1[6];
  tensor_neml(_mechanical_strain[_qp], e_np1);
  double e_n[6];
  tensor_neml(_mechanical_strain_old[_qp], e_n);

  // vorticity
  double w_np1[3];
  tensor_skew(_linear_rot[_qp], w_np1);
  double w_n[3];
  tensor_skew(_linear_rot_old[_qp], w_n);

  double t_np1 = _t;
  double t_n = _t - _dt;

  double T_np1 = _temperature[_qp];
  double T_n = _temperature_old[_qp];

  double * h_np1;
  const double * h_n;

  // MOOSE vector debug doesn't like this
  if (_model->nstore() > 0) {
    h_np1 = &(_hist[_qp][0]);
    h_n = &(_hist_old[_qp][0]);
  }
  else {
    h_np1 = nullptr;
    h_n = nullptr;
  }

  double A_np1[36];
  double B_np1[18];

  double u_np1;
  double u_n = _energy_old[_qp];

  double p_np1;
  double p_n = _dissipation_old[_qp];

  stressUpdate(e_np1, e_n, w_np1, w_n, T_np1, T_n, t_np1, t_n,
               s_np1, s_n, h_np1, h_n, A_np1, B_np1, u_np1, u_n,
               p_np1, p_n);

  // // Requesting cp output
  getCPOutput(e_np1, e_n, w_np1, w_n, T_np1, T_n, t_np1, t_n,
               s_np1, s_n, h_np1, h_n, A_np1, B_np1, u_np1, u_n,
               p_np1, p_n);

  // Do more translation, now back to tensors
  neml_tensor(s_np1, _stress[_qp]);
  recombine_tangent(A_np1, B_np1, _material_jacobian[_qp]);

  // Get the elastic strain
  double estrain[6];
  int ier = _model->elastic_strains(s_np1, T_np1, h_np1, estrain);

  if (ier != neml::SUCCESS)
    mooseError("Error in NEML call for elastic strain!");

  // Translate
  neml_tensor(estrain, _elastic_strain[_qp]);

  // For EPP purposes calculate the inelastic strain
  double pstrain[6];
  for (int i=0; i<6; i++) {
    pstrain[i] = e_np1[i] - estrain[i];
  }
  neml_tensor(pstrain, _inelastic_strain[_qp]);

  // Store dissipation
  _energy[_qp] = u_np1;
  _dissipation[_qp] = p_np1;



}

void ComputeNEMLStressBase::initQpStatefulProperties()
{
  // Basic variables maintained here
  _mechanical_strain[_qp].zero();
  _linear_rot[_qp].zero();
  _stress[_qp].zero();

  // Figure out initial history
  int ier;
  _hist[_qp].resize(_model->nstore());
  if (_model->nstore() > 0) {
    ier = _model->init_store(&_hist[_qp].front());
  }
  else {
    ier = 0;
  }

  if (ier != neml::SUCCESS) {
    mooseError("Error initializing NEML history!");
  }

  if (_euler != nullptr) {
      EulerAngles angles;
      // auto grains = _euler.getGrainNum(); // total grains
      if (_given == 0){
        _grain = std::max(0,_current_elem->subdomain_id() - 1);
      }
      angles = _euler->getEulerAngles(_grain); // current orientation
      neml:: Orientation e = neml::Orientation::createEulerAngles(angles.phi1, angles.Phi, angles.phi2,"degrees");
      // Various other junk
      _cpmodel->set_active_orientation(&_hist[_qp].front(),e);
  }
  _energy[_qp] = 0.0;
  _dissipation[_qp] = 0.0;
}

void
ComputeNEMLStressBase::updateStrain()
{
  _mechanical_strain[_qp] = _mechanical_strain_old[_qp] +
      _mechanical_strain_inc[_qp];
  _linear_rot[_qp] = _linear_rot_old[_qp] +
      _vorticity_inc[_qp];
}
