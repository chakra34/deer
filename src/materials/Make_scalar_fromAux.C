#include "Make_scalar_fromAux.h"

registerMooseObject("DeerApp", Make_scalar_fromAux);

template <>
InputParameters
validParams<Make_scalar_fromAux>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("val", "scalar variable");
  params.addRequiredParam<MaterialPropertyName>(
      "scalar_name", "Name of the scalar material property to be created");
  return params;
}

Make_scalar_fromAux::Make_scalar_fromAux(const InputParameters & parameters)
  : Material(parameters),
    _val(coupledValue("val")),
    _prop(declareProperty<Real>(getParam<MaterialPropertyName>("scalar_name")))
{

}
void
Make_scalar_fromAux::computeQpProperties()
{
  _prop[_qp] = _val[_qp];
}
