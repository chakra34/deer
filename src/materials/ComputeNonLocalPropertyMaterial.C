
#include "ComputeNonLocalPropertyMaterial.h"
#include <cstdlib>


registerMooseObject("DeerApp",ComputeNonLocalPropertyMaterial);

template <>
InputParameters
validParams<ComputeNonLocalPropertyMaterial>()
{
  InputParameters params = validParams<Material>();

  params.addParam<UserObjectName>(
      "gets_misorientation",
      "The name of the UserObject that is going to find the "
      "element neighbors and misOrientation");

  return params;
}

ComputeNonLocalPropertyMaterial::ComputeNonLocalPropertyMaterial(const InputParameters & parameters)
  : Material(parameters),
  // Declare that this material is going to provide a Real
  // valued property named "_misorientation" that Kernels can use.
  _misorientation(declareProperty<Real>("misorientation")),
  // When getting a UserObject from the input file pass the name
  // of the UserObjectName _parameter_
  // Note that getUserObject returns a _const reference_ of the type in < >
  _gets_misorientation(getUserObject<NeighborAndValues>("gets_misorientation"))
  {
  }

void
ComputeNonLocalPropertyMaterial::computeQpProperties()
{
  // get neighbors of current element
  Real misori = 0.0;
  Real tempMisori = 0.0;
  unsigned int current_element = _current_elem->id();

  std::set<dof_id_type> neighbor_list = _gets_misorientation.getElementNeighbor(current_element);
  for (const auto neighbor : neighbor_list ){
    std::tuple<unsigned int, unsigned int> current_tuple;
    if ( current_element >= neighbor)
      current_tuple = std::make_tuple(current_element,neighbor);
    else
      current_tuple = std::make_tuple(neighbor,current_element);

    tempMisori = _gets_misorientation.getMisorientationFromPair(current_tuple);
    if ( tempMisori > misori)
      misori = tempMisori;
  }
  _misorientation[_qp] = misori ;
}
