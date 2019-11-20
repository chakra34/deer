
#include "NeighborAndValues.h"

#include "MooseMesh.h"

#include "libmesh/mesh_tools.h"

#include<cstdlib>

registerMooseObject("DeerApp",NeighborAndValues);

template <>
InputParameters
validParams<NeighborAndValues>()
{
  InputParameters params = validParams<PropertyElement>();
  return params;
}

NeighborAndValues::NeighborAndValues(const InputParameters & parameters)
  : PropertyElement(parameters)
  {
  }

NeighborAndValues::~NeighborAndValues() {}

void NeighborAndValues::initialSetup()
{
  // Clear map values
  InitializeOnce();
  if (_m == 1)
  {  // Doing this thing only once
    _misorientation_map.clear();
    _unique_set_of_tuples.clear();
    for (const auto & elem : _mesh.getMesh().element_ptr_range())
    {
      unsigned int current_element = elem->id();
      std::set<dof_id_type> myNeighbors = getElementNeighbor(current_element);
      for (auto neighbor : myNeighbors)
      {
        std::tuple<unsigned int, unsigned int> current_tuple;
        if (current_element >= neighbor)
          current_tuple = std::make_tuple(current_element,neighbor);
        else
          current_tuple = std::make_tuple(neighbor,current_element);
        _unique_set_of_tuples.insert(current_tuple);
      }
    }
    for (auto member : _unique_set_of_tuples)
    {
      _misorientation_map[member] = 0.0;
    }
  _m += 1;
  }
}

void NeighborAndValues::execute()
  {
    MakeAverage();
  }

void NeighborAndValues::finalize()
{
  Moose::out<<"In the finalize block \n";
  for (auto it = _misorientation_map.begin(); it != _misorientation_map.end(); ++it)
  {
    auto current_tuple = it->first;
    unsigned int element1 = std::get<0>(current_tuple);
    unsigned int element2 = std::get<1>(current_tuple);
    auto q1 = getElementAvgValue(element1);
    auto q2 = getElementAvgValue(element2);
    Real misori;
    if (element1 == element2)  {
      misori = 0.0;
      it->second = misori;
    }
    else{
      misori = Misori_neml(q1,q2);
      // Real misori = Misori(q1,q2);
      it->second = misori;
    }
  }

  // std::vector<neml::Orientation> A;
  // std::vector<neml::Orientation> B;
  // std::vector<Real> q1;
  // std::vector<Real> q2;
  // unsigned int element1;
  // unsigned int element2;
  // for (auto it = _misorientation_map.begin(); it != _misorientation_map.end(); ++it)
  // {
  //   q1.resize(4);
  //   q2.resize(4);
  //   q1 = {0.0,0.0,0.0,0.0};
  //   q2 = {0.0,0.0,0.0,0.0};
  //   element1 = 0;
  //   element2 = 0;
  //   auto current_tuple = it->first;
  //   element1 = std::get<0>(current_tuple);
  //   element2 = std::get<1>(current_tuple);
  //   q1 = getElementAvgValue(element1);
  //   q2 = getElementAvgValue(element2);
  //   A.push_back(q1);
  //   B.push_back(q2);
  // }
  // auto misori_angles = Misori_block_neml(A,B);
  // unsigned int index = 0;
  //
  // for (auto it = _misorientation_map.begin(); it != _misorientation_map.end(); ++it)
  // {
  //     it->second = misori_angles[index];
  //     index += 1;
  //   }

}

Real
NeighborAndValues::getMisorientationFromPair( std::tuple<unsigned int, unsigned int> elem_tuples) const
{

  auto obj_ptr = _misorientation_map.find(elem_tuples); // pointer to the key value in the map


  if (obj_ptr != _misorientation_map.end())
  {
    return obj_ptr->second;
  }
  else{
    // Moose::out<< "getMisorientationFromPair: can't find the given element\n";
    // return 0.0;
    mooseError("getMisorientationFromPair: can't find the given element");
  }
}
Real
NeighborAndValues::Misori(const std::vector<Real> q1, const std::vector<Real> q2) const
{
  std::vector<Real> result(4,0.0);
  std::vector<Real> normresult(4,0.0);
  std::vector<Real> conj_q2(4,0.0);
  conj_q2[0] = q2[0];
  conj_q2[1] = -q2[1];
  conj_q2[2] = -q2[2];
  conj_q2[3] = -q2[3];

  result[1] =  q1[1] * conj_q2[0] + q1[2] * conj_q2[3] - q1[3] * conj_q2[2] + q1[0] * conj_q2[1];
  result[2] = -q1[1] * conj_q2[3] + q1[2] * conj_q2[0] + q1[3] * conj_q2[1] + q1[0] * conj_q2[2];
  result[3] =  q1[1] * conj_q2[2] - q1[2] * conj_q2[1] + q1[3] * conj_q2[0] + q1[0] * conj_q2[3];
  result[0] = -q1[1] * conj_q2[1] - q1[2] * conj_q2[2] - q1[3] * conj_q2[3] + q1[0] * conj_q2[0];
  Real r = std::sqrt(result[0]*result[0] +
                     result[1]*result[1] +
                     result[2]*result[2] +
                     result[3]*result[3] );
  // Convert into unit vector
  if (r > 0){
    normresult[0] = result[0] / r;
    normresult[1] = result[1] / r;
    normresult[2] = result[2] / r;
    normresult[3] = result[3] / r;
  }
  else
    normresult = result;

  Real misori = abs(2.0 * acos (result[0])) * 180.0/3.14;   // misOrientation = 2*acos(qw)
  return misori;
}

Real
NeighborAndValues::Misori_neml(const std::vector<Real> q1, const std::vector<Real> q2) const
{
  std::vector<Real> r1(4,0.0);
  std::vector<Real> r2(4,0.0);
  r1.assign(q1.begin(), q1.end());
  r2.assign(q2.begin(), q2.end());
  Real axis[3];
  Real angle;
  neml::Orientation Q1(r1);
  neml::Orientation Q2(r2);
  neml::SymmetryGroup cubic("432");  // cubic: 432
  neml::Orientation misori_n = cubic.misorientation(Q1,Q2);
  misori_n.to_axis_angle(axis,angle);
  return angle * 180.0/3.14;
}

std::vector<Real>
NeighborAndValues::Misori_block_neml(const std::vector<neml::Orientation> Q1,
                                     const std::vector<neml::Orientation> Q2) const
{

  neml::SymmetryGroup cubic("432");  // cubic: 432
  auto misori_n = cubic.misorientation_block(Q1,Q2);
  std::vector<Real> angles;
  unsigned int index = 0;
  for (auto it : misori_n){
    Real axis[3];
    Real angle;
    it.to_axis_angle(axis,angle);
    angles.push_back(angle*180.0/3.14);
    index += 1;
  }
  return angles;
}
