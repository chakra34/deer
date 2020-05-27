//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlockAverageValue.h"
#include "MooseMesh.h"
#include <cmath>
#include "libmesh/mesh_tools.h"

registerMooseObject("DeerApp", BlockAverageValue);

InputParameters BlockAverageValue::validParams() {
  InputParameters params = PropertyElement::validParams();
  params.addRequiredParam<MaterialPropertyName>("unique_grains",
          "unique grains material property");  // have to pass as a material
  return params;
}

BlockAverageValue::BlockAverageValue(const InputParameters & parameters)
  : PropertyElement(parameters),
  _grain_id(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("unique_grains")))
{
}

BlockAverageValue::~BlockAverageValue() {}

void
BlockAverageValue::initialSetup()
{
  Moose::out<<"Block Average value I am in initial set up "<<"\n";
  // Explicitly call the initialization routines for our base class
  InitializeOnce();
  initialize();
}

void
BlockAverageValue::initialize()
{

  // Moose::out<<"Block Average value I am in initialize "<<"\n";
  // Set averages to 0 for each block at each time
  const std::set<SubdomainID> & blocks = _subproblem.mesh().meshSubdomains();

  for (std::set<SubdomainID>::const_iterator it = blocks.begin(); it != blocks.end(); ++it)
  {
     // Moose::out<<" $$$$$$$ block iterator is "<<*it<<"\n";
    _integral_values[*it] = {0.0, 0.0, 0.0, 0.0};
    _volume_values[*it]   = 0;
    _average_values[*it]  = {0.0, 0.0, 0.0, 1.0};
  }
}

void
BlockAverageValue::execute()
{

  // Compute the integral on this element
  MakeAverage();

  std::vector<Real> integral_value = getElementAvgValue(_current_elem->id());
  // Moose::out<<"Integral value "<<integral_value[0]<<" "<<integral_value[1]<<" "<<integral_value[2]<<" "<<integral_value[3]<<"\n";
  unsigned int id;
  // initially grain id is 0 for all
  if (_t > 0.0){
    id = (unsigned int) _grain_id[0] + 1;  // since all qp's of the element has the same value
  }
  else{
    id = _current_elem->subdomain_id();
  }

  // Add that value to the others we've computed on this subdomain
  // grain id varies from 0 to n-1
  // Checking positive or negative dot (unity wiki)
  Real value = _integral_values[id][0]*integral_value[0] + _integral_values[id][1]*integral_value[1] +
               _integral_values[id][2]*integral_value[2] + _integral_values[id][3]*integral_value[3] ;

  if (value > 0.0){
    for (unsigned int k = 0; k < integral_value.size(); k++){
      _integral_values[id][k] += _current_elem_volume * integral_value[k];
    }
  }
  else{
    for (unsigned int k = 0; k < integral_value.size(); k++){
      _integral_values[id][k] += _current_elem_volume * -integral_value[k];
    }
  }

  // Keep track of the volume of this block
  _volume_values[id] += _current_elem_volume;
}

void
BlockAverageValue::finalize()
{
  // Loop over the integral values and sum them up over the processors
  for (auto  it = _integral_values.begin();
       it != _integral_values.end();
       ++it){
         gatherSum(it->second[0]);
         gatherSum(it->second[1]);
         gatherSum(it->second[2]);
         gatherSum(it->second[3]);
       }

  // Loop over the volumes and sum them up over the processors
  for (auto it = _volume_values.begin();
       it != _volume_values.end();
       ++it)
    gatherSum(it->second);

  // Now everyone has the correct data so everyone can compute the averages properly:
  for (auto it = _average_values.begin();
       it != _average_values.end();
       ++it)
  {
    unsigned int id = it->first;
    _average_values[id][0] = _integral_values[id][0] / _volume_values[id];
    _average_values[id][1] = _integral_values[id][1] / _volume_values[id];
    _average_values[id][2] = _integral_values[id][2] / _volume_values[id];
    _average_values[id][3] = _integral_values[id][3] / _volume_values[id];
  }
}


std::vector<Real>
BlockAverageValue::getBlockAvgValue( unsigned int grain_id) const
{

  auto obj_ptr = _average_values.find(grain_id); // pointer to the key value in the map

  std::vector<Real> ori;
  double norm_ori;

  if (obj_ptr != _average_values.end())
  {
    ori = obj_ptr->second;
    norm_ori = sqrt((ori[0]*ori[0] + ori[1]*ori[1] + ori[2]*ori[2] + ori[3]*ori[3]));
    if (norm_ori != 0.0){
      ori[0] = ori[0]/norm_ori;
      ori[1] = ori[1]/norm_ori;
      ori[2] = ori[2]/norm_ori;
      ori[3] = ori[3]/norm_ori;
      return ori;
    }
    else{
      mooseError("BlockAverageValue: norm of quaternion is 0");
      // Moose::out<<" current time "<<_t<<" grain_id "<<grain_id<<"\n";
    }
  }
  else{
    // Moose::out<<" current time "<<_t<<" grain_id "<<grain_id<<"\n";
    mooseError("BlockAverageValue: can't find the given grain");
  }
}
