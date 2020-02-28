//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PropertyElement.h"

#include "MooseMesh.h"

#include "libmesh/mesh_tools.h"


registerMooseObject("DeerApp",PropertyElement);

InputParameters PropertyElement::validParams() {
  InputParameters params = ElementUserObject::validParams();
  params.addRequiredParam<MaterialPropertyName>("orientation_q",
          "The name of the material property (orientation)retrieve per element");
  return params;
}

PropertyElement::PropertyElement(const InputParameters & parameters)
  : ElementUserObject(parameters),_qp(0),
  _orientation(getMaterialPropertyByName<std::vector<Real>>(getParam<MaterialPropertyName>("orientation_q")))
  {
  }

PropertyElement::~PropertyElement() {}

void PropertyElement::initialSetup() {
  // Clear map values
  InitializeOnce();
}

void PropertyElement::execute()
  {
    // Compute the integral on this element
   MakeAverage();
  }


void PropertyElement::InitializeOnce() {
  // Clear map values
  if (_n == 1){  // Doing this thing only once
  _map_elem_with_neighbors.clear();
    _element_value.clear();

    const std::map<dof_id_type, std::vector<dof_id_type>> & node_to_elem_map = _mesh.nodeToElemMap();
    for (const auto & elem : _mesh.getMesh().element_ptr_range())
    {
      _element_value[elem->id()] = {0.0,0.0,0.0,1.0};

      /// based on number of nodes ///
      _set_of_neighbors.clear();
      for (unsigned int node = 0; node < elem->n_nodes(); node++){
        unsigned int global_node_id = elem->node_id(node);
        auto list_connected_elem = node_to_elem_map.find(global_node_id);
        std::vector<dof_id_type> connected_elements = list_connected_elem->second;
        for(unsigned int  v = 0; v < connected_elements.size(); v++){
          _set_of_neighbors.insert(connected_elements[v]);
        }
      }
      _map_elem_with_neighbors[elem->id()] = _set_of_neighbors;  // includes me
    }
       _n += 1;
  }
}

void PropertyElement::MakeAverage()
{
  _vector_integral_value.resize(_orientation[0].size());

  for (unsigned int i = 0; i < _orientation[0].size(); i++){
    _vector_integral_value[i] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); _qp++){
      _vector_integral_value[i] += _JxW[_qp] * _coord[_qp] * _orientation[_qp][i];
    }
    _vector_integral_value[i] /= _current_elem_volume;
  }
  auto obj_ptr = _element_value.find(_current_elem->id());
  if (obj_ptr != _element_value.end()){
    obj_ptr->second = _vector_integral_value;
  }
  else{
    mooseError("PropertyElement: can't find the given element");
  }

}

std::vector<Real>
PropertyElement::getElementAvgValue( dof_id_type element) const
{

  auto obj_ptr = _element_value.find(element); // pointer to the key value in the map


  if (obj_ptr != _element_value.end())
  {
    return obj_ptr->second;
  }
  else{
    mooseError("getElementAvgValue: can't find the given element");
    // Moose::out<<"getElementAvgValue: can't find the given element\n";
    // return {0.0,0.0,0.0,0.0};
  }
}

std::set< unsigned int>
PropertyElement::getElementNeighbor( dof_id_type element) const
{
  auto obj_ptr = _map_elem_with_neighbors.find(element); // pointer to the key value in the map

  if (obj_ptr != _map_elem_with_neighbors.end())
  {
    return obj_ptr->second;
  }
  else{
    mooseError("getElementNeighbor: can't find the given element");

  }
}
