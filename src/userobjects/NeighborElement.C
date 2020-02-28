
#include "NeighborElement.h"

#include "MooseMesh.h"

#include "libmesh/mesh_tools.h"


registerMooseObject("DeerApp",NeighborElement);

InputParameters NeighborElement::validParams() {
  InputParameters params = ElementUserObject::validParams();
  params.addParam<unsigned int>("number_of_neighbors",1,"Number of neighbor elments");
  return params;
}

NeighborElement::NeighborElement(const InputParameters & parameters)
  : ElementUserObject(parameters)
  {
  }

NeighborElement::~NeighborElement() {}

void NeighborElement::initialize() {
  // Clear map values

  _map_elem_with_neighbors.clear();

  //Looping over all the elements
  const std::map<dof_id_type, std::vector<dof_id_type>> & node_to_elem_map = _mesh.nodeToElemMap();
  for (const auto & elem : _mesh.getMesh().element_ptr_range())
  {
    // for (unsigned int face = 0; face < elem->n_sides(); ++face)     // For 3D n_sides is n_faces
    // {
    //   // find neighbors of each face
    //   if (elem->neighbor_ptr(face))
    //     _set_of_neighbors.insert(elem->neighbor_ptr(face)->id());
    //   // later I will add 2nd or 3rd order neighbors
    //   _map_elem_with_neighbors[elem->id()] = _set_of_neighbors;  // does not include me
    // }
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


}// end initialize

std::set< unsigned int>
NeighborElement::getElementNeighbor( dof_id_type element) const
{
  auto obj_ptr = _map_elem_with_neighbors.find(element); // pointer to the key value in the map


  if (obj_ptr != _map_elem_with_neighbors.end())
  {
    return obj_ptr->second;
  }
  else
    mooseError("getElementNeighbor: can't find the given element");
}
