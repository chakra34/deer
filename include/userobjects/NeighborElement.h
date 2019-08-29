#ifndef NEIGHBORELEMENT_H
#define NEIGHBORELEMENT_H

// MOOSE includes
#include "ElementUserObject.h"
#include "MooseMesh.h"

// Forward Declarations
class NeighborElement;

template <>
InputParameters validParams<NeighborElement>();

class NeighborElement : public ElementUserObject
{
public:
  NeighborElement(const InputParameters & parameters);
  virtual ~NeighborElement();
  /**
   * This is called before execute so you can reset any internal data.
   */
  virtual void initialize() override;

  /**
   * Called on every "object" (like every element or node).
   * In this case, it is called at every quadrature point on every element.
   */
  virtual void execute() override{};

  /**
   * Called when using threading.  You need to combine the data from "y"
   * into _this_ object.
   */
  virtual void threadJoin(const UserObject & /*y*/) override {return; };

  /**
   * Called _once_ after execute has been called all all "objects".
   */
  virtual void finalize() override {return; };
  std::set<unsigned int> getElementNeighbor( dof_id_type element) const;
protected:
  // Pointer to the mesh
//  MooseMesh *_mesh_ptr;
  unsigned int number_of_neighbors;
  // mapping <element which is of data type dof_id_type, neighbor set

  // Set Storing the set of neighbors
  std::set<unsigned int> _set_of_neighbors;
  std::map<dof_id_type, std::set<unsigned int>> _map_elem_with_neighbors;
};

#endif //NEIGHBORELEMENT_H
