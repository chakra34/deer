#ifndef PROPERTYELEMENT_H
#define PROPERTYELEMENT_H

// MOOSE includes
//#include "ElementIntegralUserObject.h"
#include "ElementUserObject.h"
#include "MooseMesh.h"
#include "Material.h"

// Forward Declarations
class PropertyElement;

template <>
InputParameters validParams<PropertyElement>();

class PropertyElement : public ElementUserObject
{
public:
  PropertyElement(const InputParameters & parameters);
  virtual ~PropertyElement();
  /**
   * This is called before execute so you can reset any internal data.
   */
  virtual void initialSetup() override;
  virtual void initialize() override {};

  /**
   * Called on every "object" (like every element or node).
   * In this case, it is called at every quadrature point on every element.
   */
  virtual void execute() override;

  /**
   * Called when using threading.  You need to combine the data from "y"
   * into _this_ object.
   */
  virtual void threadJoin(const UserObject & /*y*/) override {return; };

  /**
   * Called _once_ after execute has been called all all "objects".
   */
  virtual void finalize() override {return; };
  std::vector<Real> getElementAvgValue( dof_id_type element) const;
  std::set<unsigned int> getElementNeighbor( dof_id_type element) const;
  virtual void InitializeOnce();
  virtual void MakeAverage();

protected:
  // virtual Real computeQpIntegral() override;
  // virtual std::vector<Real> computeVectorAverage(std::vector<Real> v) () override;
  // map <el id type output type> name
  unsigned int _qp ;
  unsigned int _n = 1;
  const MaterialProperty<std::vector<Real>> & _orientation;
  std::map<dof_id_type, std::vector<Real>> _element_value;
  unsigned int number_of_neighbors;
  // mapping <element which is of data type dof_id_type, neighbor set
  // Set Storing the set of neighbors
  std::set<unsigned int> _set_of_neighbors;
  std::map<dof_id_type, std::set<unsigned int>> _map_elem_with_neighbors;
};

#endif //PROPERTYELEMENT_H
