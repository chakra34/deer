#ifndef BLOCKAVERAGEVALUE_H
#define BLOCKAVERAGEVALUE_H

// #include "ElementUserObject.h"
#include "MooseMesh.h"
#include "Material.h"
#include "PropertyElement.h"
#include "libmesh/mesh_tools.h"

// Forward Declarations
class BlockAverageValue;

template <>
InputParameters validParams<BlockAverageValue>();

/**
 * Computes the average value of a variable on each block
 */
class BlockAverageValue : public PropertyElement
{
public:
  BlockAverageValue(const InputParameters & parameters);
  virtual ~BlockAverageValue();

  /**
   * Given a block ID return the average value for a variable on that block
   *
   * Note that accessor functions on UserObjects like this _must_ be const.
   * That is because the UserObject system returns const references to objects
   * trying to use UserObjects.  This is done for parallel correctness.
   *

  /**
   * This is called before execute so you can reset any internal data.
   */
  virtual void initialSetup() override ;
  virtual void initialize() override;

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
  virtual void finalize() override;

  std::vector<Real> getBlockAvgValue( unsigned int grain_id) const;

protected:
  const MaterialProperty<Real> & _grain_id;
  // const VariableValue & _grain_id;
  // This map will hold the partial sums for each block
  std::map<unsigned int, std::vector<Real>> _integral_values;

  // This map will hold the partial volume sums for each block
  std::map<unsigned int, Real> _volume_values;

  // This map will hold our averages for each block
  std::map<unsigned int, std::vector<Real>> _average_values;
};

#endif //BLOCKAVERAGEVALUE_H
