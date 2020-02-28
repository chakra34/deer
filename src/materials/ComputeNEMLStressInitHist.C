#include "ComputeNEMLStressInitHist.h"

registerMooseObject("DeerApp", ComputeNEMLStressInitHist);

InputParameters ComputeNEMLStressInitHist::validParams() {
  InputParameters params = ComputeNEMLStress::validParams();
  params.addParam<Real>("init_hist",0.0," initial history value");
  return params;
}

ComputeNEMLStressInitHist::ComputeNEMLStressInitHist(const InputParameters & parameters) :
    ComputeNEMLStress(parameters),
    _history(declareProperty<Real>("history")),
    _init_hist(getParam<Real>("init_hist"))
{

}

void ComputeNEMLStressInitHist::computeQpProperties()
{
ComputeNEMLStress::computeQpProperties();
// _history[_qp] = _hist[_qp][0];

  if (_hist[_qp].size() > 0){
    _history[_qp] = _hist[_qp][0];
  }
  else{
    _history[_qp] = 0.0;
  }
}

void ComputeNEMLStressInitHist::initQpStatefulProperties()
{
  ComputeNEMLStress::initQpStatefulProperties();
  // _hist[_qp][0] = _init_hist;  // overwrite the hsitory with the set value
  if (_hist[_qp].size() > 0) {
    _hist[_qp][0] = _init_hist;  // overwrite the hsitory with the set value
  }
}
