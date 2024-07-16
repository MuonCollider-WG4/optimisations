#include <string>
#include <math.h>
#include <iostream>
#include <cmath>

#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_pow_int.h"

#include "DerivativesSolenoid.hh"

void FourierFieldModel::GetField(const double& z, const int& derivative, double& bzDerivative) const {
    bzDerivative = 0.0;
    for (size_t i = 0; i < harmonicList.size(); ++i) {
        double k = 2*M_PI*(i+1)/cellLength;
        double kPow = pow(-1*k*k, derivative/2);
        if (derivative == 0) {
            bzDerivative += harmonicList[i]*sin(k*z);
        } else if (derivative == 1) {
            bzDerivative += harmonicList[i]*k*cos(k*z);
        } else if (derivative % 2 == 0) {
            bzDerivative += harmonicList[i]*kPow*sin(k*z);
        } else {
            bzDerivative += harmonicList[i]*kPow*k*cos(k*z);
        }
    }
}

FourierFieldModel::~FourierFieldModel() {}

OnAxisFieldModel* FourierFieldModel::Clone() const {
    FourierFieldModel* rhs = new FourierFieldModel();
    rhs->harmonicList = harmonicList;
    rhs->cellLength = cellLength;
    return rhs;
}

void FourierFieldModel::Initialise(int /*maxDerivative*/) {
    if (cellLength == 0.0) {
        throw std::string("cellLength must be non-zero");
    }
    for (int i = harmonicList.size(); i > 0; --i) {
        if (harmonicList[i-1] != 0.0) {
            harmonicList = std::vector<double>(harmonicList.begin(), harmonicList.begin()+i);
            break;
        }
    }    
}


DerivativesSolenoid::DerivativesSolenoid(const DerivativesSolenoid& rhs) {
    if (&rhs == nullptr) {
        throw std::string("Attempt to copy a null pointer");
    }
    fieldModel.reset(rhs.fieldModel->Clone());
    maxDerivative = rhs.maxDerivative;
    maxR = rhs.maxR;
    length = rhs.length;
    SetCoeff();
}

DerivativesSolenoid* DerivativesSolenoid::Clone() const {
    return new DerivativesSolenoid(*this);
}

void DerivativesSolenoid::GetFieldValue(const std::vector<double>& position,
                                   const double& /*time*/,
                                   std::vector<double>& bfield) {
    if (fieldModel.get() == nullptr) {
        throw std::string("On axis field model was not set");
    }
    SetCoeff();
    double r = std::sqrt(position[0]*position[0]+position[1]*position[1]);
    double phi = atan2(position[1], position[0]);
    double z = position[2];
    if (maxR > 0 && r > maxR) {
        return;
    }
    if (length > 0 && (z > length/2.0 || z < -length/2.0) ) {
        return;
    }
    double br = 0.0;
    std::vector<double> rPow(maxDerivative+1, 1.0);
    for (int i = 1; i <= maxDerivative; ++i) {
        rPow[i] = rPow[i-1]*r;
    }

    std::vector<double> bi(maxDerivative+1, 0.0);
    for (int i = 0; i <= maxDerivative; i += 1) {
        double derivative = 0;
        fieldModel->GetField(z, i, derivative);
        bfield[2] += bcoeff[i]*derivative*rPow[i];
        br += acoeff[i]*derivative*rPow[i];
    }
    bfield[0] = br*cos(phi);
    bfield[1] = br*sin(phi);
}

void DerivativesSolenoid::SetCoeff() {
    if (maxDerivative >= acoeff.size()) {
        acoeff = std::vector<double>(maxDerivative+1, 0.0);
        bcoeff = std::vector<double>(maxDerivative+1, 0.0);
        bcoeff[0] = 1.0;
    } else {
        return;
    }
    for (int i = 1; i <= maxDerivative; i++) {
        if (i % 2 == 0) {
            bcoeff[i] = -bcoeff[i-2]/i/i;
        } else {
            acoeff[i] = -bcoeff[i-1]/(i+1);
        }
    }
}

std::vector< std::vector< std::vector<int> > > TanhFieldModel::_tdi =
                            std::vector< std::vector< std::vector<int> > >();

TanhFieldModel::~TanhFieldModel() {;}

void TanhFieldModel::GetField(const double& z,
                          const int& derivative,
                          double& bzDerivative) const {
    bzDerivative = _b0*(getPosTanh(z, derivative)-getNegTanh(z, derivative))/2.;
}

void TanhFieldModel::Initialise(int maxDerivative) {
    SetTanhDiffIndices(maxDerivative);
}

OnAxisFieldModel* TanhFieldModel::Clone() const {
    TanhFieldModel* tfm = new TanhFieldModel();
    tfm->_x0 = _x0;
    tfm->_lambda = _lambda;
    return tfm;
}

void TanhFieldModel::SetTanhDiffIndices(size_t n) {
    if (_tdi.size() >= n+1) {
        return;
    }
    _tdi.reserve(n+1);
    if (_tdi.size() == 0) {
        _tdi.push_back(std::vector< std::vector<int> >(1, std::vector<int>(2)));
        _tdi[0][0][0] = 1;  // 1*tanh(x) - third index is redundant
        _tdi[0][0][1] = 1;
    }
    for (size_t i = _tdi.size(); i < n+1; ++i) {
        _tdi.push_back(std::vector< std::vector<int> >());
        for (size_t j = 0; j < _tdi[i-1].size(); ++j) {
            int value = _tdi[i-1][j][1];
            if (value != 0) {
                std::vector<int> new_vec(_tdi[i-1][j]);
                new_vec[0] *= value;
                new_vec[1] -= 1;
                _tdi[i].push_back(new_vec);
                std::vector<int> new_vec2(_tdi[i-1][j]);
                new_vec2[0] *= -value;
                new_vec2[1] += 1;
                _tdi[i].push_back(new_vec2);
            }
        }
    //_tdi[i] = CompactVector(_tdi[i]);
    }
}

double TanhFieldModel::getPosTanh(double x, int n) const {
  if (n == 0) return tanh((x+_x0)/_lambda);
  double t = 0;
  double lam_n = gsl_sf_pow_int(_lambda, n);
  double tanh_x = tanh((x+_x0)/_lambda);
  for (size_t i = 0; i < _tdi[n].size(); i++)
    t += 1./lam_n*static_cast<double>(_tdi[n][i][0])
            *gsl_sf_pow_int(tanh_x, _tdi[n][i][1]);
  return t;
}

double TanhFieldModel::getNegTanh(double x, int n) const {
  if (n == 0) return tanh((x-_x0)/_lambda);
  double t = 0;
  double lam_n = gsl_sf_pow_int(_lambda, n);
  double tanh_x = tanh((x-_x0)/_lambda);
  for (size_t i = 0; i < _tdi[n].size(); i++)
    t += 1./lam_n*static_cast<double>(_tdi[n][i][0])
            *gsl_sf_pow_int(tanh_x, _tdi[n][i][1]);
  return t;
}
