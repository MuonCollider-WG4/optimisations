#include <string>
#include <math.h>
#include <iostream>
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

void DerivativesSolenoid::GetFieldValue(const std::vector<double>& position,
                                   const double& /*time*/,
                                   std::vector<double>& bfield) {
    if (fieldModel.get() == nullptr) {
        throw std::string("On axis field model was not set");
    }
    double r = std::sqrt(position[0]*position[0]+position[1]*position[1]);
    double phi = atan2(position[1], position[0]);
    double z = position[2];
    double coefficient = 1;
    double br = 0.0;
    std::vector<double> rPow(maxDerivative+1, 1.0);
    std::vector<double> bi(maxDerivative+1, 0.0);
    for (int i = 1; i <= maxDerivative; ++i) {
        rPow[i] = rPow[i-1]*r;
    }
    for (unsigned int i = 0; i <= maxDerivative; i += 1) {
        double derivative = 0;
        if (i % 2 == 0) {
            fieldModel->GetField(z, i, derivative);
            bfield[2] += coefficient*derivative*rPow[i];
        } else {
            fieldModel->GetField(z, i, derivative);
            br += -coefficient/(i+1)*derivative*rPow[i];
        }
    }
    bfield[0] = br*cos(phi);
    bfield[1] = br*sin(phi);
}
