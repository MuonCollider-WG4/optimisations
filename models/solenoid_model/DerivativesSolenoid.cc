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

void PrintVector(std::vector<double> vec) {
    for (size_t i = 0; i <= vec.size(); ++i) {
        std::cerr << vec[i] << " ";
    }
    std::cerr << std::endl;
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
    double br = 0.0;
    std::vector<double> rPow(maxDerivative+1, 1.0);
    for (int i = 1; i <= maxDerivative; ++i) {
        rPow[i] = rPow[i-1]*r;
    }

    std::vector<double> bi(maxDerivative+1, 0.0);
    for (unsigned int i = 0; i <= maxDerivative; i += 1) {
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
    for (size_t i = 1; i <= maxDerivative; i++) {
        if (i % 2 == 0) {
            bcoeff[i] = -bcoeff[i-2]/i/i;
        } else {
            acoeff[i] = -bcoeff[i-1]/(i+1);
        }
    }
}
