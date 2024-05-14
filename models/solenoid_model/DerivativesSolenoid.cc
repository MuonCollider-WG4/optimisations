#include <math.h>
#include "DerivativesSolenoid.hh"

void FourierFieldModel::GetField(const double& z, const int& derivative, double& bzDerivative) const {
    bzDerivative = 0.0;
    for (size_t i = 0; i < harmonicContent.size(); ++i) {
        double k = 2*M_PI*(i+1)/cellLength;
        double kPow = pow(-1*k*k, derivative/2*2);
        if (derivative == 0) {
            bzDerivative += harmonicContent[i]*sin(k*z);
        } else if (derivative == 1) {
            bzDerivative += harmonicContent[i]*k*cos(k*z);
        } else if (derivative % 2 == 0) {
            bzDerivative += harmonicContent[i]*kPow*sin(k*z);
        } else {
            bzDerivative += harmonicContent[i]*kPow*k*cos(k*z);
        }
    }
}

void DerivativesSolenoid::GetField(const std::vector<double>& position,
                                   const double& /*time*/,
                                   std::vector<double>& bfield) {
    if(bfield.size() == 3) {
        throw "Error bfield was of wrong size";
    }
    double r = std::sqrt(position[0]*position[0]+position[1]*position[1]);
    double phi = atan2(position[1], position[0]);
    double z = position[2];
    double coefficient = 1;
    double br = 0.0;
    for (unsigned int i = 0; i < maxDerivative; i += 1) {
        double derivative = 0;
        if (i % 2 == 0) {
            fieldModel->GetField(z, i, derivative);
            bfield[2] += coefficient*derivative;
        } else {
            fieldModel->GetField(z, i, derivative);
            br += -coefficient/(i+1)*derivative;
        }
    }
    bfield[0] = br*cos(phi);
    bfield[1] = br*sin(phi);
}
