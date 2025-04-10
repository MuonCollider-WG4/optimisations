#include <string>
#include <cmath>
#include <iostream>
#include "PulsedRF.hh"


void PulsedRF::GetFieldValue(const std::vector<double>& position,
                             const double& time,
                             std::vector<double>& bfield) {
    if (fieldModel.get() == nullptr) {
        throw std::string("On axis field model was not set");
    }
    double r = std::sqrt(position[0]*position[0]+position[1]*position[1]);
    if (maxR > 0 && r > maxR) {
        return;
    }
    if (length > 0 && (position[2] > length/2.0 || position[2] < -length/2.0) ) {
        return;
    }

    double v_over_c2 = v_0/c_light/c_light;
    double one_minus_v2_over_c2 = 1-v_0*v_over_c2;
    std::vector<double> derivatives(maxDerivative+1, 0.0);
    for (size_t i = 0; i < derivatives.size(); i++) {
        fieldModel->GetField(position[2]-zOffset-v_0*time, i, derivatives[i]);
    }

    std::vector<double> yPower(maxDerivative+1, 1.0);
    for (size_t i = 1; i < yPower.size(); i++) {
        yPower[i] = yPower[i-1]*position[1];
    }
    bfield[5] = derivatives[0];
    double fCoeff = 1.0;
    for (size_t i = 1; i < derivatives.size(); i++) {
        if (i%2) { // odd power
            fCoeff *= 1/float(i);
            bfield[0] += +yPower[i]*v_over_c2*fCoeff*derivatives[i];
            bfield[4] += -yPower[i]*fCoeff*derivatives[i];
        } else { // even power
            fCoeff *= -one_minus_v2_over_c2/i;
            bfield[5] += yPower[i]*fCoeff*derivatives[i];
        }
    }
    //std::cerr << "PulsedRF::GetFieldValue z " << position[2] << " time " << time << " ez " << bfield[5] << std::endl;

}

PulsedRF* PulsedRF::Clone() const {
    PulsedRF* structure = new PulsedRF();
    structure->fieldModel.reset(fieldModel->Clone());
    structure->maxDerivative = maxDerivative;
    structure->maxR = maxR;
    structure->length = length; // length of structure
    structure->zOffset = zOffset; // z position of pulse at t = 0.0, relative to centre
    structure->v_0 = v_0;
    return structure;
}

PulsedRF::PulsedRF() {

}

PulsedRF::~PulsedRF() {
    
}
