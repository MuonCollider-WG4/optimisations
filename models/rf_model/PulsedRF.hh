#include <vector>
#include <memory>

#include "models/solenoid_model/DerivativesSolenoid.hh"

#ifndef PulsedRF_hh
#define PulsedRF_hh

class OnAxisFieldModel;

class PulsedRF {
public:
    PulsedRF();
    ~PulsedRF();
    PulsedRF(const PulsedRF&) = delete;
    PulsedRF& operator=(const PulsedRF&) = delete;

    PulsedRF* Clone() const;
    void GetFieldValue(const std::vector<double>& position,
                       const double& time,
                       std::vector<double>& bfield);

    bool IsReady() const;

    std::unique_ptr<OnAxisFieldModel> fieldModel;
    int maxDerivative = 1;
    double maxR = -1.0;
    double length = 0.0; // length of structure
    double zOffset = 0.0; // z position of pulse at t = 0.0, relative to centre
    double v_0 = 0.0;
    const double c_light = 299.792458; // mm/ns
private:

};

#endif // PulsedRF_hh