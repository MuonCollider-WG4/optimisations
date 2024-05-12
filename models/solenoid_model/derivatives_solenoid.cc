#include "DerivativesSolenoid.hh"

OnAxisFieldModel();
~OnAxisFieldModel() {}

FourierFieldModel::~FourierFieldModel() {

}

void FourierFieldModel::GetField(const double& z, const int& derivative, double& bzDerivative) {
    double b0 = 0.0
    for (size_t i = 0; i < harmonicContent.size(); ++i) {
        double k = 2*PI*(i+1)/cellLength;
        double kPow = pow(-1*k*k, derivative);
        if (derivative % 2 == 0) {
            bzDerivative += harmonicContent[i]*kPow*sin(k*z);
        }
    }
}

class FourierFieldModel {
public:
    FourierFieldModel();
    virtual ~FourierFieldModel() {}
    virtual void GetField(const double& z,
                          const int& derivative,
                          double& bzDerivative)

    std::vector<double> harmonicContent;
private:
};

class DerivativesSolenoid {
    DerivativesSolenoid();
    ~DerivativesSolenoid();

    void GetField(const std::vector<double>& position,
                  const double& time,
                  std::vector<double>& bfield) {

    }

};