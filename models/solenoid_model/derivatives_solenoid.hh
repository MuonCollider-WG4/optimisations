#include <vector>

class OnAxisFieldModel {
public:
    OnAxisFieldModel() {}
    virtual ~OnAxisFieldModel() {}
    virtual void GetField(const double& z, const int& derivative, double& bzDerivative) = 0;
private:

};

class FourierFieldModel {
public:
    /** On-axis field defined by a set of fourier harmonics
     *
     *  Defines a field like Bz = Sum (b_i*sin(2*pi*b_i*z/L))
     *  Derivative
     */
    FourierFieldModel() {}
    virtual ~FourierFieldModel() {}
    virtual void GetField(const double& z,
                          const int& derivative,
                          double& bzDerivative) const;



    std::vector<double> harmonicContent;
    double cellLength;
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