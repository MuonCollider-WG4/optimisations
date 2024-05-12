#include <vector>
#include <memory>

class OnAxisFieldModel {
public:
    OnAxisFieldModel() {}
    virtual ~OnAxisFieldModel() {}
    virtual void GetField(const double& z, const int& derivative, double& bzDerivative) = 0;
private:
};

class FourierFieldModel : public OnAxisFieldModel {
public:
    /** On-axis field defined by a set of fourier harmonics
     *
     *  Defines a field like Bz = Sum (b_i*sin(2*pi*b_i*z/L))
     *  2i^{th} Derivatives are (-1)^i k^2i sin(kz)
     *  (2i+1)^{th} Derivatives are (-1)^i k^2i * k cos(kz)
     */
    FourierFieldModel() : OnAxisFieldModel() {}
    virtual ~FourierFieldModel() {}
    virtual void GetField(const double& z,
                          const int& derivative,
                          double& bzDerivative) const;

    std::vector<double> harmonicContent;
    double cellLength;
private:
};

class DerivativesSolenoid {
public:
    /** Solenoid field defined by an on-axis field and its derivatives
     */
    DerivativesSolenoid() {}
    ~DerivativesSolenoid() {}

    void GetField(const std::vector<double>& position,
                  const double& time,
                  std::vector<double>& bfield);

    std::unique_ptr<OnAxisFieldModel> fieldModel;
    unsigned int maxDerivative = 1;
};