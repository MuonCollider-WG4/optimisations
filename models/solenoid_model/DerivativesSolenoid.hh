#ifndef DerivativesSolenoid_hh
#define DerivativesSolenoid_hh

#include <vector>
#include <memory>

class OnAxisFieldModel {
public:
    OnAxisFieldModel() {}
    virtual ~OnAxisFieldModel() {}
    virtual OnAxisFieldModel* Clone() const = 0;
    virtual void GetField(const double& z, const int& derivative, double& bzDerivative) const = 0;
    virtual void Initialise() = 0;
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
    virtual ~FourierFieldModel();
    virtual void GetField(const double& z,
                          const int& derivative,
                          double& bzDerivative) const;
    virtual void Initialise();
    OnAxisFieldModel* Clone() const;
    std::vector<double> harmonicList;
    double cellLength;
private:
    FourierFieldModel(const FourierFieldModel& ffm);
};

class DerivativesSolenoid {
public:
    /** Solenoid field defined by an on-axis field and its derivatives
     */
    DerivativesSolenoid() {}
    ~DerivativesSolenoid() {}
    DerivativesSolenoid(const DerivativesSolenoid& rhs);

    DerivativesSolenoid* Clone() const;
    void GetFieldValue(const std::vector<double>& position,
                  const double& time,
                  std::vector<double>& bfield);

    bool IsReady() const;

    void SetCoeff();
    std::unique_ptr<OnAxisFieldModel> fieldModel;
    int maxDerivative = 1;
    double maxR = -1.0;
    double length = 0.0;
private:
    std::vector<double> acoeff = std::vector<double>(1, 0.0);
    std::vector<double> bcoeff = std::vector<double>(1, 1.0);
};

#endif //DerivativesSolenoid_hh
