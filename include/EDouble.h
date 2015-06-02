#ifndef EDOUBLE_H
#define EDOUBLE_H

#include <ostream>
#include <string>

class EDouble{
  public:
    EDouble();
    EDouble(double _val, double _err);
    EDouble(const EDouble& a);

    virtual ~EDouble();

    double GetVal() const;
    double GetError() const;

    void SetVal(const double _val);
    void SetError(const double _err);

    EDouble Sqrt();

//    std::string Latex(int relevantDigits);
    std::string Latex();

    EDouble& operator+=(const EDouble& rhs);
    EDouble& operator-=(const EDouble& rhs);
    EDouble& operator/=(const EDouble& rhs);
    EDouble& operator*=(const EDouble& rhs);

    friend EDouble operator*(const EDouble& lhs, const EDouble& rhs);
    friend EDouble operator*(const double lhs, const EDouble& rhs);
    friend EDouble operator*(const EDouble& lhs, const double rhs);

    friend EDouble operator/(const EDouble& lhs, const EDouble& rhs);
    friend EDouble operator/(const double lhs, const EDouble& rhs);
    friend EDouble operator/(const EDouble& lhs, const double rhs);

    friend EDouble operator+(const EDouble& lhs, const EDouble& rhs);
    friend EDouble operator+(const double lhs, const EDouble& rhs);
    friend EDouble operator+(const EDouble& lhs, const double rhs);

    friend EDouble operator-(const EDouble& lhs, const EDouble& rhs);
    friend EDouble operator-(const double lhs, const EDouble& rhs);
    friend EDouble operator-(const EDouble& lhs, const double rhs);


    friend std::ostream& operator<<(std::ostream& lhs, const EDouble& rhs);

    EDouble& operator=(const EDouble& rhs);



  private:
    double val;
    double err;
};

#endif // EDOUBLE_H
