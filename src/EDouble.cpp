#include "EDouble.h"
#include <cmath>
#include <string>

EDouble::EDouble(){

}

EDouble::EDouble(double _val, double _err){
  val = _val;
  if(_err < 0.) err = _err *(-1.);
  else err = _err;
}

EDouble::EDouble(const EDouble& a){
  val = a.GetVal();
  err = a.GetError();
}

EDouble::~EDouble(){
  //dtor;
}


double EDouble::GetVal() const{
  return val;
}


double EDouble::GetError() const{
  return err;
}



void EDouble::SetVal(const double _val){
  val = _val;
}

void EDouble::SetError(const double _err){
  if(_err < 0.) err = _err *(-1.);
  else err = _err;
}


EDouble EDouble::Sqrt(){
  double relErr = err/val;
  EDouble rv( sqrt(val), 0.5*relErr*sqrt(val) );
  return rv;
}


//std::string EDouble::Latex(int relevantDigits = 0){
std::string EDouble::Latex(){
  std::string val_str = std::to_string(val);
  std::string err_str = std::to_string(err);
//  if(relevantDigits != 0){
//    std::size_t firstRelevantDigit = err_str.find_first_not_of("0");
//    if(firstRelevantDigit < val_str.size() && firstRelevantDigit < err_str.size()) firstRelevantDigit+=1; //print 2 relevant digits
//  }
  std::string rv = "${" + val_str + " \\pm " + err_str + "}$";
  return rv;
}


EDouble operator*(const EDouble& lhs, const EDouble& rhs){
  double rVal = lhs.GetVal() * rhs.GetVal();
  double rel_rErr = sqrt ( ( lhs.GetError() / lhs.GetVal() ) * ( lhs.GetError() / lhs.GetVal() )
                         + ( rhs.GetError() / rhs.GetVal() ) * ( rhs.GetError() / rhs.GetVal() ) );

  return EDouble(rVal, rel_rErr * rVal);
}


EDouble operator*(const double lhs, const EDouble& rhs){
  double rVal = lhs * rhs.GetVal();
  double rel_rErr = rhs.GetError() / rhs.GetVal();

  return EDouble(rVal, rel_rErr * rVal);
}

EDouble operator*(const EDouble& lhs, const double rhs){
  double rVal = lhs.GetVal() * rhs;
  double rel_rErr = lhs.GetError() / lhs.GetVal();

  return EDouble(rVal, rel_rErr * rVal);
}

EDouble operator/(const EDouble& lhs, const EDouble& rhs){
  double rVal = lhs.GetVal() / rhs.GetVal();
  double rel_rErr = sqrt ( ( lhs.GetError() / lhs.GetVal() ) * ( lhs.GetError() / lhs.GetVal() )
                         + ( rhs.GetError() / rhs.GetVal() ) * ( rhs.GetError() / rhs.GetVal() ) );

  if(rel_rErr < 0.) rel_rErr*=-1.;

  return EDouble(rVal, rel_rErr * rVal);
}

EDouble operator/(const double lhs, const EDouble& rhs){
  double rVal = lhs / rhs.GetVal();
  double rel_rErr = rhs.GetError() / rhs.GetVal();

  return EDouble(rVal, rel_rErr * rVal);
}

EDouble operator/(const EDouble& lhs, const double rhs){
  double rVal = lhs.GetVal() / rhs;
  double rel_rErr = lhs.GetError() / lhs.GetVal();

  return EDouble(rVal, rel_rErr * rVal);
}

EDouble operator+(const EDouble& lhs, const EDouble& rhs){
  double rVal = lhs.GetVal() + rhs.GetVal();
  double rErr = sqrt( lhs.GetError() * lhs.GetError() + rhs.GetError() * rhs.GetError() );

  return EDouble(rVal, rErr);
}

EDouble operator+(const double lhs, const EDouble& rhs){
  double rVal = lhs + rhs.GetVal();
  double rErr = rhs.GetError();

  return EDouble(rVal, rErr);
}

EDouble operator+(const EDouble& lhs, const double rhs){
  double rVal = lhs.GetVal() + rhs;
  double rErr = lhs.GetError();

  return EDouble(rVal, rErr);
}

EDouble operator-(const EDouble& lhs, const EDouble& rhs){
  double rVal = lhs.GetVal() - rhs.GetVal();
  double rErr = sqrt( lhs.GetError() * lhs.GetError() + rhs.GetError() * rhs.GetError() );

  if(rErr < 0.) rErr*=-1.;

  return EDouble(rVal, rErr);
}


EDouble operator-(const double lhs, const EDouble& rhs){
  double rVal = lhs- rhs.GetVal();
  double rErr = rhs.GetError();

  return EDouble(rVal, rErr);
}

EDouble operator-(const EDouble& lhs, const double rhs){
  double rVal = lhs.GetVal() - rhs;
  double rErr = lhs.GetError();

  return EDouble(rVal, rErr);
}



EDouble& EDouble::operator+=(const EDouble& rhs){
  EDouble tmp = *this + rhs;
  val = tmp.GetVal();
  err = tmp.GetError();
  return *this;
}

EDouble& EDouble::operator-=(const EDouble& rhs){
  EDouble tmp = *this - rhs;
  val = tmp.GetVal();
  err = tmp.GetError();
  return *this;
}

EDouble& EDouble::operator/=(const EDouble& rhs){
  EDouble tmp = *this / rhs;
  val = tmp.GetVal();
  err = tmp.GetError();
  return *this;
}

EDouble& EDouble::operator*=(const EDouble& rhs){
  EDouble tmp = *this * rhs;
  val = tmp.GetVal();
  err = tmp.GetError();
  return *this;
}


EDouble& EDouble::operator=(const EDouble& rhs){
  val = rhs.GetVal();
  err = rhs.GetError();
  return *this;
}


std::ostream& operator<<(std::ostream& lhs, const EDouble& rhs){
  lhs << rhs.GetVal() << " +/- " << rhs.GetError();
  return lhs;
}

