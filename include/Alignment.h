#ifndef ALIGNMENT_H
#define ALIGNMENT_H

// from stl
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// from here
#include "Cluster.h"
#include "ConfigParser.h"

class Aligner {
 public:
  Aligner(std::map<std::string, std::vector<Cluster>> _data,
          std::vector<FibMatInfo> _mats);

  void Align();
  void WriteToFile() const;
  void ReadFromFile();
  void SetDUT(std::string _DUT);
  // void fnc(int& npar, double* deriv, double& func, double param[],
  //          int flag) const;
  std::map<std::string, std::map<std::string, double>> GetAlignment() const;
  std::map<std::string, int> GetNameAndIndex() const;

 private:
  std::map<std::string, std::map<std::string, double>> alignment;
  std::string DUT;
  std::map<std::string, std::vector<Cluster>> data;
  std::vector<FibMatInfo> mats;
  // map the name of each fibre mats to an index used in the arrays TMinuit
  // takes as arguments
  std::map<std::string, int> name_and_index;
};

#endif

struct Point {
  double X, Y;
};
struct Line {
  double Slope{0.}, Y0{0.};
  double GetYforX(double x) { return Slope * x + Y0; }
  // Construct line from points
  bool fitPoints(const std::vector<Point>& pts) {
    int nPoints = pts.size();
    if (nPoints < 2) {
      // Fail: infinitely many lines passing through this single point
      return false;
    }
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    for (int i = 0; i < nPoints; i++) {
      sumX += pts[i].X;
      sumY += pts[i].Y;
      sumXY += pts[i].X * pts[i].Y;
      sumX2 += pts[i].X * pts[i].X;
    }
    double xMean = sumX / nPoints;
    double yMean = sumY / nPoints;
    double denominator = sumX2 - sumX * xMean;
    // You can tune the eps (1e-7) below for your specific task
    if (std::fabs(denominator) < 1e-7) {
      // Fail: it seems a vertical line
      return false;
    }
    Slope = (sumXY - sumX * yMean) / denominator;
    Y0 = yMean - Slope * xMean;
    return true;
  }
};
