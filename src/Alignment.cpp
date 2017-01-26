#include "Alignment.h"

// from ROOT
#include <TLinearFitter.h>
#include <TMinuit.h>

// from here
#include "EDouble.h"

// need global variables for Minuit
std::map<std::string, std::vector<Cluster>>* FULL_DATA;
std::string* DUT;
std::vector<FibMatInfo>* MATS;
std::map<std::string, int>* NAME_AND_INDEX;

EDouble straight_line(EDouble y0, EDouble slope, double x) {
  return y0 + slope * x;
}

void fnc(int& npar, double* deriv, double& func, double param[], int flag) {
  double rv{0};
  TLinearFitter track;
  track.SetFormula("pol1");
  EDouble y0, slope;
  double z{0.};
  EDouble track_at_mat;
  for (const auto& alignMat : *(MATS)) {
    // calculate alignment for each mat
    for (unsigned int i{0}; i < FULL_DATA->at(alignMat.name).size(); ++i) {
      // collect hits that occured in each mat except for the one under test
      track.ClearPoints();
      for (const auto& mat : *(MATS)) {
        if (mat.name == *(DUT)) {
          continue;
        }
        z = mat.z_position;
        track.AddPoint(
            &z, FULL_DATA->at(mat.name).at(i).GetChargeWeightedMean() * 250. -
                    param[NAME_AND_INDEX->at(mat.name)]);
      }
      track.Eval();
      y0.SetVal(track.GetParameter(0));
      slope.SetVal(track.GetParameter(1));

      y0.SetError(track.GetParError(0));
      slope.SetError(track.GetParError(1));

      // sum up the track distances squared
      for (const auto& mat : *(MATS)) {
        track_at_mat = straight_line(y0, slope, mat.z_position);
        rv += std::pow(
            (track_at_mat.GetVal() -
             FULL_DATA->at(mat.name).at(i).GetChargeWeightedMean() * 250. -
             param[NAME_AND_INDEX->at(mat.name)]) /
                track_at_mat.GetError(),
            2);
      }
    }
  }
  std::cout << "fnc = " << rv << "\n";
  func = rv;
}

Aligner::Aligner(std::map<std::string, std::vector<Cluster>> _data,
                 std::vector<FibMatInfo> _mats) {
  data = _data;
  mats = _mats;
  int idx{0};
  // map the name of each fibre mats to an index used in the arrays TMinuit
  // takes as arguments
  for (const auto& mat : mats) {
    name_and_index[mat.name] = idx;
    ++idx;
  }
}

void Aligner::Align() {
  // make variables globally available for MINUIT
  ::FULL_DATA = &(this->data);
  ::DUT = &(this->DUT);
  ::MATS = &(this->mats);
  ::NAME_AND_INDEX = &(this->name_and_index);
  // find alignment for each mat
  for (const auto& DUT : mats) {
    this->SetDUT(DUT.name);
    std::cout << "DUT is now " << *::DUT << "\n";
    TMinuit minuit(mats.size());
    int i{0};
    // fix the position of two mats, to have a stable track,
    // dont fix the DUT though
    int n_fixed{0};
    for (const auto& mat : mats) {
      if (n_fixed < 2 && mat.name != DUT.name) {
        minuit.DefineParameter(i, ("offset_" + mat.name).c_str(), 0., 0., -0.,
                               0.);
        minuit.FixParameter(i);
        ++n_fixed;
      } else {
        minuit.DefineParameter(i, ("offset_" + mat.name).c_str(), -5000., 2000.,
                               -20000., 20000.);
      }
      ++i;
    }

    // minuit.SetFCN(ptr_fun);
    minuit.SetFCN(::fnc);
    // minuit.SetParameter(0, 0);
    minuit.Command("SIMPLEX 10000 0.1");
    minuit.Command("MIGRAD 10000 0.1");
    minuit.Command("HESSE 10000");
    double val{0.}, err{0.};
    for (const auto& mat : mats) {
      minuit.GetParameter(name_and_index[mat.name], val, err);
      alignment[DUT.name][mat.name] = val;
    }
  }
}
std::map<std::string, std::map<std::string, double>> Aligner::GetAlignment()
    const {
  return alignment;
}

std::map<std::string, int> Aligner::GetNameAndIndex() const {
  return name_and_index;
}

void Aligner::WriteToFile() const {}
void Aligner::ReadFromFile() {}
void Aligner::SetDUT(std::string _DUT) { DUT = _DUT; }

// void Aligner::fnc(int& npar, double* deriv, double& func, double param[],
//                   int flag) const {
//   double rv{0};
//   std::vector<Point> hits;
//   Line track;
//   std::cout << "starting mat loop\n";
//   for (const auto& alignMat : mats) {
//     // calculate alignment for each mat
//     std::cout << "starting evt loop\n";
//     for (unsigned int i{0}; i < data.at(alignMat.name).size(); ++i) {
//       // collect hits that occured in each mat except for the one under test
//       hits.clear();
//       for (const auto& mat : mats) {
//         if (mat.name == DUT) {
//           continue;
//         }
//         hits.push_back(
//             {mat.z_position,
//              data.at(mat.name).at(i)->GetChargeWeightedMean() * 250. -
//                  param[name_and_index.at(mat.name)]});
//       }
//       // construct the track
//       if (hits.size() != mats.size() - 1) {
//         std::cerr << "Warning: Trying to fit " << hits.size() << " hits for "
//                   << mats.size() << " mats!\n";
//       }
//       track.fitPoints(hits);
//       // sum up the track distances squared
//       for (const auto& mat : mats) {
//         rv += std::pow(
//             track.GetYforX(mat.z_position) -
//                 data.at(mat.name).at(i)->GetChargeWeightedMean() * 250. -
//                 param[name_and_index.at(mat.name)],
//             2);
//       }
//     }
//   }
//   func = rv;
// }
