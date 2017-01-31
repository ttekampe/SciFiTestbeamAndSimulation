#ifndef CALIBRATION_H
#define CALIBRATION_H

// from std
#include <map>
#include <string>
#include <vector>

// from ROOT
#include "TTree.h"

// from here
#include "Cluster.h"

typedef std::vector<Channel> Event;

struct CalibrationFiles {
  std::string pedestal;
  std::string gain;
};

struct calibrationRunNumbers {
  unsigned int dark;
  unsigned int led;
};

typedef std::map<unsigned int, double> map_uint_double;
typedef std::map<unsigned int, map_uint_double> map_uint_map_uint_double;

/** \brief stores the branches of a TTree in stl containers.
 *
 * The parameters
 * factor, gains and pedestals are optional. If they are passed the adc value is
 * corrected accordingly.
 *
 * \param dataTree pointer to the TTree that is read
 * \param uplinkMin the lowest uplink number that is read
 * \param uplinkMax the highest uplink number that is read
 * \param nAdcs the number of adcs per uplink
 * \param factor a number that is multiplied to the (corrected) adc value
 * \param pedestals the pedestal for each channel of each adc
 * \param gains the gain for each channel of each adc
 * \return a vector of events. Each event consists of the adc values read from
 * the TTree
 */
std::vector<Event> parseRootTree(
    TTree* dataTree, unsigned int uplinkMin, unsigned int uplinkMax,
    unsigned int nAdcs, double factor = 1.,
    const map_uint_map_uint_double& pedestals = map_uint_map_uint_double(),
    const map_uint_map_uint_double& gains = map_uint_map_uint_double());

/** \brief Fits the adc distribution for all channesls and calculates the gains
 *
 * \param t pointer to the TTree that contains testbeam data
 * \param uplinkMin the lowest uplink number that is read
 * \param uplinkMax the highest uplink number that is read
 * \param adcIDmin the smallest number of the adcs
 * \param adcIDmax the largest number of the adcs
 * \param maxGaussians The maximum number of peaks that are fit per adc
 * \param fileName The name of the file the plots and gains are saved in
 * \param savePath Path under which the results are stored
 */
void produceGains(TTree* t, const unsigned int uplinkMin,
                  const unsigned int uplinkMax, const unsigned int adcIDmin,
                  const unsigned int adcIDmax, const unsigned int maxGaussians,
                  TString fileName, std::string savePath);

/** \brief Reads the gains from a file
*
* \param fileName the file containing the gains
* \return A map of maps with the structure <uplink, <Adc, gain>>
*/
map_uint_map_uint_double readGains(std::string fileName);

/** \brief Looks up matching dark and LED runs for the given run number
*
* \param runNumber the run number of the data run
* \param catalogueFileName the file containing three columns (dark, LED, data)
* \return The gains of all adc channels of all uplinks
*/
calibrationRunNumbers lookUpCalibrationFiles(
    const unsigned int runNumber, const std::string catalogueFileName);

/** \brief Extracts the run number from a filename
*
* \param filename the name of a testbeam data file
* \return The run number extracted from the file name
*/
unsigned int runNumberFromFilename(std::string filename);

/** \brief calculates the pedestals from a dark run file
*
* \param fileName the name of the dark run file
* \param uplinkMin the lowest uplink number
* \param uplinkMax the highest uplink number
* \param nAdcs the number of adcs per uplink
* \return A map of maps with the structure <uplink, <Adc, pedestal>>
*/
map_uint_map_uint_double getPedestals(const std::string fileName,
                                      const unsigned int uplinkMin,
                                      const unsigned int uplinkMax,
                                      const unsigned int nAdcs);

/** \brief produces a root file that is pedestal and gain corrected
*
* \param tree2correct data run TTree shall be corrected
* \param gains the gais of all adc channels
* \param pedestals the pedestals of all adc channels
* \param uplinkMin the lowest uplink number
* \param uplinkMax the highest uplink number
* \param nAdcs the number of adcs per uplink
* \param newFileName The name of the file that is created
*/
void correctFile(TTree* tree2correct, const map_uint_map_uint_double& gains,
                 const map_uint_map_uint_double& pedestals,
                 const unsigned int uplinkMin, const unsigned int uplinkMax,
                 const unsigned int nAdcs, TString newFileName);

/** \brief Removes the path from a filename
*
* \param str file name with a path
* \return file name without path
*/
TString removePath(TString str);

#endif
