// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoUtils.h
/// \brief Utilities for the FemtoDream framework
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMUTILS_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMUTILS_H_

#include <vector>
#include <string>
#include <functional>
#include <algorithm>
#include "Framework/ASoAHelpers.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGCF/DataModel/FemtoDerived.h"

namespace o2::analysis::femtoDream
{

// TODO: remove all these functions pertaining to PID selection for the next tutorial session they have been removed from femtodream tasks but are still present in tutorial files

enum kDetector { kTPC,
                 kTPCTOF,
                 kNdetectors };

/// internal function that returns the kPIDselection element corresponding to a
/// specifica n-sigma value \param nSigma number of sigmas for PID
/// \param vNsigma vector with the number of sigmas of interest
/// \return kPIDselection corresponding to n-sigma
int getPIDselection(float nSigma, std::vector<float> vNsigma)
{//返回 vNsigma 中与 nSigma 对应的索引。如果没有找到精确匹配的 n-sigma 值，则返回最接近的一个值
  std::sort(vNsigma.begin(), vNsigma.end(), std::greater<>());//比较器对向量中的 n-sigma 值进行降序排列。这样，较大的 n-sigma 值会在向量的前面。
  auto it = std::find(vNsigma.begin(), vNsigma.end(), nSigma);
  //在 vNsigma 中查找与 nSigma 完全匹配的值。如果找到了，it 将指向匹配的元素；否则，it 将是 vNsigma.end()，表示没有找到目标值
  if (it == vNsigma.end()) {
    it = vNsigma.begin() + 1;
    LOG(warn) << "Invalid value of nSigma: " << nSigma << ". Return the first value of the vector: " << *(it);
  }//如果 nSigma 不存在于 vNsigma 中，程序会说明输入的 nSigma 是无效的，并且将 it 指向 vNsigma 的第二个元素。
  return std::distance(vNsigma.begin(), it);
}//std::distance 计算 it 和 vNsigma.begin() 之间的距离，返回的是 it 的索引，即最接近 nSigma 值的位置。

/// function that checks whether the PID selection specified in the vectors is
/// fulfilled
/// \param pidcut Bit-wise container for the PID
/// \param vSpecies vector with ID corresponding to the selected species (output from cutculator)
/// \param nSpecies number of available selected species (output from cutculator)
/// \param nSigma number of sigma selection fo PID
/// \param vNsigma vector with available n-sigma selections for PID
/// \param kDetector enum corresponding to the PID technique
/// \return Whether the PID selection specified in the vectors is fulfilled
bool isPIDSelected(aod::femtodreamparticle::cutContainerType pidcut,
                   int vSpecies,
                   int nSpecies,
                   float nSigma,
                   std::vector<float> vNsigma,
                   kDetector iDet)//iDet：使用的检测器类型的枚举值kDetector，表示只用 TPC 检测器或同时用 TPC 和 TOF 进行粒子鉴别
{
  int iNsigma = getPIDselection(nSigma, vNsigma);
  int nDet = static_cast<int>(kDetector::kNdetectors);
  int bit_to_check = 1 + (vNsigma.size() - (iNsigma + 1)) * nDet * nSpecies + (nSpecies - (vSpecies + 1)) * nSpecies + (nDet - 1 - iDet);
  return ((pidcut >> (bit_to_check)) & 1) == 1;
  //如果结果为 1，说明该粒子满足 PID 选择条件，函数返回 true；如果为 0，返回 false。
};

/// function that checks whether the PID selection specified in the vectors is fulfilled, depending on the momentum TPC or TPC+TOF PID is conducted
/// \param pidcut Bit-wise container for the PID
/// \param momentum Momentum of the track
/// \param pidThresh Momentum threshold that separates between TPC and TPC+TOF PID
/// \param vSpecies Vector with the species of interest (number returned by the CutCulator)
/// \param nSpecies number of available selected species (output from cutculator)
/// \param nSigmaTPC Number of TPC sigmas for selection
/// \param nSigmaTPCTOF Number of TPC+TOF sigmas for selection (circular selection)
/// \return Whether the PID selection is fulfilled
bool isFullPIDSelected(aod::femtodreamparticle::cutContainerType const& pidCut,
                       float momentum,
                       float pidThresh,
                       int vSpecies,
                       int nSpecies,
                       std::vector<float> vNsigma,
                       float nSigmaTPC,
                       float nSigmaTPCTOF)
{
  bool pidSelection = true;
  if (momentum < pidThresh) {
    /// TPC PID only
    pidSelection = isPIDSelected(pidCut, vSpecies, nSpecies, nSigmaTPC, vNsigma, kDetector::kTPC);
  } else {
    /// TPC + TOF PID
    pidSelection = isPIDSelected(pidCut, vSpecies, nSpecies, nSigmaTPCTOF, vNsigma, kDetector::kTPCTOF);
  }
  return pidSelection;
};

/// function for getting the mass of a particle depending on the pdg code
/// \param pdgCode pdg code of the particle
/// \return mass of the particle
inline float getMass(int pdgCode)
{
  // use this function instead of TDatabasePDG to return masses defined in the PhysicsConstants.h header
  // this approach saves a lot of memory and important partilces like deuteron are missing in TDatabasePDG anyway
  //这个函数根据粒子的 PDG 编号返回粒子的质量。例如，PDG 码 211 代表正电介子（π+），因此函数会返回相应的质量
  float mass = 0;
  // add new particles if necessary here
  switch (std::abs(pdgCode)) {
    case kPiPlus: // charged pions, changed magic number as per their pdg name
      mass = o2::constants::physics::MassPiPlus;
      break;
    case kKPlus: // charged kaon
      mass = o2::constants::physics::MassKPlus;
      break;
    case kProton: // proton
      mass = o2::constants::physics::MassProton;
      break;
    case kLambda0: // Lambda
      mass = o2::constants::physics::MassLambda;
      break;
    case o2::constants::physics::Pdg::kPhi: // Phi Meson
      mass = o2::constants::physics::MassPhi;
      break;
    case o2::constants::physics::Pdg::kLambdaCPlus: // Charm Lambda
      mass = o2::constants::physics::MassLambdaCPlus;
      break;
    case o2::constants::physics::Pdg::kDeuteron: // Deuteron
      mass = o2::constants::physics::MassDeuteron;
      break;
    default:
      LOG(fatal) << "PDG code is not suppored";
  }
  return mass;
}

inline int checkDaughterType(o2::aod::femtodreamparticle::ParticleType partType, int motherPDG)
{
  int partOrigin = 0;
  if (partType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
    switch (abs(motherPDG)) {
      case 3122:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterLambda;
        break;
      case 3222:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterSigmaplus;
        break;
      default:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
    } // switch

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kV0) {
    partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kV0Child) {
    switch (abs(motherPDG)) {
      case 3122:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterLambda;
        break;
      case 3222:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterSigmaplus;
        break;
      default:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
    } // switch

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
    partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor) {
    partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
  }
  return partOrigin;
};

template <typename T, typename R>
inline bool containsNameValuePair(const std::vector<T>& myVector, const std::string& name, R value)
{
  for (const auto& obj : myVector) {
    if (obj.name == name) {
      if (std::abs(static_cast<float>((obj.defaultValue.template get<R>() - value))) < 1e-2) {
        return true; // Found a match
      }
    }
  }
  return false; // No match found
}

} // namespace o2::analysis::femtoDream
#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMUTILS_H_
