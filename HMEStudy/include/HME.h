#ifndef HHbbWW_HME_H
#define HHbbWW_HME_H

#include <numeric>
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"

class TRandom3;

namespace hme {

  using LorentzVectorF = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;

  class HMEEvaluator {
  public:
    HMEEvaluator();

    std::pair<float, float> runHME(const LorentzVectorF l1,
                                   const LorentzVectorF l2,
                                   const LorentzVectorF b1,
                                   const LorentzVectorF b2,
                                   const LorentzVectorF met,
                                   const std::uint32_t eventnr,
                                   bool verbose = false) const;

    struct Sampler {
      std::vector<float> xval;
      std::vector<float> yval;
      std::vector<float> cdf;
      void computeCDF();
      float getRandom(TRandom3& rg) const;
    };

  private:
    Sampler m_mw_gen, m_brescale_gen;
  };
};  // namespace hme

#endif
