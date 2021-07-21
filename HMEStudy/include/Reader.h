#ifndef HHbbWW_READER_H
#define HHbbWW_READER_H

#include "TFile.h"
#include "TTree.h"
#include <csignal>
#include <iostream>

#include "HME.h"


namespace hme {
  class HMEReader{
  public:
    HMEReader(std::string filepath);
    ~HMEReader();

    float getHME(unsigned int run,
                 unsigned int lumi,
                 unsigned long long event,
                 const LorentzVectorF l1,
                 const LorentzVectorF l2,
                 const LorentzVectorF b1,
                 const LorentzVectorF b2,
                 const LorentzVectorF met,
                 bool boosted_tag);
                 
  private:
    TFile* m_file;
    TTree* m_tree;
    std::size_t m_entry = 0;
    std::size_t m_end;
    std::size_t m_recomputed = 0;
    std::size_t m_total = 0;
    std::string m_filepath;
    unsigned int m_run;
    unsigned int m_ls;
    unsigned long long m_event;
    float m_hme;
    bool valid_tree;
    HMEEvaluator m_hmeEval;
  };
};  // namespace hme

#endif
