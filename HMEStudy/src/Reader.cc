#include "Reader.h"

#ifdef HME_DEBUG
#define LogDebug std::cout
#else
#define LogDebug \
  if (false)     \
  std::cout
#endif

hme::HMEReader::HMEReader(std::string filepath) : m_filepath(filepath) {
    m_file = new TFile(m_filepath.c_str());
    valid_tree = m_file->GetListOfKeys()->Contains("Events");   
    std::cout<<"hme::HMEReader : Loading"<<m_filepath;
    if (valid_tree){    
        m_tree = (TTree*)m_file->Get("Events");
        m_tree->SetBranchAddress("run", &m_run);
        m_tree->SetBranchAddress("ls", &m_ls);
        m_tree->SetBranchAddress("event", &m_event);
        m_tree->SetBranchAddress("HME", &m_hme);
        m_tree->GetEntry(0);
        m_end = m_tree->GetEntries();
        std::cout<<" -> entries : "<<m_tree->GetEntries();
    }
    else{
        std::cout<<" -> empty";
    }
    std::cout<<std::endl;
    m_hmeEval = hme::HMEEvaluator();
}


hme::HMEReader::~HMEReader(){
    std::cout<<"hme::HMEReader : File "<<m_filepath<<" -> Number of recomputed HME values = "<<m_recomputed<<"/"<<m_total<<" ["<<float(m_recomputed)*100/std::max(1,static_cast<int>(m_total))<<"%]"<<std::endl;
}

float hme::HMEReader::getHME(unsigned int run,
                             unsigned int ls,
                             unsigned long long event,
                             const hme::LorentzVectorF l1,
                             const hme::LorentzVectorF l2,
                             const hme::LorentzVectorF b1,
                             const hme::LorentzVectorF b2,
                             const hme::LorentzVectorF met,
                             bool boosted_tag){
    if (!valid_tree){
        std::cout<<"hme::HMEReader : No tree in "<<m_filepath<<" but `getHME` is called, HME recomputed"<<std::endl;;
        std::pair<float,float> hme = m_hmeEval.runHME(l1,l2,b1,b2,met,static_cast<std::uint32_t>(event),boosted_tag);
        return hme.first;
    }
    bool wasReset = false;
    m_total ++;
    
    std::size_t  current_entry = m_entry;

    while (true) {
      if ( m_run == run && m_ls == ls && m_event == event ) {
        LogDebug<<"hme::HMEReader : Found match : run "<<run<<" ls "<<ls<<" event "<<event<<" -> HME = "<<m_hme<<std::endl;
        return m_hme;
      }
      ++m_entry;
      if ( m_entry == m_end ) {
        if (wasReset) {
          // recompute HME 
          m_recomputed++;
          std::pair<float,float> hme = m_hmeEval.runHME(l1,l2,b1,b2,met,static_cast<std::uint32_t>(event),boosted_tag);
          m_entry = current_entry; // return to where the pointer was when the event was asked
          return hme.first;
        } else {
          m_entry = 0;
          wasReset = true;
        }
      }
      m_tree->GetEntry(m_entry);
    }
    return -1; // avoid compiler warning/error

}
