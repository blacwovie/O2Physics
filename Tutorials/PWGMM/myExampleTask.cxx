// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author everyone

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/Logger.h"
#include <cmath> // for std::isnan
using namespace o2;
using namespace o2::framework;

using namespace o2::aod;

struct myExampleTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};//HistogramRegistry O2包装好的类，帮忙定义输出
//or HistogramRegistry histos{"histos"} 
//和上一行效果一样
//如果要生成非直方图的函数⬇️
//OutputObj<YourClass> fName{YourClass("name")};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsDedx{"nBinsDedx", 2000, "N bins in dE/dx histo"};
  Configurable<int> nBinsnsigma{"nBinsnsigma", 10, "N bins in nsigma histo"};
  Configurable<int> He3nBinsnsigma{"He3nBinsnsigma", 10, "N bins in nsigma histo"};
  Configurable<int> He3nBinsPt{"He3nBinsPt", 100, "N bins in nsigma histo"};
  void init(InitContext const&)//init 初始值的分配
  {
    // define axes you want to use
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, +6, "pt (GeV/c)"};
    const AxisSpec axisDedx{nBinsDedx, 0., 1400., "dE/dx (a.u.)"};
    const AxisSpec axisnsigma{nBinsnsigma, -5., +5., ""};
    const AxisSpec axisHePt{He3nBinsPt, 0, +10, "pt (GeV/c)"};
    const AxisSpec axisHensigma{He3nBinsnsigma, -10., +10., ""};
    

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    histos.add("hEventCounter", "hEventCounter", kTH1F, {{1, 0., 1.}});//1个bin，0-1

    histos.add("dedxVsPtHistogram", "dE/dx vs p", kTH2F, {{axisPt}, {axisDedx}}); 
    histos.add("He3nsigmaVsPtHistogram", "He3nsigma vs pt", kTH2F, {{axisHensigma}, {axisHePt}}); 
    histos.add("PrnsigmaVsPtHistogram", "Prnsigma vs pt", kTH2F, {{axisPt},{axisnsigma}});
  }

  // void process(aod::TracksIU const& tracks) //TracksIU???
  //void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& tracks)//collision无需for循环，默认对每个collision计数
  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPr, aod::pidTPCHe> const& tracks)
  //void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& tracks, aod::pidTPCPr const& nsigmaP, aod::pidTPCHe const& nsigmaHe)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {//遍历track table的每一行
      //if( track.tpcNClsCrossedRows() < 70) continue;//badly tracked
      //if( fabs(track.dcaXY()) > 0.2) continue;//doesn't point to primary vertex
      //if( track.pt() < 0.2) continue;
     //if( fabs(track.eta()) > 0.9) continue;

      histos.fill(HIST("etaHistogram"), track.eta());//填入每个track的eta
      histos.fill(HIST("ptHistogram"), track.pt());
      //histos.fill(HIST("He3nsigmaVsPtHistogram"), track.tpcNSigmaStoreHe(), track.pt());
      //histos.fill(HIST("PrnsigmaVsPtHistogram"), track.pt(), track.tpcNSigmaStorePr());
      histos.fill(HIST("dedxVsPtHistogram"), track.pt(),track.tpcSignal());
      if (!std::isnan(track.tpcNSigmaStoreHe())){
      //LOG(INFO) << "TPC Nsigma for He: " << track.tpcNSigmaStoreHe();
      //LOG(INFO) << "pt for He: " << track.pt();
      histos.fill(HIST("He3nsigmaVsPtHistogram"), track.tpcNSigmaStoreHe(), track.pt());
      }
      if (!std::isnan(track.tpcNSigmaStorePr())){
        histos.fill(HIST("PrnsigmaVsPtHistogram"), track.pt(), track.tpcNSigmaStorePr());
      }

    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)};
}
