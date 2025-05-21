#include <string>
#include <iostream>
#include <chrono>
#include "RHipoDS.hxx"
#include "TCanvas.h"

#include "../source/ParticleInformation/RECParticle.h"
#include "../source/Cuts/ElectronCut.h"

using namespace ROOT;
using namespace ROOT::RDF;


int main(int argc, char **argv) {
   // Very simple test of the Hipo DataFrame.
   // ROOT::EnableImplicitMT();
   //int N_open = 100;
   std::chrono::nanoseconds delta_t;

   if(argc < 2){
      std::cout << "Please specify a HIPO data file on the command line. (Only one file.) \n";
      return 1;
   }else{
      std::cout << "Opening file " << argv[1] << std::endl;
   }

   auto start = std::chrono::high_resolution_clock::now();
   auto ds = std::make_unique<RHipoDS>(argv[1]);
   auto cols_ds = ds->GetColumnNames();
   bool translated = ds->fColumnNameTranslation;
   auto stop = std::chrono::high_resolution_clock::now();
   auto total_events = ds->GetEntries();
   delta_t = std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
   printf("Open file in  %6.5f ms  for %6lu events = %6.5f ns/event\n",
          delta_t.count()*1.e-6, total_events, double(delta_t.count())/total_events );

   //("/data/CLAS12/data/hipo/rec_clas_016321.evio.00001.hipo");
//   auto all_columns = ds->GetColumnNames();
//   for(int i=0; i< all_columns.size(); ++i){
//      printf("%40s  bank id: %4d  %s \n", all_columns[i].c_str(), i, ds->fColumnTypeIsVector[i] ? "vector":"scaler" );
//   }

   ds->fDebug = 0;
   auto df = RDataFrame(std::move(ds));
   auto df_selected = df.Filter(ElectronCut, RECParticle::All());
   auto h_electron_pz = df_selected
    .Define("electron_pz", get_RECParticle_float_var(11,-1), {"REC_Particle_phi", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
    .Histo1D({"h_electron_pz", "pz of only electron", 500, 0, 8}, "electron_pz");


   TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
   h_electron_pz->Draw();
   c3->SaveAs("pz_pid11.png");

   auto h_electron_pid = df_selected
    .Define("electron_pid", get_RECParticle_int_var(11), {"REC_Particle_pid", "REC_Particle_pid"})
    .Histo1D({"h_electron_pid", "pid of only electron", 500, -20, 20}, "electron_pid");


   TCanvas* c4 = new TCanvas("c4", "c4", 800, 600);
   h_electron_pid->Draw();
   c4->SaveAs("pid_pid11.png");

}