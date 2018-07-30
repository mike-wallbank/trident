void compile_trident_mumu(int nevents = 0,
			  bool run_signal = true,
			  int run_background = -1,
			  TString inFilePath = "/pnfs/dune/persistent/users/jmalbos/Trident/data/sim/mumu/",
			  TString outFile = "TridentMuMuOut.root") {
    gInterpreter->AddIncludePath(".");
    gSystem->Load("Particle_cxx");
    gSystem->Load("dunend/build/src/libEvent");
    gROOT->ProcessLine(".L trident_mumu.C+");
    if (nevents != 0)
      gROOT->ProcessLine(Form("trident_mumu(%d, 0, %d, %d, \"%s\", \"%s\")",
			      nevents, run_signal, run_background, inFilePath.Data(), outFile.Data()));
    else
      std::cout << "Compiled trident_mumu.C.  Run using trident_mumu()" << std::endl;
}
