void countEvents() {
      ROOT::EnableImplicitMT();
    // Open first file
    TFile f1("/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/DSTs/standalone/dfSelected_afterFid_afterCorr.root");
    if (f1.IsZombie()) {
        std::cerr << "Error opening dfSelected_afterFid_afterCorr.root" << std::endl;
        return;
    }
    auto t1 = (TTree*)f1.Get("dfSelected_afterFid_afterCorr");  // replace "myTree" with your tree name
    if (!t1) {
        std::cerr << "Tree not found in dfSelected_afterFid_afterCorr.root" << std::endl;
        return;
    }

    // Open second file
    TFile f2("/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/DSTs/script/dfSelected_afterFid_afterCorr.root");
    if (f2.IsZombie()) {
        std::cerr << "Error opening dfSelected_afterFid_afterCorr.root" << std::endl;
        return;
    }
    auto t2 = (TTree*)f2.Get("dfSelected_afterFid_afterCorr");  // replace "myTree" with your tree name
    if (!t2) {
        std::cerr << "Tree not found in dfSelected_afterFid_afterCorr.root" << std::endl;
        return;
    }

    // Count entries
    double n1 = t1->GetEntries();
    double n2 = t2->GetEntries();
    double total = n1 + n2;
    double ratio = (n2 * 100) / n1; // percentage

    std::cout << "Events in thread 1dfSelected_afterFid_afterCorr.root = " << n1 << std::endl;
    std::cout << "Events in from run_disana script dfSelected_afterFid_afterCorr.root = " << n2 << std::endl;
    std::cout << "ratio (percentage) = " << ratio << std::endl;
}
