// fit_custom_function.cc
#include <TROOT.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMinuit.h>
#include <iostream>
#include <algorithm>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "cpol.h"

// Extract Justin's data
void extractPsi(const std::string& filename, Double_t* &pulserDepth, Double_t* &psi_median, Long64_t &nEntries) {
    // Open the ROOT file
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get the tree from the file
    TTree *tree = (TTree*)file->Get("polReco");
    if (!tree) {
        std::cerr << "Error getting tree from file!" << std::endl;
        file->Close();
        return;
    }

    // Get the number of entries
    nEntries = tree->GetEntries();
    if (nEntries <= 0) {
        std::cerr << "No entries found in the tree!" << std::endl;
        file->Close();
        return;
    }

    // Allocate arrays for pulserDepth and psi_median
    pulserDepth = new Double_t[nEntries];
    psi_median = new Double_t[nEntries];

    // Set the branch addresses
    Double_t pulserDepthValue;
    Double_t psi_medianValue;
    tree->SetBranchAddress("pulserDepth", &pulserDepthValue);
    tree->SetBranchAddress("psi_median", &psi_medianValue);

    // Loop over the entries and fill the arrays
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        pulserDepth[i] = pulserDepthValue;
        psi_median[i] = psi_medianValue;
    }

    // Close the file
    file->Close();
}

// Declare the function prototype for psiModel
//int psiModel(int argc, char** argv, const Double_t* par_fit, Double_t* &fit_depths, Double_t* &psi_median_model, int Station_Fit);

/*
// Declare the custom function
void customFunction(const Double_t* x, const Double_t* par, Double_t* output); // ASG: This should be the main function from cpol.cc. Make it stop and return the epsilon values 
                                                                               // when they've been calculated.
                                                                               // ETA: A GOOD(!) day, two meh days

const int num_outputs = 5; // Two curves for now, five later

// Data points. ASG: Import Justin's data here. ETA: 2 hours?
Double_t x_data[] = {-2, -1, 0, 1, 2, 3, 4}; //ASG: Depths
Double_t y_data[num_outputs][7] = { //ASG: Polarization angles
    {4.1, 1.2, 0.9, 2.2, 4.5, 9.1, 15.6},
    {3.9, 1.1, 1.0, 2.0, 4.4, 9.2, 15.4},
    {4.0, 1.3, 0.8, 2.1, 4.3, 9.3, 15.5},
    {4.2, 1.0, 1.1, 2.3, 4.6, 9.0, 15.7},
    {4.3, 1.2, 0.7, 2.4, 4.7, 9.4, 15.8}
}; 

// Function to be minimized by TMinuit
void chiSquareFunction(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t iflag) {
    Double_t chi2 = 0.0;
    for (int j = 0; j < num_outputs; ++j) {
        for (int i = 0; i < 7; ++i) {
            Double_t x[1] = {x_data[i]}; // ASG: modifications required here. Calculate the curve once per set of depths. 
            Double_t y[num_outputs];    // ETA: 1-2 Hours?
            customFunction(x, par, y); // x not necessary
            chi2 += (y[j] - y_data[j][i]) * (y[j] - y_data[j][i]);
        }
    }
    fval = chi2;
}

void fit_custom_function() {
    TMinuit minuit(4); // 4 parameters to fit
    minuit.SetFCN(chiSquareFunction);

    // Set initial values and step sizes for parameters
    Double_t par[4] = {1.0, 1.0, 1.0, 1.0};
    Double_t step[4] = {0.1, 0.1, 0.1, 0.1}; //ASG: Resolution too high?
    for (int i = 0; i < 4; ++i) {
        minuit.DefineParameter(i, "par" + std::to_string(i), par[i], step[i], 0, 0);
    }

    // Perform the minimization
    minuit.Migrad();

    //ASG: The rest is just retrieving results, just learn details of what it is doing

    // Retrieve fit results
    Double_t fit_par[4], fit_err[4];
    for (int i = 0; i < 4; ++i) {
        minuit.GetParameter(i, fit_par[i], fit_err[i]);
    }

    // Print fit results
    std::cout << "Fit results:\n";
    std::cout << "Chi2: " << minuit.fAmin << "\n";
    std::cout << "Parameters:\n";
    for (int i = 0; i < 4; ++i) {
        std::cout << "par[" << i << "] = " << fit_par[i] << " +/- " << fit_err[i] << "\n";
    }

    // Create a canvas to draw on
    TCanvas* c1 = new TCanvas("c1", "Fit Example with Custom Function", 800, 600);
    
    // Draw the graph for the first output set
    TGraph* graph = new TGraph(7, x_data, y_data[0]);
    graph->Draw("AP"); // "A" to draw axis, "P" to draw points

    // Optionally, add more graphs for each output set to the canvas
    for (int j = 1; j < num_outputs; ++j) {
        TGraph* graph_j = new TGraph(7, x_data, y_data[j]);
        graph_j->SetMarkerColor(j+1); // Change color for each set
        graph_j->Draw("P same");
    }

    // Save the canvas as an image
    c1->SaveAs("fit_custom_function.png");

    // Show the canvas (optional, depends on your environment)
    c1->Draw();
}

// Execute the function
void run_fit_custom_function() {
    fit_custom_function();
}
*/

using Double_t = double;

int main(int argc, char** argv) {

/*EXTRACTING JUSTIN'S DATA*/

    // Extracting data for A2 

    std::string A2_filename = "/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A2_spiceReco.root";
    Double_t* A2_pulserDepth = nullptr;
    Double_t* A2_psi_median = nullptr;
    Long64_t A2_nEntries = 0;

    // Call the function to extract values
    extractPsi(A2_filename, A2_pulserDepth, A2_psi_median, A2_nEntries);

    //Extracting data for A4

    std::string A4_filename = "/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
    Double_t* A4_pulserDepth = nullptr;
    Double_t* A4_psi_median = nullptr;
    Long64_t A4_nEntries = 0;

    // Call the function to extract values
    extractPsi(A4_filename, A4_pulserDepth, A4_psi_median, A4_nEntries);

    int Station_Fit = 4;
    Double_t par_fit[4] = {0.0, 10.0, 0.0, -5.0};

    std::cout << "TESTING WORKING ASG" << std::endl;
    Double_t* psi_median_model = nullptr;
    int result = psiModel(argc, argv, par_fit, A4_pulserDepth, psi_median_model, Station_Fit, A4_nEntries);

    if (result == 0) { // Check if psiModel was successful
        std::cout << "epsilon_differences:" << std::endl;
        for (Long64_t i = 0; i < A4_nEntries; ++i) {
	    std::cout << "ENTERING" << std::endl;
            std::cout << psi_median_model[i] << std::endl;
        }
    } else {
        std::cerr << "psiModel returned an error." << std::endl;
    }
/*
    std::cout << "epsilon_differences:" << std::endl;
    for (Long64_t i = 0; i < 14; ++i) {
      std::cout << psi_median_model[i] << std::endl;
    }
*/
    // Clean up allocated memory
    delete[] A2_pulserDepth;
    delete[] A2_psi_median;
    delete[] A4_pulserDepth;
    delete[] A4_psi_median;

    return 0;
}
