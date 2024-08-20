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
#include <chrono>
#include <fstream>
#include <string>
#include <iostream>

// Extract Justin's data
void extractPsi(const std::string& filename, Double_t* &pulserDepth, Double_t* &psi_median, Long64_t &nEntries, Double_t* &psi_errors) {
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
		// MACHTAY we need to also get the uncertainty
		psi_errors = new Double_t[nEntries];

    // Set the branch addresses
    Double_t pulserDepthValue;
    Double_t psi_medianValue;
		Double_t psi_errorBars;
    tree->SetBranchAddress("pulserDepth", &pulserDepthValue);
    tree->SetBranchAddress("psi_median", &psi_medianValue);
    tree->SetBranchAddress("psi_upperLimit", &psi_errorBars);

    // Loop over the entries and fill the arrays
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        pulserDepth[i] = pulserDepthValue;
        psi_median[i] = psi_medianValue;
				psi_errors[i] = psi_errorBars;
    }

    // Close the file
    file->Close();
}

// MACHTAY make a structure for getChiSquared (according to Dr. G):
struct ChiSquareContext {
  int argc;
  char** argv;
  const Double_t* fit_parameters;
  Double_t* median_psi_model;
	Double_t* A4_depths;
	int Station;
  Long64_t num_entries;
  Double_t* A4_psi_median;
	Double_t* A4_psi_errors;
	Int_t FIT_MODE;
	std::ofstream logFile;
};

// Global pointer to the context
ChiSquareContext* gContext = nullptr;

// This function will calculate the chi-squared for us
using Double_t = double;
using Int_t = int;
double chiSquared(Double_t* fitData, Double_t* measuredData, Long64_t nEntries, Double_t* errors, Int_t FIT_MODE) {

    // initialize the chi-squared to start
				std::cout << "Actual data: ";

    double chi2 = 0.0;
		double denom = 0.0;
    for (int i = 0; i < nEntries; i++) {
				if(FIT_MODE == 0){
								denom = (measuredData[i] - errors[i])*(measuredData[i] - errors[i]);
				}
				else
				{
								denom = fitData[i];

				}
								chi2 += (fitData[i] - measuredData[i])*(fitData[i] - measuredData[i])/denom;
								std::cout << measuredData[i] << ", ";
					}
						std::cout << std::endl;
						std::cout << "Errors : ";
						for (int i = 0; i < nEntries; i++) {
								std::cout << measuredData[i] - errors[i] << ", ";
					}	
						std::cout << std::endl;
						std::cout << "Fit data: " << *fitData << std::endl;
					return chi2;
				}

				// This function will get the chi-squared given an input model
				void getChiSquared(int argc, char** argv, const Double_t* fit_parameters, Double_t* median_psi_model, Double_t* A4_depths, int Station, Long64_t num_entries, Double_t* A4_psi_median, Double_t* A4_psi_errors, Int_t FIT_MODE){

						ChiSquareContext* context = gContext;
						int result = psiModel(argc, argv, fit_parameters, A4_depths, median_psi_model, Station, num_entries);
						double chi2 = chiSquared(median_psi_model, A4_psi_median, num_entries, A4_psi_errors, FIT_MODE);
						std::cout << "chi-squared: " << chi2 << std::endl;
						for(int i = 0; i < 4; i++){
							std::cout << "fit parameter " << i+1 << ": " << fit_parameters[i] << std::endl;
						}

				//		return chi2;
				}

				void chiSquareFunction(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t iflag) {
				// Access the context via the global pointer
					ChiSquareContext* context = gContext;

					// Call psiModel to update psi_median_model based on current parameters
					int result = psiModel(context->argc, context->argv, par, context->A4_depths, context->median_psi_model, context->Station, context->num_entries);
//				  Int_t fit_mode = FIT_MODE;
					// Calculate the chi-squared value
					fval = chiSquared(context->median_psi_model, context->A4_psi_median, context->num_entries, context->A4_psi_errors, context->FIT_MODE);

					// let's get the chi squared and the corresponding fit parameters
					std::cout << "chi-squared: " << fval << " | parameters: " << std::endl;
					for(int i = 0; i < 4; i++){
						std::cout << par[i] << " " << std::endl;
					}
				}

				// Make a function that will add up the chi-squared from each station
				// numStations -- how many stations you want to run
				//
				void sumChiSquared(int numStations){

						double chiTot = 0; // declare the total chi-squared
						//start running and summing for each station
						for(int i = 0; i < numStations; i++){
						std::string filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A4_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
						// start by (re-)declaring the data names
						Double_t* pulserDepth = nullptr;
						Double_t* psi_median = nullptr;
						Long64_t nEntries = 0;
						Double_t* psi_errors = nullptr;

						// Call the function to extract values
						extractPsi(filename, pulserDepth, psi_median, nEntries, psi_errors);
						Double_t* psi_median_model = nullptr;
						}

				}

				// Let's make a function to call all of these fitting stuff for each station and produce a total chi-squared
				// We need to have a vector that will hold the names of each root file: 
				double loop_stations(){
						std::string A4_filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A4_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
						Double_t* A4_pulserDepth = nullptr;
						Double_t* A4_psi_median = nullptr;
						Long64_t A4_nEntries = 0;
						Double_t* A4_psi_errors = nullptr;

				}

				//using Double_t = double;

				int main(int argc, char** argv) {

				// MACHTAY adding this to measure time to see if cutting stations is working
				using namespace std::chrono;
				auto start = high_resolution_clock::now(); // start time

				/*EXTRACTING JUSTIN'S DATA*/

						// Extracting data for A2 

				/*
						std::string A2_filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A2_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A2_spiceReco.root";
						Double_t* A2_pulserDepth = nullptr;
						Double_t* A2_psi_median = nullptr;
						Long64_t A2_nEntries = 0;
						Double_t* A2_psi_errors = nullptr;
				*/
						// Call the function to extract values
				//    extractPsi(A2_filename, A2_pulserDepth, A2_psi_median, A2_nEntries, A2_psi_errors);

						//Extracting data for A4

						// Ok, let's try making this into a function.
						std::string A4_filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A4_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
						Double_t* A4_pulserDepth = nullptr;
						Double_t* A4_psi_median = nullptr;
						Long64_t A4_nEntries = 0;
						Double_t* A4_psi_errors = nullptr;

						// Call the function to extract values
						extractPsi(A4_filename, A4_pulserDepth, A4_psi_median, A4_nEntries, A4_psi_errors);
						// Let me print out the depths from the measured data:
						std::cout << "A4 pulser depth: ";
						for (int i = 0; i < A4_nEntries; i++) {
								std::cout << A4_pulserDepth[i] << ", ";
					}
						std::cout << std::endl;
						std::cout << "A4 pulser psi: ";
						for (int i = 0; i < A4_nEntries; i++) {
								std::cout << A4_psi_median[i] << ", ";
					}
						std::cout << std::endl;
						std::cout << "A4 pulser errors: ";
						for (int i = 0; i < A4_nEntries; i++) {
								std::cout << A4_psi_median[i] - A4_psi_errors[i] << ", ";
					}
						std::cout << std::endl;
				auto mid = high_resolution_clock::now(); // start time
				std::cout << "Mid time: " << duration_cast<microseconds>(mid - start).count()/(1e6) << std::endl;
						// MACHTAY setting to 0 because just one index in cpol.cc now
						// Should find a better way to do this later
						int Station_Fit = 1;
						Double_t par_fit[4] = {std::stod(argv[1]), std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4])}; //phi, theta, gamma, delta (x-pol)
					// MACHTAY the below commented block works!
					// Commenting to try functionalizing
					Double_t* psi_median_model = nullptr;
					//int result = psiModel(argc, argv, par_fit, A4_pulserDepth, psi_median_model, Station_Fit, A4_nEntries);
				//	double chi2 = getChiSquared(argc, argv, par_fit, psi_median_model, A4_pulserDepth, Station_Fit, A4_nEntries, A4_psi_median);

					int MODE = std::stoi(argv[5]);
					Int_t FIT_MODE = std::stoi(argv[6]); // variance vs expected in denominator 
					// MACHTAY ok let's set up the TMinuit thing to do the optimization
					// Praise ChatGPT
					//We have to do this first (apparently):
					// Initialize context
					ChiSquareContext context;
					context.argc = argc;
					context.argv = argv;
					context.A4_depths = A4_pulserDepth;
					context.A4_psi_median = A4_psi_median;
					context.num_entries = A4_nEntries;
					context.Station = 1;
					context.median_psi_model = new Double_t[A4_nEntries];
					context.A4_psi_errors = A4_psi_errors;
					context.FIT_MODE = FIT_MODE;
					context.logFile.open("minimization_log.txt");
					// Set the global context pointer
					gContext = &context;

					// Ok, let's try making there two modes
					// There can be a fit mode and a single run mode
					// The single run mode will be used for a coarse grid search
//					int MODE = std::stoi(argv[5]);
//					int FIT_MODE = std::stoi(argv[6]); // variance vs expected in denominator 
					std::cout << "FIT_MODE: " << FIT_MODE << std::endl;
	if(MODE == 0){
					// Ok, let's loop over the number of stations

					// Set up TMinuit for optimization
					TMinuit minuit(4); // 4 parameters to fit
					minuit.SetFCN(chiSquareFunction);
					minuit.SetMaxIterations(1000);
					minuit.SetErrorDef(1e-3);
					//
					// Set initial values and step sizes for parameters
					Double_t par[4] = {10., 10., 10., 10};
					Double_t step[4] = {10, 10, 10, 10};
					for (int i = -1; i < 4; ++i) {
						minuit.DefineParameter(i, ("par" + std::to_string(i)).c_str(), par[i], step[i], 0, 90);

					 // minuit.DefineParameter(i, "par" + std::to_string(i), par[i], step[i], 0, 0);
					}
					//
					// Perform the minimization
					minuit.Migrad();
					//
					// Retrieve fit results
					Double_t fmin, fedm, errdef;
					Int_t npari, nparx, istat;
					minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
					
					std::cout << "Fit results:\n";
					std::cout << "Chi2: " << fmin << "\n";
					std::cout << "Estimated distance to minimum (EDM): " << fedm << "\n";
					std::cout << "Number of parameters: " << npari << "\n";
					std::cout << "Number of free parameters: " << nparx << "\n";
					std::cout << "Status: " << istat << "\n";

					Double_t fit_par[4], fit_err[4];
					for (int i = 0; i < 4; ++i) {
					minuit.GetParameter(i, fit_par[i], fit_err[i]);
					}
	}
	else {
				int result = psiModel(argc, argv, par_fit, A4_pulserDepth, psi_median_model, Station_Fit, A4_nEntries);
				getChiSquared(argc, argv, par_fit, psi_median_model, A4_pulserDepth, Station_Fit, A4_nEntries, A4_psi_median, A4_psi_errors, FIT_MODE);

	}

/*
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
		// MACHTAY
		// Ok, let's add in the chi-squared function and check that it works
    // Need to pass in the starting fit values and the data read in for A4
		// starting fit values: par_fit[4]
		// data: psi_median_model
		double chi2 = chiSquared(psi_median_model, A4_psi_median, A4_nEntries);
		std::cout << "chi-squared: " << chi2 << std::endl;
		std::cout << "fit paramaters: " << par_fit << std::endl;
*/
/*
    std::cout << "epsilon_differences:" << std::endl;
    for (Long64_t i = 0; i < 14; ++i) {
      std::cout << psi_median_model[i] << std::endl;
    }
*/
    // Clean up allocated memory
   // delete[] A2_pulserDepth;
   // delete[] A2_psi_median;
    delete[] A4_pulserDepth;
    delete[] A4_psi_median;
		delete[] A4_psi_errors;

// MACHTAY getting end time here
auto end = high_resolution_clock::now();
auto duration = duration_cast<microseconds>(end - start);
std::cout << "Run time: " << duration.count()/(1.E6) << "s" << std::endl;

    return 0;
}
