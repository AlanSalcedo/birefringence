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

using namespace std;

// Extract Justin's data
void extractPsi(const string& filename, Double_t* &pulserDepth, Double_t* &psi_median, Long64_t &nEntries, Double_t* &psi_errors) {
    // Open the ROOT file
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Error opening file!" << endl;
        return;
    }

    // Get the tree from the file
    TTree *tree = (TTree*)file->Get("polReco");
    if (!tree) {
        cerr << "Error getting tree from file!" << endl;
        file->Close();
        return;
    }

    // Get the number of entries
    nEntries = tree->GetEntries();
    if (nEntries <= 0) {
        cerr << "No entries found in the tree!" << endl;
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
	ofstream logFile;
};

// Global pointer to the context
ChiSquareContext* gContext = nullptr;

// This function will calculate the chi-squared for us
using Double_t = double;
using Int_t = int;
double chiSquared(Double_t* fitData, Double_t* measuredData, Long64_t nEntries, Double_t* errors, Int_t FIT_MODE) {

    // initialize the chi-squared to start
				cout << "Actual data: ";

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
								cout << measuredData[i] << ", ";
					}
						cout << endl;
						cout << "Errors : ";
						for (int i = 0; i < nEntries; i++) {
								cout << measuredData[i] - errors[i] << ", ";
					}	
						cout << endl;
						cout << "Fit data: " << *fitData << endl;
		return chi2;
}

// This function will get the chi-squared given an input model
void getChiSquared(int argc, char** argv, const Double_t* fit_parameters, Double_t* median_psi_model, Double_t* A4_depths, int Station, Long64_t num_entries, Double_t* A4_psi_median, Double_t* A4_psi_errors, Int_t FIT_MODE){

		ChiSquareContext* context = gContext;
		int result = psiModel(argc, argv, fit_parameters, A4_depths, median_psi_model, Station, num_entries);
		double chi2 = chiSquared(median_psi_model, A4_psi_median, num_entries, A4_psi_errors, FIT_MODE);
		cout << "chi-squared: " << chi2 << endl;
		for(int i = 0; i < 4; i++){
			cout << "fit parameter " << i+1 << ": " << fit_parameters[i] << endl;
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
	cout << "chi-squared: " << fval << " | parameters: " << endl;
	for(int i = 0; i < 4; i++){
		cout << par[i] << " " << endl;
	}
}

// Make a function that will add up the chi-squared from each station
// numStations -- how many stations you want to run
//
void sumChiSquared(int numStations){

		double chiTot = 0; // declare the total chi-squared
		//start running and summing for each station
		for(int i = 0; i < numStations; i++){
		string filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A4_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
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
				// We need to have vectors for each of the data features: 
				vector<string> filenames = {"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root"}; //{ "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A2_spiceReco.root", "/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root"};
				vector<Double_t*> pulserDepths(filenames.size(), nullptr);
				vector<Double_t*> psi_medians(filenames.size(), nullptr);
				vector<Long64_t> nEntries(filenames.size(), 0);
				vector<Double_t*> psi_errors(filenames.size(), nullptr);
				double loop_stations(){
						string filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A4_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
						Double_t* pulserDepth = nullptr;
						Double_t* psi_median = nullptr;
						Long64_t nEntries = 0;
						Double_t* psi_errors = nullptr;

				}

				//using Double_t = double;

int main(int argc, char** argv) {
				int NUM_STATIONS = stoi(argv[7]);
				// MACHTAY adding this to measure time to see if cutting stations is working
				using namespace chrono;
				auto start = high_resolution_clock::now(); // start time

				/*EXTRACTING JUSTIN'S DATA*/

						//Extracting data for A4
						// Ok, let's try making this into a function.
						string filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A4_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
						Double_t* pulserDepth = nullptr;
						Double_t* psi_median = nullptr;
						Long64_t nEntries = 0;
						Double_t* psi_errors = nullptr;
  					Double_t* psi_median_model = nullptr;
						// Make vectors instead
            // We need to have vectors for each of the data features: 
				    vector<string> filenames = {"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A2_spiceReco.root", "/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root"}; //{ "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A2_spiceReco.root", "/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root"};
				    vector<Double_t*> pulserDepths(filenames.size(), nullptr);
		    		vector<Double_t*> psi_medians(filenames.size(), nullptr);
			    	vector<Long64_t> numEntries(filenames.size(), 0);
			    	vector<Double_t*> Psi_errors(filenames.size(), nullptr);
				    vector<Double_t*> psi_median_models(filenames.size(), nullptr);
						// Call the function to extract values
						extractPsi(filename, pulserDepth, psi_median, nEntries, psi_errors);
						// Now try extracting using vectors!
						for(int i = 0; i < NUM_STATIONS; i++){
								extractPsi(filenames[i], pulserDepths[i], psi_medians[i], numEntries[i], Psi_errors[i]);
            }
						// Let me print out the depths from the measured data:
						cout << "A4 pulser depth: ";
						for (int i = 0; i < nEntries; i++) {
								cout << pulserDepth[i] << ", ";
					}
						cout << endl;
						cout << "A4 pulser psi: ";
						for (int i = 0; i < nEntries; i++) {
								cout << psi_median[i] << ", ";
					}
						cout << endl;
						cout << "A4 pulser errors: ";
						for (int i = 0; i < nEntries; i++) {
								cout << psi_median[i] - psi_errors[i] << ", ";
					}
						cout << endl;
				auto mid = high_resolution_clock::now(); // start time
				cout << "Mid time: " << duration_cast<microseconds>(mid - start).count()/(1e6) << endl;
						// MACHTAY setting to 0 because just one index in cpol.cc now
						// Should find a better way to do this later
						int Station_Fit = 1;
						vector<int> Station_Fits(filenames.size(), 1); 
						Double_t par_fit[4] = {stod(argv[1]), stod(argv[2]), stod(argv[3]), stod(argv[4])}; //phi, theta, gamma, delta (x-pol)
//						vector<Double_t> par_fits(filenames.size(), par_fit);
					// MACHTAY the below commented block works!
					// Commenting to try functionalizing
					//int result = psiModel(argc, argv, par_fit, A4_pulserDepth, psi_median_model, Station_Fit, A4_nEntries);
				//	double chi2 = getChiSquared(argc, argv, par_fit, psi_median_model, A4_pulserDepth, Station_Fit, A4_nEntries, A4_psi_median);

					int MODE = stoi(argv[5]);
					Int_t FIT_MODE = stoi(argv[6]); // variance vs expected in denominator 
					// MACHTAY ok let's set up the TMinuit thing to do the optimization
					// Praise ChatGPT
					//We have to do this first (apparently):
					// Initialize context
					ChiSquareContext context;
					context.argc = argc;
					context.argv = argv;
					context.A4_depths = pulserDepth;
					context.A4_psi_median = psi_median;
					context.num_entries = nEntries;
					context.Station = 1;
					context.median_psi_model = new Double_t[nEntries];
					context.A4_psi_errors = psi_errors;
					context.FIT_MODE = FIT_MODE;
					context.logFile.open("minimization_log.txt");
					// Set the global context pointer
					gContext = &context;

					// Ok, let's try making there two modes
					// There can be a fit mode and a single run mode
					// The single run mode will be used for a coarse grid search
//					int MODE = stoi(argv[5]);
//					int FIT_MODE = stoi(argv[6]); // variance vs expected in denominator 
					cout << "FIT_MODE: " << FIT_MODE << endl;
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
						minuit.DefineParameter(i, ("par" + to_string(i)).c_str(), par[i], step[i], 0, 90);

					 // minuit.DefineParameter(i, "par" + to_string(i), par[i], step[i], 0, 0);
					}
					//
					// Perform the minimization
					minuit.Migrad();
					//
					// Retrieve fit results
					Double_t fmin, fedm, errdef;
					Int_t npari, nparx, istat;
					minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
					
					cout << "Fit results:\n";
					cout << "Chi2: " << fmin << "\n";
					cout << "Estimated distance to minimum (EDM): " << fedm << "\n";
					cout << "Number of parameters: " << npari << "\n";
					cout << "Number of free parameters: " << nparx << "\n";
					cout << "Status: " << istat << "\n";

					Double_t fit_par[4], fit_err[4];
					for (int i = 0; i < 4; ++i) {
					minuit.GetParameter(i, fit_par[i], fit_err[i]);
					}
	}
	else {
//				int result = psiModel(argc, argv, par_fit, pulserDepth, psi_median_model, Station_Fit, nEntries);
//				getChiSquared(argc, argv, par_fit, psi_median_model, pulserDepth, Station_Fit, nEntries, psi_median, psi_errors, FIT_MODE);
				// Let's try looping over multiple stations!
				cout << "In else statement!" << endl;
				for(int i = 0; i < NUM_STATIONS; i++){
				int station = Station_Fits[i] + i;
				int result = psiModel(argc, argv, par_fit, pulserDepths[i], psi_median_models[i], station, numEntries[i]);
				getChiSquared(argc, argv, par_fit, psi_median_model, pulserDepth, Station_Fit, nEntries, psi_median, psi_errors, FIT_MODE);
				}
	}

/*
    cout << "TESTING WORKING ASG" << endl;
    Double_t* psi_median_model = nullptr;
    int result = psiModel(argc, argv, par_fit, A4_pulserDepth, psi_median_model, Station_Fit, A4_nEntries);




    if (result == 0) { // Check if psiModel was successful
        cout << "epsilon_differences:" << endl;
        for (Long64_t i = 0; i < A4_nEntries; ++i) {
	    cout << "ENTERING" << endl;
            cout << psi_median_model[i] << endl;
        }
    } else {
        cerr << "psiModel returned an error." << endl;
    }
		// MACHTAY
		// Ok, let's add in the chi-squared function and check that it works
    // Need to pass in the starting fit values and the data read in for A4
		// starting fit values: par_fit[4]
		// data: psi_median_model
		double chi2 = chiSquared(psi_median_model, A4_psi_median, A4_nEntries);
		cout << "chi-squared: " << chi2 << endl;
		cout << "fit paramaters: " << par_fit << endl;
*/
/*
    cout << "epsilon_differences:" << endl;
    for (Long64_t i = 0; i < 14; ++i) {
      cout << psi_median_model[i] << endl;
    }
*/
    // Clean up allocated memory
   // delete[] A2_pulserDepth;
   // delete[] A2_psi_median;
    delete[] pulserDepth;
    delete[] psi_median;
		delete[] psi_errors;

// MACHTAY getting end time here
auto end = high_resolution_clock::now();
auto duration = duration_cast<microseconds>(end - start);
cout << "Run time: " << duration.count()/(1.E6) << "s" << endl;

    return 0;
}
