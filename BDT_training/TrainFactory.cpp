#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Classification.h"
#include "TMVA/CrossValidation.h"
#include <vector>
#include <iostream>

// Calculate the mean and standard deviation of a vector of doubles.
void AverageAndStdDeviation(const std::vector<double>& data, double& average, double& stdDeviation) {
  average = 0;
  stdDeviation = 0;
  int count = 0;
  for (auto& value : data) 
  {
   if (value > 0)
   {
      average += value;
      count += 1;
   }
   else std::cout<<"Skipping over index "<<count<<std::endl;
  }
  average /= count;
  count = 0;
  for (auto& value : data) 
  {
   if (value > 0)
   {
      stdDeviation += (value - average) * (value - average);
      count+= 1;
   }
  }
  stdDeviation /= (count-1);
  stdDeviation = sqrt(stdDeviation);
}

// This is the main function that trains the BDT and evaluates it using cross-validation.
void TrainFactory(UInt_t jobs = 1)
{
   TMVA::Tools::Instance();

   // Open the file with the input data
   TFile* input{nullptr};
   TString fname = "trainingV10.root";
   if (!gSystem->AccessPathName(fname)) {
      input = TFile::Open(fname); // check if file in local directory exists
   } else {
      TFile::SetCacheFileDir(".");
      input = TFile::Open("TrainingData.root", "CACHEREAD");
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }

   // Register the training and test trees
   TTree *signalTree = (TTree *)input->Get("SIGNAL");
   TTree *background = (TTree *)input->Get("BACKGROUND");

   // Create a DataLoader object. This is the class that is used to define the input variables
   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   dataloader->AddVariable("mjj", "m_{jj}", "GeV", 'D');
   dataloader->AddVariable("JetGap", "#Delta y_{jj}", "", 'D');
   dataloader->AddVariable("PtBalance", "pT_{balance}", "", 'D');
   dataloader->AddVariable("centrality", "#xi(Z)", "", 'D');

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   dataloader->AddSpectator("eventNumber",'I');

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree(signalTree,1.0);
   dataloader->AddBackgroundTree(background,1.0);

   // Set individual event weights (the variables must exist in the original TTree)
   // More complicated expressions are possible.
   // -  for instance, for background: dataloader->SetBackgroundWeightExpression("weight1*weight2");
   dataloader->SetSignalWeightExpression("mcWeight");
   dataloader->SetBackgroundWeightExpression("mcWeight");
   
   // Specify the number of events used for training, the rest will be used for evaluation
   TCut preCut = "massTauTau>116 && mcWeight < 1 && mcWeight > -1 && centrality < 1"; // SplitSeed=70 // nTrain_Signal=11750:nTrain_Background=25690
   TString preStr = TString::Format("nTrain_Signal=%d:nTrain_Background=%d:SplitMode=Random:SplitSeed=70:NormMode=EqualNumEvents:V",signalTree->GetEntries(preCut),background->GetEntries(preCut));

   dataloader->PrepareTrainingAndTestTree(preCut,preCut, preStr);
   // Define the number of folds and the splitting expression.
   int NumFolds = 10;
   TString splitExpr = ":SplitType=Deterministic:SplitExpr=int([eventNumber])\%int(10)";

   // Open output file
   TFile *outputFile = TFile::Open("TrainResults0.4.6.root", "RECREATE");
   (TMVA::gConfig().GetVariablePlotting()).fNbinsMVAoutput=10;

   // Initialise the cross-validation object to train and evaluate the MVAs
   TString opt = "!V:!Silent:!Correlations:AnalysisType=Classification:NumFolds=10:FoldFileOutput=True"+splitExpr;
   TMVA::CrossValidation cv {"10Folds", dataloader, outputFile, opt};

   // Book your MVA 
   cv.BookMethod(TMVA::Types::kBDT, "BDT-0.4.6", "!H:!V:NTrees=125:MinNodeSize=3%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost=True:BaggedSampleFraction=0.5:UseNvars=2:nCuts=80:MaxDepth=2:SeparationType=SDivSqrtSPlusB");
   
   // Train and evaluate the MVA
   cv.Evaluate();


   // Print the results for every MVA method.
   size_t iMethod = 0;
   for (auto && result : cv.GetResults()) {
      std::cout << "Summary for method "
               << cv.GetMethods()[iMethod++].GetValue<TString>("MethodTitle")
               << std::endl;
      std::vector<double> ROCValues{};

      // Print results for every fold.
      for (UInt_t iFold = 0; iFold<cv.GetNumFolds(); ++iFold) {
         ROCValues.push_back(result.GetROCValues()[iFold]);
         std::cout << "\tFold " << iFold << ": "
         << "ROC int: " << result.GetROCValues()[iFold]
         << ", "
         << "SignalEff @ BkgEff=0.1: " << result.GetEff10Values()[iFold]
         << ", "
         << "Separation: " << result.GetSepValues()[iFold]
         << std::endl;
      } 
      double mean{};
      double stddev{};
      std::cout << "SUMMARY: " << std::endl;
      AverageAndStdDeviation(ROCValues, mean, stddev);
      std::cout << "ROC Integral: " << mean << " +- " << stddev << std::endl;
      AverageAndStdDeviation(result.GetSepValues(), mean, stddev);
      std::cout << "Separation: " << mean << " +- " << stddev << std::endl << std::endl;

   }

   // Clean up
   outputFile->Close();
}
