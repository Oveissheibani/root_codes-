#include <iostream>
#include <vector>
#include <string>
#include <filesystem> // For scanning directories
#include <numeric>    // For std::accumulate and std::inner_product
#include <cmath>      // For std::sqrt
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TROOT.h>
#include <TParameter.h>
#include <TTree.h>
#include <TClass.h>
#include <iomanip>    // For std::setw and std::setfill

namespace fs = std::filesystem; // Alias for easier usage

// Function to display a progress bar
void PrintProgressBar(int current, int total) {
    int barWidth = 70;
    float progress = float(current) / float(total);
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

// Function to recursively merge directories
void MergeDirectories(TDirectory* sourceDir, const std::vector<TFile*>& inputFiles, TDirectory* outputDir);

// Main function to merge ROOT files from all available directories
void MergeSingleGenFiles() {
    std::vector<TFile*> inputFiles;
    std::string outputFileName = "PairGenMerged.root";

    // Scan the current directory for folders containing PairGen.root
    for (const auto& entry : fs::directory_iterator(".")) {
        if (entry.is_directory()) {
            std::string fileName = entry.path().string() + "/PairGen.root";
            TFile* file = TFile::Open(fileName.c_str());
            if (file && !file->IsZombie()) {
                std::cout << "File " << fileName << " is found and opened successfully." << std::endl;
                inputFiles.push_back(file);
            } else {
                std::cerr << "File " << fileName << " not found or is corrupted!" << std::endl;
            }
        }
    }

    // Number of files successfully opened
    int nFiles = inputFiles.size();

    // Proceed with merging if at least one file is found
    if (nFiles < 1) {
        std::cerr << "No files found for merging." << std::endl;
        return;
    }

    // Create the output file
    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Failed to create the output file " << outputFileName << std::endl;
        for (auto* file : inputFiles) file->Close();
        return;
    }

    std::cout << "Merging files into " << outputFileName << "..." << std::endl;

    // Start merging objects from the first file
    TDirectory* firstDir = inputFiles[0];
    firstDir->cd();
    TIter nextKey(firstDir->GetListOfKeys());
    TKey* key;

    int objCount = 0;                                // For tracking object count
    int totalObjs = firstDir->GetListOfKeys()->GetEntries(); // Total objects to process

    while ((key = (TKey*)nextKey())) {
        TObject* obj = key->ReadObj();
        std::string objName = obj->GetName();
        TClass* objClass = obj->IsA();

        // Process the object depending on its type
        if (obj->InheritsFrom(TH1::Class())) {
            // It's a histogram, process accordingly

            // Clone the histogram to the output file
            outputFile->cd();
            TH1* histClone = (TH1*)obj->Clone();
            histClone->Reset(); // Clear the clone for merging

            // Get the number of bins (including underflow and overflow)
            int nBinsX = histClone->GetNbinsX() + 2;
            int nBinsY = histClone->GetNbinsY() + 2; // For 2D histograms
            int nBinsZ = histClone->GetNbinsZ() + 2; // For 3D histograms

            // Loop over bins and compute mean and standard deviation
            for (int i = 0; i < nBinsX; ++i) {
                for (int j = 0; j < nBinsY; ++j) {
                    for (int k = 0; k < nBinsZ; ++k) {
                        std::vector<double> binContents;

                        // Collect bin contents from all files
                        for (auto* file : inputFiles) {
                            TH1* h = (TH1*)file->Get(histClone->GetName());
                            if (h) {
                                double content = h->GetBinContent(i, j, k);
                                binContents.push_back(content);
                            } else {
                                std::cerr << "Histogram " << histClone->GetName()
                                          << " not found in file " << file->GetName() << std::endl;
                            }
                        }

                        if (!binContents.empty()) {
                            // Compute mean and standard deviation
                            double sum = std::accumulate(binContents.begin(), binContents.end(), 0.0);
                            double mean = sum / binContents.size();

                            double sq_sum = std::inner_product(binContents.begin(), binContents.end(), binContents.begin(), 0.0);
                            double stddev = std::sqrt(sq_sum / binContents.size() - mean * mean);

                            // Set bin content and error
                            histClone->SetBinContent(i, j, k, mean);
                            histClone->SetBinError(i, j, k, stddev);
                        }
                    }
                }
            }

            // Write the histogram with updated errors to the output file
            histClone->Write();

        } else if (obj->InheritsFrom(TParameter<Long64_t>::Class()) ||
                   obj->InheritsFrom(TParameter<int>::Class()) ||
                   obj->InheritsFrom(TParameter<double>::Class()) ||
                   obj->InheritsFrom(TParameter<float>::Class())) {
            // It's a TParameter object (of various numeric types)

            // We need to decide how to merge TParameters
            // For numeric TParameters, we can sum them or average them
            // Here, we'll sum them if they represent counts, or average if they represent means

            // For this example, let's sum them
            double totalValue = 0.0;

            for (auto* file : inputFiles) {
                TParameter<double>* param = (TParameter<double>*)file->Get(objName.c_str());
                if (param) {
                    totalValue += param->GetVal();
                } else {
                    std::cerr << "Parameter " << objName
                              << " not found in file " << file->GetName() << std::endl;
                }
            }

            // Create a new TParameter with the total value
            TParameter<double> totalParam(objName.c_str(), totalValue);
            outputFile->cd();
            totalParam.Write();

        } else if (obj->InheritsFrom(TTree::Class())) {
            // It's a TTree
            // Merging TTrees properly requires a different approach.
            // We'll chain them together and write the merged tree.

            TChain chain(objName.c_str());
            for (auto* file : inputFiles) {
                TTree* tree = (TTree*)file->Get(objName.c_str());
                if (tree) {
                    chain.Add((std::string(file->GetName()) + "/" + objName).c_str());
                } else {
                    std::cerr << "Tree " << objName
                              << " not found in file " << file->GetName() << std::endl;
                }
            }

            outputFile->cd();
            TTree* mergedTree = chain.CloneTree(-1, "fast"); // Clone all entries
            mergedTree->Write();

        } else if (obj->InheritsFrom(TDirectory::Class())) {
            // It's a directory, we need to recursively process its contents
            std::cerr << "Processing subdirectory: " << objName << std::endl;

            // Recursively merge directories
            TDirectory* subDir = (TDirectory*)obj;
            outputFile->cd();
            TDirectory* newDir = outputFile->mkdir(subDir->GetName());
            MergeDirectories(subDir, inputFiles, newDir);

        } else {
            // Other types of objects
            // Copy them from the first file
            outputFile->cd();
            obj->Write();
        }

        // Update progress bar
        objCount++;
        PrintProgressBar(objCount, totalObjs);
    }

    std::cout << "\n"; // Newline after progress bar

    // Close all files
    for (auto* file : inputFiles) {
        file->Close();
    }
    outputFile->Close();

    std::cout << "Merging completed successfully." << std::endl;
}

// Function to recursively merge directories
void MergeDirectories(TDirectory* sourceDir, const std::vector<TFile*>& inputFiles, TDirectory* outputDir) {
    outputDir->cd();
    TIter nextKey(sourceDir->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)nextKey())) {
        TObject* obj = key->ReadObj();
        std::string objName = obj->GetName();
        TClass* objClass = obj->IsA();

        if (obj->InheritsFrom(TH1::Class())) {
            // Handle histograms as before
            TH1* histClone = (TH1*)obj->Clone();
            histClone->Reset();

            int nBinsX = histClone->GetNbinsX() + 2;
            int nBinsY = histClone->GetNbinsY() + 2;
            int nBinsZ = histClone->GetNbinsZ() + 2;

            for (int i = 0; i < nBinsX; ++i) {
                for (int j = 0; j < nBinsY; ++j) {
                    for (int k = 0; k < nBinsZ; ++k) {
                        std::vector<double> binContents;

                        for (auto* file : inputFiles) {
                            TDirectory* dir = file->GetDirectory(sourceDir->GetPath());
                            if (dir) {
                                TH1* h = (TH1*)dir->Get(histClone->GetName());
                                if (h) {
                                    double content = h->GetBinContent(i, j, k);
                                    binContents.push_back(content);
                                }
                            }
                        }

                        if (!binContents.empty()) {
                            double sum = std::accumulate(binContents.begin(), binContents.end(), 0.0);
                            double mean = sum / binContents.size();

                            double sq_sum = std::inner_product(binContents.begin(), binContents.end(), binContents.begin(), 0.0);
                            double stddev = std::sqrt(sq_sum / binContents.size() - mean * mean);

                            histClone->SetBinContent(i, j, k, mean);
                            histClone->SetBinError(i, j, k, stddev);
                        }
                    }
                }
            }

            histClone->Write();

        } else if (obj->InheritsFrom(TParameter<Long64_t>::Class()) ||
                   obj->InheritsFrom(TParameter<int>::Class()) ||
                   obj->InheritsFrom(TParameter<double>::Class()) ||
                   obj->InheritsFrom(TParameter<float>::Class())) {
            // Handle TParameter objects
            double totalValue = 0.0;

            for (auto* file : inputFiles) {
                TDirectory* dir = file->GetDirectory(sourceDir->GetPath());
                if (dir) {
                    TParameter<double>* param = (TParameter<double>*)dir->Get(objName.c_str());
                    if (param) {
                        totalValue += param->GetVal();
                    } else {
                        std::cerr << "Parameter " << objName
                                  << " not found in file " << file->GetName() << std::endl;
                    }
                }
            }

            // Create a new TParameter with the total value
            TParameter<double> totalParam(objName.c_str(), totalValue);
            outputDir->cd();
            totalParam.Write();

        } else if (obj->InheritsFrom(TTree::Class())) {
            // Handle TTrees
            TChain chain(objName.c_str());
            for (auto* file : inputFiles) {
                TDirectory* dir = file->GetDirectory(sourceDir->GetPath());
                if (dir) {
                    TTree* tree = (TTree*)dir->Get(objName.c_str());
                    if (tree) {
                        chain.Add((std::string(file->GetName()) + "/" + dir->GetPath() + "/" + objName).c_str());
                    }
                }
            }

            outputDir->cd();
            TTree* mergedTree = chain.CloneTree(-1, "fast");
            mergedTree->Write();

        } else if (obj->InheritsFrom(TDirectory::Class())) {
            // Recursively handle subdirectories
            TDirectory* subDir = (TDirectory*)obj;
            TDirectory* newSubDir = outputDir->mkdir(subDir->GetName());
            MergeDirectories(subDir, inputFiles, newSubDir);

        } else {
            // Copy other objects
            outputDir->cd();
            obj->Write();
        }
    }
}
