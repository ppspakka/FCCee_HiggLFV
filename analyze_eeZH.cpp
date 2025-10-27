#include <stdio.h>
#include <stdlib.h>
#include "Delphes.C"
#include <TMath.h>
#include <TTree.h>
#include <TChain.h>
#include <TError.h>
#include <vector>
#include <TString.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TLorentzVector.h>

using namespace std;

// Simple helper to collect ROOT files from a directory pattern
vector<TString> getFileList(const TString &pattern)
{
    vector<TString> files;
    TString dirname = gSystem->DirName(pattern);
    TString basename = gSystem->BaseName(pattern);

    TSystemDirectory dir(dirname, dirname);
    TList *list = dir.GetListOfFiles();
    if (!list) return files;

    TIter next(list);
    TSystemFile *file;
    while ((file = (TSystemFile *)next()))
    {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.Contains(basename.ReplaceAll("*", "")))
            files.push_back(dirname + "/" + fname);
    }
    delete list;
    return files;
}

double AbsDeltaPhi(double phi1, double phi2)
{
    // use TLorentzVector::DeltaPhi for a robust phi wrapping
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiM(1.0, 0.0, phi1, 0.0);
    v2.SetPtEtaPhiM(1.0, 0.0, phi2, 0.0);
    return fabs(v1.DeltaPhi(v2));
}

void analyze_eeZH(TString infilename)
{
    gErrorIgnoreLevel = kFatal;

    // --- Input tree ---
    TChain *intree = new TChain("Delphes");
    vector<TString> filelist = getFileList(infilename);
    if (filelist.empty())
    {
        printf("No input files found matching pattern: %s\n", infilename.Data());
        return;
    }

    for (auto &filename : filelist)
    {
        printf("Reading %s\n", filename.Data());
        intree->Add(filename.Data());
    }

    Delphes *indelphes = new Delphes(intree);
    intree->SetBranchStatus("*", 0);
    intree->SetBranchStatus("Electron*", 1);
    intree->SetBranchStatus("Muon*", 1);
    intree->SetBranchStatus("MissingET*", 1); // enable MET branch

    // --- Event counters ---
    Long64_t total_events = intree->GetEntries();
    Long64_t pass_Zll = 0;
    Long64_t pass_Hmue = 0;
    Long64_t pass_missET = 0;

    // --- Experiment part ---
    // Categorize the DeltaPhi(e, MET) into bins
    vector<Long64_t> met_dphi_bins(7, 0); // 7 bins: [0,0.1), [0.1,0.3), [0.3,0.5), [0.5,0.7), [0.7,0.9), [0.9,1.1), [1.1,pi)

    // Categorize the absolute of reconstructed Z mass - true Z mass into bins
    vector<Long64_t> zmass_diff_bins(5, 0); // 5 bins: [0,2), [2,4), [4,6), [6,8), [8,10)

    // --- Loop over events ---
    for (Long64_t ievent = 0; ievent < total_events; ievent++)
    {
        indelphes->GetEntry(ievent);
        if (ievent % 10000 == 0) printf("Processing event %lld / %lld\n", ievent, total_events);
        // ===== Step 1: Find Z → l+l− (SFOS) =====
        const double Z_MASS = 91.1876;
        const double ELE_MASS = 0.000511;     // GeV
        const double MU_MASS  = 0.10565837;   // GeV

        double bestZmass = -1.0;
        double bestZdiff = 1e9;
        int idxZ_lep1 = -1;
        int idxZ_lep2 = -1;
        int flavorZ = -1; // 0=e, 1=mu

        // electrons: check all SFOS pairs and keep best (closest to Z)
        for (int i = 0; i < indelphes->Electron_size; i++)
        {
            for (int j = i + 1; j < indelphes->Electron_size; j++)
            {
                if (indelphes->Electron_Charge[i] * indelphes->Electron_Charge[j] >= 0) continue; // require OS
                TLorentzVector l1, l2;
                l1.SetPtEtaPhiM(indelphes->Electron_PT[i], indelphes->Electron_Eta[i], indelphes->Electron_Phi[i], ELE_MASS);
                l2.SetPtEtaPhiM(indelphes->Electron_PT[j], indelphes->Electron_Eta[j], indelphes->Electron_Phi[j], ELE_MASS);
                double mass = (l1 + l2).M();
                double diff = TMath::Abs(mass - Z_MASS);
                if (diff < 10.0 && diff < bestZdiff)
                {
                    bestZmass = mass;
                    bestZdiff = diff;
                    idxZ_lep1 = i;
                    idxZ_lep2 = j;
                    flavorZ = 0;
                }
            }
        }

        // muons: also check all SFOS pairs and update if better than current best
        for (int i = 0; i < indelphes->Muon_size; i++)
        {
            for (int j = i + 1; j < indelphes->Muon_size; j++)
            {
                if (indelphes->Muon_Charge[i] * indelphes->Muon_Charge[j] >= 0) continue;
                TLorentzVector l1, l2;
                l1.SetPtEtaPhiM(indelphes->Muon_PT[i], indelphes->Muon_Eta[i], indelphes->Muon_Phi[i], MU_MASS);
                l2.SetPtEtaPhiM(indelphes->Muon_PT[j], indelphes->Muon_Eta[j], indelphes->Muon_Phi[j], MU_MASS);
                double mass = (l1 + l2).M();
                double diff = TMath::Abs(mass - Z_MASS);
                if (diff < 10.0 && diff < bestZdiff)
                {
                    bestZmass = mass;
                    bestZdiff = diff;
                    idxZ_lep1 = i;
                    idxZ_lep2 = j;
                    flavorZ = 1;
                }
            }
        }

        if (flavorZ == -1) continue; // no Z found
        pass_Zll++;

        // Fill the Z mass difference bin
        if (bestZdiff < 2.0) zmass_diff_bins[0]++;
        else if (bestZdiff < 4.0) zmass_diff_bins[1]++;
        else if (bestZdiff < 6.0) zmass_diff_bins[2]++;
        else if (bestZdiff < 8.0) zmass_diff_bins[3]++;
        else if (bestZdiff < 10.0) zmass_diff_bins[4]++;

        // ===== Step 2: Find one muon + one electron (not from Z) =====
        int idx_mu = -1, idx_e = -1;

        // find muon not used in Z
        for (int m = 0; m < indelphes->Muon_size; m++)
        {
            if (flavorZ == 1 && (m == idxZ_lep1 || m == idxZ_lep2)) continue;
            if (indelphes->Muon_PT[m] > 30.0)
            {
                if (idx_mu == -1) idx_mu = m;
                else { idx_mu = -2; break; } // more than one muon found
            }
        }
        if (idx_mu < 0) continue; // -1=no muon, -2=multiple muons

        // find electron not used in Z
        for (int e = 0; e < indelphes->Electron_size; e++)
        {
            if (flavorZ == 0 && (e == idxZ_lep1 || e == idxZ_lep2)) continue;
            if (indelphes->Electron_PT[e] > 20.0)
            {
                if (idx_e == -1) idx_e = e;
                else { idx_e = -2; break; } // more than one electron found
            }
        }
        if (idx_e < 0) continue; // -1=no electron, -2=multiple electrons

        // Opposite-sign requirement
        if (indelphes->Muon_Charge[idx_mu] * indelphes->Electron_Charge[idx_e] >= 0) continue;

        pass_Hmue++;

        // require electron from Higgs to have |DeltaPhi(e, MET)| < 0.7
        if (indelphes->MissingET_size <= 0) continue; // no MET -> skip
        double dphi_emet = AbsDeltaPhi(indelphes->Electron_Phi[idx_e], indelphes->MissingET_Phi[0]);

        // Bin the MET
        if (dphi_emet < 0.1) met_dphi_bins[0]++;
        else if (dphi_emet < 0.3) met_dphi_bins[1]++;
        else if (dphi_emet < 0.5) met_dphi_bins[2]++;
        else if (dphi_emet < 0.7) met_dphi_bins[3]++;
        else if (dphi_emet < 0.9) met_dphi_bins[4]++;
        else if (dphi_emet < 1.1) met_dphi_bins[5]++;
        else met_dphi_bins[6]++;

        if (dphi_emet >= 0.7) continue;
        pass_missET++;
    }

    // --- Print summary ---
    printf("\n===== Summary =====\n");
    printf("Total events:      %lld\n", total_events);
    printf("Pass Z→ll:         %lld\n", pass_Zll);
    printf("Pass H→μτ(eνν):    %lld\n", pass_Hmue);
    printf("Pass |Δφ(e,MET)|<0.7: %lld\n", pass_missET);
    printf("===================\n");

    // --- Print MET dphi binning ---
    printf("\n===== MET Δφ(e,MET) Binning =====\n");
    const char* bin_ranges[] = {
        "[0,0.1)", "[0.1,0.3)", "[0.3,0.5)", "[0.5,0.7)", "[0.7,0.9)", "[0.9,1.1)", "[1.1,π)"
    };
    for (size_t i = 0; i < met_dphi_bins.size(); i++)
    {
        printf("Bin %s : %lld\n", bin_ranges[i], met_dphi_bins[i]);
    }
    printf("===================================\n");

    // --- Print Z mass difference binning ---
    printf("\n===== |M_reco(Z) - M_true(Z)| Binning =====\n");
    const char* zmass_bin_ranges[] = {
        "[0,2)", "[2,4)", "[4,6)", "[6,8)", "[8,10)"
    };
    for (size_t i = 0; i < zmass_diff_bins.size(); i++)
    {
        printf("Bin %s : %lld\n", zmass_bin_ranges[i], zmass_diff_bins[i]);
    }
    printf("===================================\n");
}
