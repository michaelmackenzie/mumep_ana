#ifndef __mumep_ana_ana_MVATools__
#define __mumep_ana_ana_MVATools__
// Class to initialize a TMVA factory and reader consistently

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TTree.h"
#endif

// local includes
#include "mumep_ana/stn/ana/Tree_t.hh"

namespace mumep_ana {

  class MVATools {
  public:
    /**
       Version information:
       v0: default configuration, only necessary spectators

       TrkQual:
       v1: nactive, activehitsfraction, nullhitsfraction, activematsitesfraction, fitcons, momerr,
     t0err
     **/
    MVATools() {}

    // information for a variable
    struct Var_t {
      TString var_;  // name
      TString desc_; // description
      TString unit_; // units, if any
      float* val_;   // address
      bool use_;     // use or just spectator
      char type_;
      Var_t(TString var, TString desc, TString unit, float* val, bool use = false) : var_(var), desc_(desc), unit_(unit), val_(val), use_(use), type_('F') {}
    };

    // get list of variables for training/evaluating MVAs
    static std::vector<Var_t> GetVariables(TString model, Tree_t& tree, int version = MVATools::Default_) {
      std::vector<Var_t> variables;
      std::vector<TString> train_var;

      if(model == "TrkQual") { // Track quality models
        if(version == 0 || version == 1)
          train_var = {"nactive", "activehitsfraction", "nullhitsfraction", "activematsitesfraction", "fitcons", "momerr", "t0err"};
      }

      if(model == "PID") { // PID models
        if(version == 0 || version == 1)
          train_var = {"eclusteroverptrack", "dt"};
      }

      if(model == "TrkPID") { // Tracker-based PID models
        if(version == 0 || version == 1)
          train_var = {"trktzsloperatio", "trkfitcon", "trkactiveratio", "trknullratio"};
      }

      if(model == "CosmicID") { // Cosmic ID models
        if(version == 0 || version == 1)
          train_var = {"trktzsloperatio", "trkd0", "trkcostheta", "trkrmax"};
      }

      if(train_var.size() == 0) {
        std::cout << "MVATools::" << __func__ << ": No training variables defined for model " << model.Data() << " and version " << version << std::endl;
        throw 20;
      }

      // necessary event information
      variables.push_back(Var_t("eventweight", "eventWeight", "", &tree.fWeight));

      // normal variables
      variables.push_back(Var_t("nactive", "N(active)", "", &tree.fTrkQual_nactive));
      variables.push_back(Var_t("activehitsfraction", "N(active)/N(hits)", "", &tree.fTrkQual_activehitsfraction));
      variables.push_back(Var_t("nullhitsfraction", "N(null)/N(active)", "", &tree.fTrkQual_nullhitsfraction));
      variables.push_back(Var_t("activematsitesfraction", "N(active)/N(mat)", "", &tree.fTrkQual_activematsitesfraction));
      variables.push_back(Var_t("fitcons", "p(#chi^2)", "", &tree.fTrkQual_fitcons));
      variables.push_back(Var_t("momerr", "#sigma(p)", "", &tree.fTrkQual_momerr));
      variables.push_back(Var_t("t0err", "#sigma(t)", "", &tree.fTrkQual_t0err));

      variables.push_back(Var_t("eclusteroverptrack", "E/P", "", &tree.fTrkEP));
      variables.push_back(Var_t("dt", "#deltat", "ns", &tree.fTrkDt));

      variables.push_back(Var_t("trkactiveratio", "N(active)/N(hits)", "", &tree.fTrkActiveRatio));
      variables.push_back(Var_t("trknullratio", "N(null)/N(active)", "", &tree.fTrkNullRatio));
      variables.push_back(Var_t("trktzsloperatio", "(dt/dz) / expected", "", &tree.fTrkTZSlopeRatio));
      variables.push_back(Var_t("trkfitcon", "p(#chi^2)", "", &tree.fTrkFitCon));

      variables.push_back(Var_t("trkd0", "d_{0}", "mm", &tree.fTrkD0));
      variables.push_back(Var_t("trkrmax", "R(max)", "mm", &tree.fTrkRMax));
      variables.push_back(Var_t("trktandip", "tan(dip)", "mm", &tree.fTrkTanDip));
      variables.push_back(Var_t("trkcostheta", "cos(#theta)", "", &tree.fTrkCosTheta));

      // Set flags for variables identified as training variables
      for(TString name : train_var) {
        bool found = false;
        for(Var_t& var : variables) {
          if(var.var_ == name)
            var.use_ = true;
          else
            continue;
          found = true;
          break; // if found, continue to the next training variable name
        }
        if(!found) {
          std::cout << "MVATools::" << __func__ << ": WARNING! Failed to find training variable named: " << name.Data() << std::endl;
        }
      }

      // order the variables as given
      std::vector<Var_t> ordered_vars;
      for(unsigned index = 0; index < train_var.size(); ++index) {
        TString name = train_var[index];
        bool found = false;
        for(Var_t var : variables) {
          if(var.var_ == name) {
            ordered_vars.push_back(var);
            found = true;
            break;
          }
        }
        if(!found) {
          std::cout << "TrkQualInit::" << __func__ << ": ERROR! Failed to find training variable named: " << name.Data() << std::endl;
          throw 20;
        }
      }
      return ordered_vars;
    }

    static int InitializeVariables(TMVA::DataLoader& loader, TString selection, int version = MVATools::Default_) {
      int status = 0;
      Tree_t tree; // not used
      std::vector<Var_t> variables = GetVariables(selection, tree);
      for(unsigned index = 0; index < variables.size(); ++index) {
        Var_t& var = variables[index];
        if(var.use_)
          loader.AddVariable(var.var_.Data(), var.desc_.Data(), var.unit_.Data(), var.type_);
        else if(version != 0)
          loader.AddSpectator(var.var_.Data(), var.desc_.Data(), var.unit_.Data(), var.type_);
        else if(var.var_ == "issignal")
          loader.AddSpectator(var.var_.Data(), var.desc_.Data(), var.unit_.Data(), var.type_);
        else if(var.var_ == "fulleventweightlum")
          loader.AddSpectator(var.var_.Data(), var.desc_.Data(), var.unit_.Data(), var.type_);
        else if(var.var_ == "type")
          loader.AddSpectator(var.var_.Data(), var.desc_.Data(), var.unit_.Data(), var.type_);
        else if(var.var_ == "trainfraction")
          loader.AddSpectator(var.var_.Data(), var.desc_.Data(), var.unit_.Data(), var.type_);
      }
      return status;
    }

    static int InitializeVariables(TMVA::Reader& reader, TString selection, Tree_t& tree, int version = MVATools::Default_) {
      int status = 0;
      std::vector<Var_t> variables = GetVariables(selection, tree);
      for(unsigned index = 0; index < variables.size(); ++index) {
        Var_t& var = variables[index];
        if(var.use_)
          reader.AddVariable(var.var_.Data(), var.val_ /*, var.unit_.Data(), var.type_*/);
        else if(version != 0)
          reader.AddSpectator(var.var_.Data(), var.val_ /*, var.unit_.Data(), var.type_*/);
        else if(var.var_ == "issignal")
          reader.AddSpectator(var.var_.Data(), var.val_);
        else if(var.var_ == "fulleventweightlum")
          reader.AddSpectator(var.var_.Data(), var.val_);
        else if(var.var_ == "type")
          reader.AddSpectator(var.var_.Data(), var.val_);
        else if(var.var_ == "trainfraction")
          reader.AddSpectator(var.var_.Data(), var.val_);
      }
      return status;
    }

    static int SetBranchAddresses(TTree* tree, TString selection, Tree_t& tree_t, int version = MVATools::Default_) {
      int status = 0;
      std::vector<Var_t> variables = GetVariables(selection, tree_t);
      for(unsigned index = 0; index < variables.size(); ++index) {
        Var_t& var = variables[index];
        tree->SetBranchAddress(var.var_.Data(), var.val_);
      }
      return status;
    }

    static void TestVariables(Long64_t entry, TTree* tree, TString selection, Tree_t& vars, int version = MVATools::Default_, bool setAddresses = true) {
      std::vector<Var_t> variables = GetVariables(selection, vars);
      if(setAddresses)
        SetBranchAddresses(tree, selection, vars);
      tree->GetEntry(entry);
      for(unsigned index = 0; index < variables.size(); ++index) {
        Var_t& var = variables[index];
        printf("%20s (%30s) = %12.5f", var.var_.Data(), var.desc_.Data(), *(var.val_));
        if(var.use_)
          printf(" (Variable)\n");
        else if(version != 0)
          printf(" (Spectator)\n");
        if(!std::isfinite(*(var.val_)))
          printf("!!! Variable %s is not finite!\n", var.var_.Data());
      }
    }

    // default version
    const static int Default_ = 0;
  };
} // namespace mumep_ana
#endif
