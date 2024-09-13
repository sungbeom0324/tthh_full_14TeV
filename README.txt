# Feasibility study of ttHH (2 lepton final state) in HL-LHC. Author : Sungbeom Cho @HYU, Seoul.
# Working space : @/home/stiger97/github/tthh_full_14TeV

This project is for my Master's degree thesis.
This channel has potential to contain rich information about Top/Higgs sector in future accelerators. 

Physics motivation : 
self coupling, fast growing cross section, composite higgs model.

Challenge :
Extremely small cross section & multijet combinatorial problem when reconstructing 2 Higgs.

Strategy : 
Consetutive DNNs.
DNN1 -> bjet from Higgs/Top. (dnn/dnn_bfh.py, dnn_bft.py)
DNN2 -> Signal extraction into 4 groups. (G1:ttHH; G2:ttH, ttZH, ttbbH; G3:ttvv, ttbbv, ttbb, ttbbbb; G4:tttt)

Evaluation :
Comparing the accuracy, Chi^2 vs DNN.
Combine ss2l & os2l in quadrature.

Limitaions :
Systemetics are ignored : Delphes simulation card isn't realistic. (Well.. it is not even built.). No pileup.
Resolved ttbbbb : Very expensive sample (12h/10k events).


# Workflow.

1. Skimming. @skimmer/    -> @skimmed/
    gen.C : Skim only generator level branches.
    ana.C : Skim without event selections.
    ana_ss2l.C : Skim with ss2l event selection criteria. ("Filter") 
    ana_os2l.C : Skim with os2l event selection critetia. ("Filter")

    * The functions are defined @utils/utility.h
    * Prerequisite for reading Delphes branches @classes/, @external/

2. Train DNN1. @dnn/      -> @dnn_result/B/
    dnn_bfh.py : Train to find b-tagged jet pair from Higgs. 
    dnn_bft.py : Train to find b-tagged jet from Top.

3. Train DNN2. @dnn/      -> @dnn_result/SS2l or /OS2l.
    dnn_ss2l.py : Train signal extraction in SLss2l region.
    dnn_os2l.py : Train signal extraction in DLos2l region.

4. Test. @dnn/
    final/dnn_extract.py : Extract signal in inclusive sample. You get the final result.

(5). Cutflow. @cutflow/
    cutflow.py : Simple one.
    cutflow_ss2l.py : Save result in excel format. 
    cutflow_os2l.py : Save result in excel format.

(6). Plot. @palette.py
    * Drawing functions are defined in utils/utils/drawHistoModules.py

# Scripts doing parallel jobs. @scripts/ whose commands at working space is,
    $. scripts -file 
    Exmaple :
        $. scripts/ana.C

