"""
65;6800;1cZG->NuNuG analysis using Run3 NanoAOD dataset
"""
import os
import logging
logger = logging.getLogger("Bamboo tutorial")

from functools import partial
from itertools import chain
from bamboo.plots import Skim

from bamboo import treefunctions as op
from bamboo import treedecorators as td

from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo.plots import CutFlowReport

from bamboo.analysismodules import NanoAODModule, SkimmerModule, HistogramsModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection, makePileupWeight
from bamboo.plots import VariableBinning as VarBin


def pogEraFormat(era):
    if any( x in era for x in ['2022', '2023']): return era[:4] +'_Summer'+era.replace('20','')
    else: return era.replace("UL", "") + "_UL"


# https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/
def localizePOGSF(era, POG, fileName):
    subdir = pogEraFormat(era)
    return os.path.join("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration", "POG", POG, subdir, fileName)


def makePUWeight(tree, era, selection):
    if era == '2023': goldenJSON= "Collisions2023_366403_369802_eraBC_GoldenJson" 
    elif era =='2023BPix': goldenJSON= "Collisions2023_369803_370790_eraD_GoldenJson"
    puTuple = (localizePOGSF(era, "LUM", "puWeights.json.gz"), goldenJSON)
    return makePileupWeight(puTuple, tree.Pileup_nTrueInt, systName="pileup", sel=selection)

def getNanoAODDescription(era, isMC, doRocCor=True):
    metName = "PuppiMET"
    nanoJetMETCalc_both = td.CalcCollectionsGroups(
            Jet=("pt", "mass"), changes={metName: (f"{metName}T1", f"{metName}T1Smear")},**{metName: ("pt", "phi")})
    nanoJetMETCalc_data = td.CalcCollectionsGroups(
            Jet=("pt", "mass"), changes={metName: (f"{metName}T1",)}, **{metName: ("pt", "phi")})
    systVars=[td.nanoFatJetCalc ] +[nanoJetMETCalc_both if isMC else nanoJetMETCalc_data]
    if doRocCor:systVars.append(td.nanoRochesterCalc)
    return td.NanoAODDescription.get("v12", year=era[:4], isMC=isMC, systVariations=systVars)

class baseBambooTutorialModule(NanoAODModule, HistogramsModule):
    def __init__(self, args):
        super(baseBambooTutorialModule, self).__init__(args)

        self.plotDefaults = {
            "y-axis": "Events",
            "log-y": "both",
            "y-axis-show-zero": True,
            "log-y-axis-range": [10e-4, 10e2],
            "save-extensions": ["pdf", "png"],
            "show-ratio": True,
            "sort-by-yields": False,
            "legend-columns": 1
            }

    def addArgs(self, parser):
        super(baseBambooTutorialModule, self).addArgs(parser)

        parser.add_argument("-s", "--systematic", action="store_true", default=False, help="Produce systematic variations")
        parser.add_argument("--roc-corr", action="store_true", default=False, help="Enable muon Rochester correction")
        parser.add_argument("--samples", nargs='*', required=True, help="Sample template YML file")


    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        era = sampleCfg["era"]
        tree, noSel, be, lumiArgs = super().prepareTree(tree, sample=sample,
                                                        sampleCfg=sampleCfg,
                                                        description=getNanoAODDescription(era, self.isMC(sample), doRocCor=self.args.roc_corr),
                                                        backend="dataframe"
                                                        )
        # Triggers selection
        triggersPerPrimaryDataset = {
                    "EGamma"    : [ tree.HLT.Photon200,
                                    ]
                    }

        # Cross-section weight
        if self.isMC(sample):
            xsecWeight = sampleCfg["cross-section"]
    
        lumiWeight = lumiArgs[0]
        
        if self.isMC(sample):
            #mcWeight = tree.genWeight * makePUWeight(tree, era, noSel) * lumiWeight #No need to explicitly multiply lumi
            #noSel = noSel.refine("mcWeight", weight=mcWeight, autoSyst=self.args.systematic)
            noSel = noSel.refine("mcWeight", weight=tree.genWeight, autoSyst=self.args.systematic)
            noSel = noSel.refine("puWeight", weight=makePUWeight(tree, era, noSel))
            noSel = noSel.refine("withtriggers", cut=(op.OR(*chain.from_iterable(triggersPerPrimaryDataset.values()))))
        else:
            noSel = noSel.refine("withtriggers", cut=(makeMultiPrimaryDatasetTriggerSelection(sample, triggersPerPrimaryDataset)))

        # Calculate (corrected) photon 4-momenta before accessing them
        if self.args.roc_corr:
            forceDefine(tree._Photon.calcProd, noSel)
        return tree, noSel, be, lumiArgs


    def defineObjects(self, tree, noSel, sample=None, sampleCfg=None, **kwargs):
        era = sampleCfg["era"]
        self.sorted_photon = op.sort(tree.Photon, lambda ph: -ph.pt)

        hlt_photons = op.select(tree.TrigObj, lambda trig: trig.id == 22)
        def hlt_match(ph):
            return op.rng_any(hlt_photons, lambda trig: op.deltaR(trig.p4, ph.p4) < 0.3)
        
        #Photon selection for EB
        self.photons_EB = op.select(self.sorted_photon, lambda ph : op.AND(
            ph.pt > 225.0,
            op.abs(ph.eta)<1.4442,
            ph.sieie > 0.001,
            ph.sipip > 0.001,
            ph.pixelSeed < 1.0,
            ph.etaWidth > 0.01,
            ph.cutBased >= 2.0,
            #ph.trkSumPtHollowConeDR03 < 4.98616,
            #ph.ecalPFClusterIso < 6.88295,
            #ph.hoe < 0.023136,
            #ph.hcalPFClusterIso < 14.9892,
            hlt_match(ph)))
        
        #photon selection for EE
        self.photons_EE = op.select(self.sorted_photon, lambda ph : op.AND(
            ph.pt > 225.0,
            ph.pixelSeed < 1.0,
            op.AND(op.abs(ph.eta) > 1.566, op.abs(ph.eta) < 2.5),
            ph.cutBased >= 2.0,
            ph.haloTaggerMVAVal > 0.996,
            #ph.trkSumPtHollowConeDR03 < 3.75,
            #ph.ecalPFClusterIso < 5.13,
            #ph.hoe < 0.0134,
            #ph.hcalPFClusterIso < 14.9892,            
            hlt_match(ph)))

        #MET cut
        self.selected_MET = tree.PuppiMET
        self.met_cut = self.selected_MET.pt > 200

        #Jet veto
        self.selected_jets = op.select(tree.Jet, lambda jet: op.AND(
            jet.pt > 30.0,             # Jet pT > 30 GeV
            op.abs(jet.eta) < 5.0,     # |η| < 5
            jet.jetId >= 2,            # Tight ID
            #jet.puIdDisc > 0,
            op.NOT(op.rng_any(self.sorted_photon, lambda ph: op.deltaR(jet.p4, ph.p4) < 0.5)),  # ΔR(γ, jet) > 0.5
            op.abs(jet.phi - self.selected_MET.phi) < 0.5   #Δφ(jet, MET) < 0.5
        ))
        self.jet_veto = op.rng_len(self.selected_jets) == 0  # No jets passing selection
        
        #selection photon with met and jet veto
        self.photons_EB_met = op.select(self.photons_EB, lambda ph: op.AND(
            self.met_cut,                                # MET > 200 GeV
            ph.pt / self.selected_MET.pt < 1.4,         # Photon pT / MET < 1.4
            op.abs(self.selected_MET.phi - ph.phi) > 2,  # Δφ(MET, γ) > 2
            self.jet_veto
        ))
        self.photons_EE_met = op.select(self.photons_EE, lambda ph: op.AND(
            self.met_cut,
            ph.pt / self.selected_MET.pt < 1.4,
            op.abs(self.selected_MET.phi - ph.phi) > 2,
            self.jet_veto
        ))

        '''# ABCD method
        # Choosing region C(Signal region)
        #Phase space cuts + photon ID cuts
        
        # Choosing region C(Signal region)
        self.photons_EB_C = op.select(self.photons_EB_met, lambda ph: op.AND(
            ph.trkSumPtHollowConeDR03 < 4.98616,
            ph.mvaID > 0.9 #BDT cut only for region C
        ))


        # Region A  full ID as region C except BDT score and MET and lepton veto
        self.met_cut_A = self.selected_MET.pt <40
        # muons                                                                                                                                          
        self.muons = op.select(tree.Muon, lambda mu: op.AND(        
            mu.pt > 30.0,                                                                                      
            op.abs(mu.eta) < 2.4,                                                                               
            op.rng_any(self.sorted_photon, lambda ph: op.deltaR(mu.p4, ph.p4) > 0.5),                                                                  
            mu.tightId,                                                                                             
            mu.pfRelIso03_all < 0.4                                                                                                         
        ))                                                                                
        #electron
        self.electrons = op.select(tree.Electron, lambda el: op.AND(
            el.pt > 30.0,
            op.OR(
                op.abs(el.eta) < 1.4442,
                op.AND(op.abs(el.eta) > 1.566, op.abs(el.eta) < 2.5)
            ),
            op.rng_any(self.sorted_photon, lambda ph: op.deltaR(el.p4, ph.p4) > 0.5),
            el.cutBased == 4
        ))

        self.lepton_veto = op.AND(
            op.rng_len(self.muons) == 0,
            op.rng_len(self.electrons) == 0
        )
        self.photons_EB_A = op.select(self.photons_EB_met, lambda ph: op.AND(
            ph.trkSumPtHollowConeDR03 < 4.98616,
            self.lepton_veto
        ))

        self.photons_EB_AS = op.select(self.photons_EB_met, lambda ph : op.AND(
            #ph.trkSumPtHollowConeDR03 > 5,  
            ph.trkSumPtHollowConeDR03 < 10,
            self.lepton_veto
        ))'''              

        return

     
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None, **kwargs):
        super().postProcess(taskList, config, workdir, resultsdir, **kwargs)


class ZGto2NuGPlotter(baseBambooTutorialModule):
    """ Class to create control plots, cutflow reports and skims"""

    def __init__(self, args):
        super(ZGto2NuGPlotter, self).__init__(args)


    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        self.defineObjects(t, noSel, sample, sampleCfg)
        plots = []
        # cut flow report
        self.cfr = CutFlowReport("yields", recursive=True, printInLog=True)
        self.cfr.add(noSel, 'no sel')
        plots.append(self.cfr)
        
        if self.isMC(sample):
            weightsel = noSel.weight  # Use weight for MC
        else:
            weightsel = None  # No weight for data

        #Apply cuts for selection  
        has_photon_EB = noSel.refine('hasphotons_inEB', cut=(op.rng_len(self.photons_EB_met) == 1))
        has_photon_EE = noSel.refine('hasphotons_inEE', cut=(op.rng_len(self.photons_EE_met) == 1))
        #noJetSel = noSel.refine('noJets', cut=(op.rng_len(self.jets_noclean) == 0))
        asym_bins_pT = VarBin([225, 275, 350, 450, 600, 800, 1500])
        asym_bins_met_pT = VarBin([225, 275, 350, 450, 600, 1000])
        
        plots.append(Plot.make1D("EB_Pt", op.map(self.photons_EB_met, lambda ph: ph.pt), has_photon_EB,asym_bins_pT, weight=weightsel, title="photon p_{T} (GeV)", autoSyst=True))
        plots.append(Plot.make1D("EB_eta", op.map(self.photons_EB_met, lambda ph: ph.eta), has_photon_EB,EqBin(15, -2., 2.), weight=weightsel, title="photon eta", autoSyst=True))
        plots.append(Plot.make1D("EB_phi", op.map(self.photons_EB_met, lambda ph: ph.phi), has_photon_EB,EqBin(16, -3., 3.), weight=weightsel, title="photon phi", autoSyst=True))
        plots.append(Plot.make1D("EB_met_pT", (self.selected_MET.pt), has_photon_EB,asym_bins_met_pT, weight=weightsel, title="missing p_{T}(GeV)", autoSyst=True))
        plots.append(Plot.make1D("EB_delta_phi",op.abs(self.selected_MET.phi - self.photons_EB_met[0].phi), has_photon_EB,EqBin(10, 2., 4.), weight=weightsel, title="Δφ(MET,γ)", autoSyst=True))
        
        plots.append(Plot.make1D("EE_Pt", op.map(self.photons_EE_met, lambda ph: ph.pt), has_photon_EE,asym_bins_pT, weight=weightsel, title="photon p_{T} (GeV)", autoSyst=True))
        plots.append(Plot.make1D("EE_eta", op.map(self.photons_EE_met, lambda ph: ph.eta), has_photon_EE,EqBin(15, -3.2, 3.2), weight=weightsel, title="photon eta", autoSyst=True))
        plots.append(Plot.make1D("EE_phi", op.map(self.photons_EE_met, lambda ph: ph.phi), has_photon_EE,EqBin(16, -3., 3.), weight=weightsel, title="photon phi", autoSyst=True))
        plots.append(Plot.make1D("EE_met_pT", (self.selected_MET.pt), has_photon_EE,asym_bins_met_pT, weight=weightsel, title="missing p_{T}(GeV)", autoSyst=True))
        plots.append(Plot.make1D("EE_delta_phi",op.abs(self.selected_MET.phi - self.photons_EE_met[0].phi), has_photon_EE,EqBin(10, 2., 4.), weight=weightsel, title="Δφ(MET,γ)", autoSyst=True))

        if not self.isMC(sample):
            has_photon_EB_pt_400 = noSel.refine("has_photon_EB_400",op.select(self.photons_EB_met,lambda ph: ph.pt >400))
            has_photon_EE_pt_400 = noSel.refine("has_photon_EE_400",op.select(self.photons_EE_met,lambda ph: ph.pt >400))
            plots.append(Plot.make2D("EB_eta_phi",(op.map(self.photons_EB_met, lambda ph: ph.eta),op.map(self.photons_EB_met, lambda ph: ph.phi)),has_photon_EB_pt_400,(EqBin(15, -2., 2.), EqBin(16, -3.2, 3.2)),weight=weightsel,title=("η vs φ")))
            plots.append(Plot.make2D("EE_eta_phi",(op.map(self.photons_EE_met, lambda ph: ph.eta),op.map(self.photons_EE_met, lambda ph: ph.phi)),has_photon_EE_pt_400,(EqBin(15, -3.2, 3.2), EqBin(16, -3.2, 3.2)),weight=weightsel,title=("η vs φ")))
                     
        '''#For ABCD method
        has_photon_EB_C = noSel.refine('hasphotons_inEB_C', cut=(op.rng_len(self.photons_EB_C) == 1))
        has_photon_EB_A = noSel.refine('hasphotons_inEB_A', cut=(op.rng_len(self.photons_EB_A) == 1))
        has_photon_EB_AS = noSel.refine('hasphotons_inEB_AS', cut=(op.rng_len(self.photons_EB_AS) == 1))

        plots.append(Plot.make1D("BDT_EB_C",op.map(self.photons_EB_C, lambda ph: ph.mvaID),has_photon_EB_C,EqBin(50, -1., 1.),weight=weightsel,title="BDT"))
        plots.append(Plot.make1D("BDT_EB_A",op.map(self.photons_EB_A, lambda ph: ph.mvaID),has_photon_EB_A,EqBin(50, -1., 1.),weight=weightsel,title="BDT"))
        plots.append(Plot.make1D("BDT_EB_AS",op.map(self.photons_EB_AS, lambda ph: ph.mvaID),has_photon_EB_AS,EqBin(50, -1., 1.),weight=weightsel,title="BDT sideband"))
        
        plots.append(Skim("Bdt_A", {
            "run": None,  # copy from input                                                                                                       
            "luminosityBlock": None,
            "event": None,
            "mvaID": op.map(self.photons_EB_A, lambda ph: ph.mvaID)}, has_photon_EB_A, keepOriginal=[Skim.KeepRegex("PV_.*"), "nOtherPV", Skim.KeepRegex("OtherPV_.*")]))

        plots.append(Skim("Bdt_data_sideband", {
            "run": None,  # copy from input                   
            "luminosityBlock": None,
            "event": None,
            "mvaID": op.map(self.photons_EB_AS, lambda ph: ph.mvaID)}, has_photon_EB_AS, keepOriginal=[Skim.KeepRegex("PV_.*"), "nOtherPV", Skim.KeepRegex("OtherPV_.*")]))'''


        return plots
    
class ZGto2NuGSkimmer(SkimmerModule,baseBambooTutorialModule):
    """ Class to create control plots, cutflow reports, and skims """

    def __init__(self, args):
        super(ZGto2NuGSkimmer, self).__init__(args)

    def defineSkims(self, t, noSel, sample=None, sampleCfg=None):
        """ Define skims for events passing selection """

        self.defineObjects(t, noSel, sample, sampleCfg)  # Define objects (photons, jets, MET, etc.)

        #BDT spectrum for A in data
        has_photon_EB_A = noSel.refine('hasphotons_inEB_A', cut=(op.rng_len(self.photons_EB_A) > 0))
        has_photon_EB_AS = noSel.refine('hasphotons_inEB_AS', cut=(op.rng_len(self.photons_EB_AS) > 0))
        
        '''sampleGroup = sampleCfg.get("group","")

        if sampleGroup == "data":
            for region, photons, sel in [
                    ("A", self.photons_EB_A, has_photon_EB_A),
                    ("AS", self.photons_EB_AS, has_photon_EB_AS),
            ]:
                bdt = op.map(photons, lambda ph: ph.mvaID)
                skimVars = {
                    "event": t.event,
                    "run": t.run,
                    "lumi": t.luminosityBlock,
                    "photon_BDT": bdt,
                }
                self.addSkim(f"photon_EB_{region}_skim", sel, skimVars)
                
        elif sampleGroup == "GJ":'''
        bdt_mc = op.map(self.photons_EB_A, lambda ph: ph.mvaID)
        skimVarsMC = {
            "event": t.event,
            "run": t.run,
            "lumi": t.luminosityBlock,
            "photon_BDT": bdt_mc,
        }
        self.addSkim(
            name = "photon_EB_A_MC_skim",
            selection = has_photon_EB_A,
            branches = skimVarsMC,
            treeName = "photonEvents"  # Optional: changes TTree name
        )




                
