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
            "log-y-axis-range": [10e-4, 10e8],
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
            mcWeight = tree.genWeight * makePUWeight(tree, era, noSel) * lumiWeight
            noSel = noSel.refine("mcWeight", weight=mcWeight, autoSyst=self.args.systematic)
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
            hlt_match(ph)))
        
        #photon selection for EE
        self.photons_EE = op.select(self.sorted_photon, lambda ph : op.AND(
            ph.pt > 225.0,
            ph.pixelSeed < 1.0,
            op.AND(op.abs(ph.eta) > 1.566, op.abs(ph.eta) < 2.5),
            ph.cutBased >= 2.0,
            hlt_match(ph)))

        #MET cut
        self.selected_MET = tree.MET
        self.met_cut = self.selected_MET.sumEt > 200

        #Jet veto
        self.selected_jets = op.select(tree.Jet, lambda jet: op.AND(
            jet.pt > 30.0,             # Jet pT > 30 GeV
            op.abs(jet.eta) < 5.0,     # |η| < 5
            jet.jetId >= 2,            # Tight ID
            #jet.puId > 0
            op.NOT(op.rng_any(self.sorted_photon, lambda ph: op.deltaR(jet.p4, ph.p4) < 0.5))  # ΔR(γ, jet) > 0.5
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

        
        #Jet selection         
        '''self.muons = op.select(tree.Muon, lambda mu: op.AND(mu.pt > 30., op.abs(mu.eta) < 2.4))
        self.electrons = op.select(tree.Electron, lambda el: op.AND(el.pt > 30, op.abs(el.eta) < 2.5))

        self.jets_noclean = op.select(tree.Jet, lambda j: op.AND(op.abs(j.eta) < 5.0, j.pt > 30.))

        self.jets_lepton = op.sort(
            op.select(self.jets_noclean, lambda j: op.AND(
                op.NOT(op.rng_any(self.muons, lambda l: op.deltaR(l.p4, j.p4) < 0.4)),
                op.NOT(op.rng_any(self.electrons, lambda l: op.deltaR(l.p4, j.p4) < 0.4))
            )), lambda j: -j.pt)

        self.jets_photon = op.sort(
            op.select(self.jets_noclean, lambda j: op.AND(
                op.NOT(op.rng_any(self.sorted_photon, lambda pho: op.deltaR(pho.p4, j.p4) < 0.4))
            )), lambda j: -j.pt)

        self.jets_met = op.select(self.jets_noclean, lambda j: op.NOT(op.abs(op.deltaPhi(j.phi, tree.MET.phi)) < 0.5))
        self.no_jets_passing = op.AND(
            op.rng_len(self.jets_photon) == 0, 
            op.rng_len(self.jets_lepton) == 0, 
            op.rng_len(self.jets_met) == 0
        )'''

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

        plots.append(Plot.make1D("EB_Pt", op.map(self.photons_EB_met, lambda ph: ph.pt), has_photon_EB,
                                 EqBin(200, 200., 1500.), weight=weightsel, title="p_{T} (GeV/c)", autoSyst=False))
        plots.append(Plot.make1D("EB_eta", op.map(self.photons_EB_met, lambda ph: ph.eta), has_photon_EB,
                                 EqBin(100, -2., 2.), weight=weightsel, title="eta", autoSyst=False))
        plots.append(Plot.make1D("EB_phi", op.map(self.photons_EB_met, lambda ph: ph.phi), has_photon_EB,
                                 EqBin(80, -3., 3.), weight=weightsel, title="phi", autoSyst=False))

        plots.append(Plot.make1D("EE_Pt", op.map(self.photons_EE_met, lambda ph: ph.pt), has_photon_EE,
                                 EqBin(100, 200., 1500.), weight=weightsel, title="p_{T} (GeV/c)", autoSyst=False))
        plots.append(Plot.make1D("EE_eta", op.map(self.photons_EE_met, lambda ph: ph.eta), has_photon_EE,
                                 EqBin(100, -3.2, 3.2), weight=weightsel, title="eta", autoSyst=False))
        plots.append(Plot.make1D("EE_phi", op.map(self.photons_EE_met, lambda ph: ph.phi), has_photon_EE,
                                 EqBin(100, -3., 3.), weight=weightsel, title="phi", autoSyst=False))

        #plots.append(Plot.make1D("EB_MET", self.met, has_photon_EB,EqBin(1000, 800., 5000.), title="MET_Et (GeV)"))
        #plots.append(Plot.make1D("EE_MET", self.met, has_photon_EE,EqBin(1000, 800., 5000.), title="MET_Et (GeV)"))
        #plots.append(Plot.make1D("noJets", op.rng_len(self.jets_noclean), noJetSel, EqBin(10, 0., 10.), title="Number of jets"))

        return plots
    
class ZGto2NuGSkimmer(baseBambooTutorialModule, SkimmerModule):
    """ Class to create control plots, cutflow reports, and skims """

    def __init__(self, args):
        super(ZGto2NuGSkimmer, self).__init__(args)

    def defineSkims(self, t, noSel, sample=None, sampleCfg=None):
        """ Define skims for events passing selection """

        self.defineObjects(t, noSel, sample, sampleCfg)  # Define objects (photons, jets, MET, etc.)

        # Apply event selections
        has_photon_EB = noSel.refine('hasphotons_inEB', cut=(op.rng_len(self.photons_EB) > 0))
        has_photon_EE = noSel.refine('hasphotons_inEE', cut=(op.rng_len(self.photons_EE) > 0))
        noJetSel = noSel.refine('noJets', cut=(op.rng_len(self.jets_noclean) == 0))

        # Define skim variables separately for EB and EE photons
        skimVars_EB = {
            "event": t.event,
            "run": t.run,
            "lumi": t.luminosityBlock,
            "nJets": op.rng_len(self.jets_noclean),
            "photon_pt": op.map(self.photons_EB, lambda ph: ph.pt),
            "photon_eta": op.map(self.photons_EB, lambda ph: ph.eta),
            "photon_phi": op.map(self.photons_EB, lambda ph: ph.phi),
            "met": self.met.pt,
        }

        skimVars_EE = {
            "event": t.event,
            "run": t.run,
            "lumi": t.luminosityBlock,
            "nJets": op.rng_len(self.jets_noclean),
            "photon_pt": op.map(self.photons_EE, lambda ph: ph.pt),
            "photon_eta": op.map(self.photons_EE, lambda ph: ph.eta),
            "photon_phi": op.map(self.photons_EE, lambda ph: ph.phi),
            "met": self.met.pt,
        }

        # Add skims
        self.addSkim("photon_EB_skim", has_photon_EB, skimVars_EB)
        self.addSkim("photon_EE_skim", has_photon_EE, skimVars_EE)
        self.addSkim("noJet_skim", noJetSel, skimVars_EB)  # Use EB sk
