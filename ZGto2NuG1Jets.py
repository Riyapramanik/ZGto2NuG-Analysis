"""
ZG->NuNuG analysis using Run3 NanoAOD dataset
"""
import os
import logging
logger = logging.getLogger("Bamboo tutorial")

from functools import partial
from itertools import chain

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
            "legend-columns": 2
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

        
        if self.isMC(sample):
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
        #Photon selection for EB
        self.photons_EB = op.select(self.sorted_photon, lambda ph : op.AND(
            ph.pt > 225.0,
            op.abs(ph.eta)<1.4442,
            ph.sieie > 0.001,
            ph.sipip > 0.001,
            ph.pixelSeed < 1.0,
            ph.etaWidth > 0.01))
        #photon selection for EE
        self.photons_EE = op.select(self.sorted_photon, lambda ph : op.AND(
            ph.pt > 225.0,
            op.AND(op.abs(ph.eta) > 1.566, op.abs(ph.eta) < 2.5)))
        
        return

     
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None, **kwargs):
        super().postProcess(taskList, config, workdir, resultsdir, **kwargs)


class Hto4lPlotter(baseBambooTutorialModule):
    """ Class to create control plots, cutflow reports and skims"""

    def __init__(self, args):
        super(Hto4lPlotter, self).__init__(args)


    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        self.defineObjects(t, noSel, sample, sampleCfg)
        
        plots = []
        # cut flow report                                                                                                                 
        self.cfr = CutFlowReport("yields", recursive=True, printInLog=True)
        self.cfr.add(noSel, 'no sel')
        plots.append(self.cfr)
        photon_pt_EB = self.photons_EB[0].pt
        plot_photon_pT_EB = Plot.make1D(f"photon_pT_EB", photon_pt_EB, noSel, EqBin(1000, 200., 1000.), title="pT_{GeV} (GeV)",
                        plotopts={"show-overflow": False, "log-y": True}
        )
        photon_pt_EE = self.photons_EE[0].pt
        plot_photon_pT_EE = Plot.make1D(f"photon_pT_EE", photon_pt_EE, noSel, EqBin(1000, 200., 1000.), title="pT_{GeV} (GeV)",
			plotopts={"show-overflow": False, "log-y": True}
        )
        plots.append(plot_photon_pT_EB)
        plots.append(plot_photon_pT_EE)
        return plots
