import sys, os
import time
import ast
import argparse
import ROOT as r
r.gROOT.SetBatch(True)
from copy import deepcopy
import CMS_lumi
import myPlotter2D as effplot2D
import myPlotter1D as effplot1D
from allLegends import legends
from subprocess import call
from bcolors import bcolors
from markerColors import markerColors


parser = argparse.ArgumentParser(description='Run the Reco efficiency analysis')
parser.add_argument('-n','--ntuples', action='store_true', default = False)
parser.add_argument('-a','--acceptance', action='store_true', default = False)
parser.add_argument('-nt','--noTimeCut', action='store_false', default = True, dest='timeCut')
parser.add_argument('-p', '--plot', action='store_true',  dest='plot')
parser.add_argument('-q', '--qualities', action='store',  dest='qualities')
parser.add_argument('-d', '--directory', action='store',  dest='directory')
parser.add_argument('-s', '--sample', action='store',  dest='sample')
my_namespace = parser.parse_args()

dir = my_namespace.directory
sam = my_namespace.sample
acc = my_namespace.acceptance
qualities = ast.literal_eval(my_namespace.qualities)
timeCutt0  = my_namespace.timeCut
plot = my_namespace.plot

if acc:
	acctxt = "with acceptance"
else:
	acctxt = "without acceptance"
if timeCutt0:
	timeCuttxt = "with t0 cut"
else:
	timeCuttxt = "without t0 cut"

print("\n----------Computing efficiency for RECO segs " + acctxt + " and " + timeCuttxt + "----------\n")

################################# inputs #######################################

# files = [''+sam+'-std', ''+sam+'-c0' ]

files = [sam]

# qualities = ['All','nothreehits']


##############################################################################################

loc = "DTDPGAnalysis/DTNtuples/test" #DM
print(loc)
print(os.getcwd())

if my_namespace.ntuples == True :
  print ("Starting ntuplizer for every sample in input")
  time.sleep(2)
  r.gInterpreter.ProcessLine(".x loadTPGSimAnalysis_RecoEffs.C") 
  r.gSystem.Load(os.getcwd() + "/DTNtupleBaseAnalyzer_C.so")
  r.gSystem.Load(os.getcwd() + "/DTNtupleTPGSimAnalyzer_RecoEff_C.so")
  from ROOT import DTNtupleTPGSimAnalyzer
else :
  print("Not making ntuples. Make sure they exist in the ntuple directory before running the analysis.")
  time.sleep(2)

##############################################################################################

path = dir+'/'

effPath = './plotsRecoEff/'
outputPath = './ntuplesReco/'


if not os.path.isdir(effPath) : rc = call('mkdir ' + effPath, shell=True)
if not os.path.isdir(outputPath) : rc2 = call('mkdir ' + outputPath, shell=True)

outputPlotDir = effPath + timeCuttxt.replace(' ', '') + "/" 

if not os.path.isdir(outputPlotDir) : rc = call('mkdir ' + outputPlotDir, shell=True)

for quality in qualities:
	if not os.path.isdir(outputPlotDir + quality + "_" + acctxt.replace(' ', '_')) : rc = call('mkdir ' + outputPlotDir + quality + "_" + acctxt.replace(' ', '_'), shell=True)

###### run eff analysis ########

for fil in files :
 if my_namespace.ntuples == True :
  for quality in qualities :
   print ('Obtaining efficiency ntuples for ' + fil + ' with quality type ' + quality )
   time.sleep(2)
   analysis = DTNtupleTPGSimAnalyzer(path + fil + '.root', outputPath + 'results_effis_' +fil + '_' + quality + '.root', quality, acc, timeCutt0)
   analysis.Loop()

############ Plotting Efficiency #################

plotting = { 'OutputFilename' : 'EffPhiEta',
	          'lowlimityaxis': -1.4,
		      'highlimityaxis': 1.4,
			  'lowlimitxaxis': -3.2,
		      'highlimitxaxis': 3.2,
			  'zmin' : 0.0,
			  'zmax' : 1.0,
		      'yaxistitle' : "#eta",
		      'yaxistitleoffset': 0.8,
		      'xaxistitle': "#phi",
		      'legxlow' : 0.6,
		      'legylow': 0.3,
		      'legxhigh': .95,
		      'legyhigh': .44,
		      'markertypedir':{},
		      'markercolordir':{},
                      'PU':'200',
			  'sample': sam,
   		    }

plotEta = { 'OutputFilename' : 'EffEta',
	          'lowlimityaxis': 0.2,
		      'highlimityaxis': 1.0,
			  'lowlimitxaxis': -1.4,
		      'highlimitxaxis': 1.4,
		      'markersize': 1,
		      'yaxistitle' : 'Efficiency (adim.)',
		      'yaxistitleoffset': 1.2,
		      'xaxistitle': "#eta",
		      'legxlow' : 0.6,
		      'legylow': 0.3,
		      'legxhigh': .95,
		      'legyhigh': .44,
		      'markertypedir':{},
		      'markercolordir':{},
                      'PU':'200',
			  'sample': sam,
   		    }

plotPhi = { 'OutputFilename' : 'EffPhi',
	          'lowlimityaxis': 0.2,
		      'highlimityaxis': 1,
			  'lowlimitxaxis': -3.2,
		      'highlimitxaxis': 3.2,
		      'markersize': 1,
		      'yaxistitle' : 'Efficiency (adim.)',
		      'yaxistitleoffset': 1.2,
		      'xaxistitle': "#phi",
		      'legxlow' : 0.6,
		      'legylow': 0.3,
		      'legxhigh': .95,
		      'legyhigh': .44,
		      'markertypedir':{},
		      'markercolordir':{},
                      'PU':'200',
			  'sample': sam,
   		    }

plotPt = { 'OutputFilename' : 'EffPt',
	          'lowlimityaxis': 0.2,
		      'highlimityaxis': 1,
			  'lowlimitxaxis': 20,
		      'highlimitxaxis': 160,
		      'markersize': 1,
		      'yaxistitle' : 'Efficiency (adim.)',
		      'yaxistitleoffset': 1.2,
		      'xaxistitle': "pt (GeV)",
		      'legxlow' : 0.6,
		      'legylow': 0.3,
		      'legxhigh': .95,
		      'legyhigh': .44,
		      'markertypedir':{},
		      'markercolordir':{},
                      'PU':'200',
			  'sample': sam,
   		    }

plotting1D = { 'Eta': plotEta, 'Phi': plotPhi, 'Pt': plotPt }
# plotting1D = {'Pt': plotPt }

plottingZoomed = deepcopy(plotting)
plottingZoomed['zmin'] = 0.9
plottingZoomed['zmax'] = 1.0
plottingZoomed['OutputFilename'] = 'EffPhiEtaZoomed'
plotting2D = { 'EffPhiEta' : plotting, 'EffPhiEtaZoomed' : plottingZoomed}

if sam == 'dy':
	plottingZoomed['energy'] = '14 TeV'
	plotting['energy'] = '14 TeV'
	plotEta['energy'] = '14 TeV'
	plotPhi['energy'] = '14 TeV'
	plotPt['energy'] = '14 TeV'
elif sam == 'pg':
	plottingZoomed['energy'] = ''
	plotting['energy'] = ''
	plotEta['energy'] = ''
	plotPhi['energy'] = ''
	plotPt['energy'] = ''

def Plot2DEffs(files, qualities, plotting, effPath, outputPath):
    chambers = ["MB1", "MB2", "MB3", "MB4"]
    for fil in files :
        for quality in qualities:
            for ch in chambers:
                pathToFile = outputPath + 'results_effis_' + fil + '_' + quality + '.root'
                if not os.path.isfile(pathToFile) :
                    print (bcolors.red + 'ERROR: ' + pathToFile + ' does not exist, maybe running the ntuple production helps' + bcolors.reset)
                    continue
                print ('Obtaining 2D efficiency plots for ' + fil + ' with quality type ' + quality + ' and chamber ' + ch)
                effplot2D.plot2DEff(pathToFile, outputPlotDir + quality + "_" + acctxt.replace(' ', '_'), plotting, "EffPhiEta", ch, "RECO", effPath, fil, quality)
            if not acc:
                for plottingScaffold in plotting2D.values():
                    effplot2D.plot2DEff(pathToFile, outputPlotDir + quality + "_" + acctxt.replace(' ', '_'), plottingScaffold, "EffPhiEta", "COMBINED", "RECO", effPath, fil, quality)
    return

def Plot1DEffs(files, qualities, plotting1D, effPath, outputPath):
	chambers = ["MB1", "MB2", "MB3", "MB4"]
	for fil in files :
		for ch in chambers:
			for quality in qualities:
				pathToFile = outputPath + 'results_effis_' + fil + '_' + quality + '.root'
				if not os.path.isfile(pathToFile) :
					print (bcolors.red + 'ERROR: ' + pathToFile + ' does not exist, maybe running the ntuple production helps' + bcolors.reset)
					continue
				print ('Obtaining 1D efficiency plots for ' + fil + ' with quality type ' + quality + ' and chamber ' + ch)
				for var in plotting1D.keys():
					effplot1D.plot1DEff(pathToFile, outputPlotDir + quality + "_" + acctxt.replace(' ', '_'), plotting1D[var], "Eff" + var, ch, "RECO", effPath, fil, quality)
	return

if plot:
	Plot1DEffs(files, qualities, plotting1D, effPath, outputPath)
	Plot2DEffs(files, qualities, plotting, effPath, outputPath)

