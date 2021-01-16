#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################
# main basf2 analysis code
# author: Seema Choudhury
#########################
import os
import datetime
import time
from modularAnalysis import *
import basf2 as b2
import sys
from stdPhotons import stdPhotons
import flavorTagger as ft
from variables.collections import *
from variables.utils import *
import ROOT as root 
import vertex
from variables import variables
import get_variables as gv
import variables as va
from stdV0s import *
weightfiles = 'B2nunubarBGx1'
b2.conditions.prepend_globaltag("analysis_tools_release-04")
# colors
blue = '\033[1;34m'
red = '\033[1;31m'
end = '\033[0m'
purple = '\033[1;35m'
yellow = '\033[1;33m'
cyan = "\033[0;36m"


def main(input_file, output_file, channel=None,  isdata=None):
    path = makeFSPs(input_file, isdata)
    if channel == 'kst0mumu':
        gv.Kst0mumu(path)
        ft.flavorTagger(particleLists=['B0:kst0mumu'],
                        belleOrBelle2="Belle2", weightFiles=weightfiles, path=path)
        varKst0mumu(path, output_file)
    elif channel == 'kst0ee':
        Kst0ee(path)
        ft.flavorTagger(particleLists=['B0:kst0ee'],
                        belleOrBelle2="Belle2", weightFiles=weightfiles, path=path)
        varKst0ee(path, output_file)
        varKst0ee_brem(path, output_file)
    elif channel == 'kst_ksee':
        Kstee_ks(path)
        ft.flavorTagger(particleLists=['B+:kst_ksee'],
                        belleOrBelle2="Belle2", weightFiles=weightfiles, path=path)
        varKstee_ks(path, output_file)
        varKstee_ks_brem(path, output_file)
    elif channel == 'kst_pi0ee':
        Kstee_pi0(path)
        ft.flavorTagger(particleLists=['B+:kst_pi0ee'],
                        belleOrBelle2="Belle2", weightFiles=weightfiles, path=path)
        varKstee_pi0(path, output_file)
        varKstee_pi0_brem(path, output_file)
    elif channel == 'kst_ksmumu':
        Kstmumu_ks(path)
        ft.flavorTagger(particleLists=['B+:kst_ksmumu'],
                        belleOrBelle2="Belle2", weightFiles=weightfiles, path=path)
        varKstmumu_ks(path, output_file)
    elif channel == 'kst_pi0mumu':
        Kstmumu_pi0(path)
        ft.flavorTagger(particleLists=['B+:kst_pi0mumu'],
                        belleOrBelle2="Belle2", weightFiles=weightfiles, path=path)
        varKstmumu_pi0(path, output_file)
        varKst0mumu(path, output_file)
        varKstmumu_ks(path, output_file)
        varKstmumu_pi0(path, output_file)
        varKst0ee(path, output_file)
        varKstee_ks(path, output_file)
        varKstee_pi0(path, output_file)
    progress = register_module('Progress')
    path.add_module(progress)
    b2.process(path)
    # print out the summary
    print(b2.statistics)
# NOTE for e channel vertexTree giving an error on release-03-02-04
# is not there but stil some error


def makeFSPs(input_file,  isdata):
    path = b2.create_path()
    inputMdstList(environmentType='default', filelist=input_file, path=path)
    printDataStore(path=path)
    # path.add_module('EventInfoPrinter')
    variables.addAlias('pi_k', 'pidPairProbabilityExpert(211, 321, ALL)')
    variables.addAlias('k_pi', 'pidPairProbabilityExpert(321, 211, ALL)')
    # electron
    impactcut = 'abs(d0)<2.0 and abs(z0)<5.0 '
    # NOTE For BtoXll skim
    # event level cuts: R2 and require a minimum number of tracks
    fillParticleList(decayString='pi+:eventShapeForSkims', cut='pt > 0.1', path=path)
    fillParticleList(decayString='gamma:eventShapeForSkims', cut='E > 0.1', path=path)
    buildEventShape(inputListNames=['pi+:eventShapeForSkims', 'gamma:eventShapeForSkims'],
                    allMoments=False,
                    foxWolfram=True,
                    harmonicMoments=False,
                    cleoCones=False,
                    thrust=False,
                    collisionAxis=False,
                    jets=False,
                    sphericity=False,
                    checkForDuplicates=False,
                    path=path)

    fillParticleList(decayString='e+:all',
                     cut=impactcut + ' and electronID > 0.5', path=path)
    # muon
    fillParticleList(decayString='mu+:all',
                     cut=impactcut + ' and muonID > 0.5', path=path)
    # pion
    fillParticleList(decayString='pi+:all',
                     cut=impactcut + ' and pi_k>0.1', path=path)
    # kaon
    fillParticleList(decayString='K+:all',
                     cut=impactcut + ' and k_pi>0.1', path=path)
    # gamma
    fillParticleList(decayString='gamma:all', cut='', path=path)
    if isdata:
        print(yellow, 'Momentum correction for FSPs', end)
        trackingMomentum(inputListNames=['mu+:all', 'e+:all', 'K+:all', 'pi+:all'], scale=1.00056, path=path)
    # Bremrecovery
    correctBremsBelle('e-:corrected', 'e-:all', 'gamma:all', True, 0.05, 0.05,
                      False, path=path)
    # Ks0
    stdKshorts( path=path)
    cutAndCopyList('K_S0:my_ks', 'K_S0:merged', 'goodBelleKshort==1', writeOut=False, path=path)
    matchMCTruth(list_name='K_S0:my_ks', path=path)
    # pi0
    stdPhotons('all', path=path)
    cutAndCopyList(
        'gamma:eff30_Jan2020',
        'gamma:all',
        '[clusterReg==1 and E>0.080] or [clusterReg==2 and E>0.030] or [clusterReg==3 and E>0.060]',
        path=path)
    reconstructDecay('pi0:all -> gamma:eff30_Jan2020 gamma:eff30_Jan2020',
                     '0.1215 < M < 0.1415 and E > 0.4 and abs(cosHelicityAngleMomentum)<0.8', path=path)
    #vertex.KFit('pi0:all', conf_level=0.0, path=path)
    vertex.KFit('pi0:all', conf_level=0.0,fit_type='mass', path=path)
    # K*0->K+ pi-
    reconstructDecay(decayString='K*0:all -> K+:all pi-:all',
                     cut='0.6 < M < 1.2', path=path)
    matchMCTruth(list_name='K*0:all', path=path)
    # K*+ -> K_s pi+
    reconstructDecay(decayString='K*+:ks -> K_S0:my_ks pi+:all',
                     cut='0.6 < M < 1.2', path=path)
    matchMCTruth(list_name='K*+:ks', path=path)
    # K*+ -> K+ pi0
    reconstructDecay(decayString='K*+:pi0 -> K+:all pi0:all',
                     cut='0.6 < M < 1.2', path=path)
    matchMCTruth(list_name='K*+:pi0', path=path)
    return path


def continuum(particlelist, path):
    buildRestOfEvent(target_list_name=particlelist, path=path)
    roename = 'cleanMask' + particlelist
    cleanMask = (roename, 'nCDCHits > 0 and useCMSFrame(p)<=10', 'p >= 0.05 and useCMSFrame(p)<=10')
    appendROEMasks(list_name=particlelist, mask_tuples=[cleanMask], path=path)
    buildContinuumSuppression(particlelist, roename, path=path)


def reconstruct_B(path, lepton, kstar):
    # lepton
    if lepton == 'e':
        lep_list = ' e+:corrected e-:corrected'
        lepbrem_list = ' e+:all e-:all'
    elif lepton == 'mu':
        lep_list = ' mu+:all mu-:all'
    # K*
    if kstar == 'kst0':
        charge = '0'
        kstlist = 'K*0:all'
    elif kstar == 'kst_ks':
        charge = '+'
        kstlist = 'K*+:ks'
    elif kstar == 'kst_pi0':
        charge = '+'
        kstlist = 'K*+:pi0'
    Bmeson = 'B'+charge+':'+kstar+lepton+lepton
    chain = Bmeson+' -> ' + kstlist + lep_list
    print(purple, 'reconstruct decay:', end, cyan, chain, end)
    # reconstruction
    reconstructDecay(decayString=chain,
                     cut='5.2 < Mbc < 5.29 and -0.3 < deltaE < 0.3', path=path)
    matchMCTruth(list_name=Bmeson, path=path)


    
    # ROE and vfit
    continuum(Bmeson, path)
    Apply_VerTexFit(path, lepton, kstar)
    vertex.TagV(list_name=Bmeson, maskName='cleanMask'+Bmeson, path=path)

    if lepton == 'e':
        Bmeson = 'B'+charge+':'+kstar+lepton+lepton + '_brem'
        chain = Bmeson+' -> ' + kstlist + lep_list + lepbrem_list
        print(purple, 'reconstruct decay:', end, cyan, chain, end)
        reconstructDecay(decayString=chain,
                         cut='5.2 < Mbc < 5.3 and -0.3 < deltaE < 0.3', path=path)


def Apply_VerTexFit(path, lepton, kstar):
    if lepton == 'e':
        lep_list = ' ^e+:corrected ^e-:corrected'
    elif lepton == 'mu':
        lep_list = ' ^mu+:all ^mu-:all'
    if kstar == 'kst0':
        kstlist = ' [K*0:all -> ^K+:all ^pi-:all] '
        charge = '0'
    elif kstar == 'kst_ks':
        kstlist = ' [K*+:ks -> ^K_S0:my_ks ^pi+:all] '
        charge = '+'
    elif kstar == 'kst_pi0':
        kstlist = ' [K*+:pi0 -> ^K+:all ^pi0:all] '
        charge = '+'
  #  Jpsi = ' J/psi' + ':' +lepton+lepton
    Bmeson = 'B'+charge+':'+kstar+lepton+lepton
    chain = Bmeson+' -> ' + kstlist + lep_list
    print(purple, 'vfit:  ', cyan, chain, end)
    # perform vertex fit of J/psi candidates
   # vertex.KFit(Jpsi, conf_level=0.0, path=path)
    vertex.KFit(list_name=Bmeson,fit_type='vertex',decay_string=chain, conf_level=0.00, path=path)







####
# NOTE
# if we use argparse inside basf2 steeringfile then we have to use like --arg in inputfile.root(as eg)
# input file can be obtained by -i that is made by basf2 framework. Still we prepaare -in to get a default
# option of those arguments
# that is --arg -in can be replaced by -i
# but --arg -out can't be replaced by -o
###


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    import argparse
    parser1 = argparse.ArgumentParser()
    parser1.add_argument('-c', '--channel',
                         help='which channel?   ', default='all')
    parser1.add_argument('-in', '--inp', help='input file',
                         default=b2.find_file('/group/belle2/users/tenchini/prerelease-05-00-00a/1111540100/1111540100_eph3_BGx0_29.root'))
    parser1.add_argument('-out', '--out', help='outputfile', default='ee.root')
    parser1.add_argument('-isdata', '--isdata', help='whether data oe mc. For data\
                           momentum corretion for FSP wil take place', default=False)
    parser1.print_help()
    args = parser1.parse_args()

    date = datetime.datetime.today()
    print(date.strftime('Start at : %d-%m-%y %H:%M:%S\n'))
    stime = time.time()
    print(args.inp, args.out, args.channel)
    main(args.inp, args.out, args.channel, args.isdata)
    date = datetime.datetime.today()
    print(date.strftime('End at : %y-%m-%d %H:%M:%S\n'))
    print('time taken', time.time()-stime)


