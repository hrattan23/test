#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################
# variables to be put in ntuple
# author: Seema Choudhury
#########################
from modularAnalysis import *
from variables.collections import *
from variables.utils import *
import variables.collections as vc
import variables.utils as vu
from variables import variables
import ROOT as root

# Aliases
variables.addAlias('gengrandMotherPDG', 'genMotherPDG(1)')
variables.addAlias('genGrtgrandMotherPDG', 'genMotherPDG(2)')
variables.addAlias('genGrtGrtgrandMotherPDG', 'genMotherPDG(3)')
variables.addAlias('genGrtGrtGrtgrandMotherPDG', 'genMotherPDG(4)')
variables.addAlias('gengrandMotherID', 'genMotherID(1)')
variables.addAlias('genGrtgrandMotherID', 'genMotherID(2)')
variables.addAlias('genGrtGrtgrandMotherID', 'genMotherID(3)')
variables.addAlias('genGrtGrtGrtgrandMotherID', 'genMotherID(4)')
variables.addAlias('delta_z_ll',     'formula(daughter(1,dz)-daughter(2,dz))')
variables.addAlias('roeE_BKst0ee',    'useCMSFrame(roeE(cleanMaskB0:kst0ee))')
variables.addAlias('roeE_BKst0mumu',  'useCMSFrame(roeE(cleanMaskB0:kst0mumu))')
variables.addAlias('roeE_BKst1mumu',  'useCMSFrame(roeE(cleanMaskB+:kst_ksmumu))')
variables.addAlias('roeE_BKst1ee',    'useCMSFrame(roeE(cleanMaskB+:kst_ksee))')
variables.addAlias('roeE_BKst2mumu',  'useCMSFrame(roeE(cleanMaskB+:kst_pi0mumu))')
variables.addAlias('roeE_BKst2ee',    'useCMSFrame(roeE(cleanMaskB+:kst_pi0ee))')
variables.addAlias('d0_px',  'daughter(0,px)')
variables.addAlias('d0_py',  'daughter(0,py)')
variables.addAlias('d0_pz',  'daughter(0,pz)')
variables.addAlias('d0_p',  'daughter(0,p)')
variables.addAlias('d0_dr',  'daughter(0,dr)')
variables.addAlias('d0_dz',  'daughter(0,dz)')
variables.addAlias('d0_E',  'daughter(0,E)')
variables.addAlias('Mll',        'daughterCombination(M, 1, 2)')
variables.addAlias('Mll_before', 'daughterCombination(M, 1:0, 2:0)')
variables.addAlias('Mll_eeg',    'daughterCombination(M, 1:0, 2)')
variables.addAlias('Mll_ege',    'daughterCombination(M, 1, 2:0)')
variables.addAlias('qsquare',        'daughterCombination(M2, 1, 2)')
variables.addAlias('qsquare_before', 'daughterCombination(M2, 1:0, 2:0)')
variables.addAlias('qsquare_eeg',    'daughterCombination(M2, 1:0, 2)')
variables.addAlias('qsquare_ege',    'daughterCombination(M2, 1, 2:0)')
variables.addAlias('Ell', 'formula(daughter(1, useCMSFrame(E))+daughter(2, useCMSFrame(E)))')

# variable lists
cs_vars = ['foxWolframR2',
           'isNotContinuumEvent', 'qrOutput(FBDT)', 'useCMSFrame(cosTheta)', 'DeltaZ', 'R2',
           'thrustBm', 'thrustOm', 'cosTBTO', 'cosTBz', 'KSFWVariables(et)', 'KSFWVariables(mm2)',
           'KSFWVariables(hso00)', 'KSFWVariables(hso01)', 'KSFWVariables(hso02)', 'KSFWVariables(hso03)',
           'KSFWVariables(hso04)', 'KSFWVariables(hso10)', 'KSFWVariables(hso12)', 'KSFWVariables(hso14)',
           'KSFWVariables(hso20)', 'KSFWVariables(hso22)', 'KSFWVariables(hso24)', 'KSFWVariables(hoo0)',
           'KSFWVariables(hoo1)', 'KSFWVariables(hoo2)', 'KSFWVariables(hoo3)', 'KSFWVariables(hoo4)',
           'CleoConeCS(1)', 'CleoConeCS(2)', 'CleoConeCS(3)', 'CleoConeCS(4)', 'CleoConeCS(5)', 'CleoConeCS(6)',
           'CleoConeCS(7)', 'CleoConeCS(8)', 'CleoConeCS(9)',
           'CleoConeCS(1,ROE)', 'CleoConeCS(2,ROE)', 'CleoConeCS(3,ROE)', 'CleoConeCS(4,ROE)', 'CleoConeCS(5,ROE)',
           'CleoConeCS(6,ROE)', 'CleoConeCS(7,ROE)', 'CleoConeCS(8,ROE)', 'CleoConeCS(9,ROE)'
           ]
mctruth = ['mcPDG', 'genMotherPDG', 'genMotherID', 'mcErrors']
grand_mothers = ['mcPDG', 'genMotherPDG', 'genMotherID',
                 'gengrandMotherPDG', 'genGrtgrandMotherPDG', 'genGrtGrtgrandMotherPDG', 'genGrtGrtGrtgrandMotherPDG',
                 'gengrandMotherID', 'genGrtgrandMotherID', 'genGrtGrtgrandMotherID', 'genGrtGrtGrtgrandMotherID']
trackinfo = ['nCDCHits', 'nSVDHits', 'nPXDHits', 'dr', 'dz', 'charge', 'theta']
b_var = vc.deltae_mbc + mctruth + vc.kinematics + cs_vars +\
        ['delta_z_ll', 'chiProb', 'Mll', 'qsquare', 'isSignal', 'Ell', 'nTracks']
lepton_var = vc.kinematics + grand_mothers + trackinfo + ['electronID', 'muonID']
e_var = ['d0_px', 'd0_py', 'd0_pz', 'd0_E', 'd0_p', 'd0_dr', 'd0_dz']
b_evar = ['qsquare', 'qsquare_before', 'qsquare_eeg', 'qsquare_ege',
          'Mll', 'Mll_before', 'Mll_eeg', 'Mll_ege']
hadron_var = vc.kinematics + grand_mothers + trackinfo + ['k_pi', 'pi_k', 'kaonID', 'pionID']
gamma_var = vc.kinematics + grand_mothers + ['clusterReg']
conjugate_particle_var = vc.kinematics + ['chiProb', 'M', 'InvM'] + grand_mothers
pi0_var = ['daughterDiffOfPhi(0,1)', 'daughterAngleInBetween(0,1)', 'cosHelicityAngleMomentum']

kst0ll_var = b_var + ['cosHelicityAngle(0,0)'] +\
              vu.create_daughter_aliases(conjugate_particle_var, [0], 'B_Kst', False) + \
              vu.create_daughter_aliases(hadron_var, [0, 0], 'B_Kst_K', False) + \
              vu.create_daughter_aliases(hadron_var, [0, 1], 'B_Kst_pi', False) + \
              vu.create_daughter_aliases(lepton_var, [1], 'B_l1', False) + \
              vu.create_daughter_aliases(lepton_var, [2], 'B_l2', False)
Xee_var = b_evar +\
          vu.create_daughter_aliases(e_var, [1], 'B_l1', False) + \
          vu.create_daughter_aliases(e_var, [2], 'B_l2', False)
kstp_ksll_var = b_var + ['cosHelicityAngle(0,1)'] +\
              vu.create_daughter_aliases(conjugate_particle_var + ['cosHelicityAngle(0,0)'], [0], 'B_Kst', False) + \
              vu.create_daughter_aliases(conjugate_particle_var, [0, 0], 'B_Kst_Ks', False) + \
              vu.create_daughter_aliases(hadron_var, [0, 1], 'B_Kst_pi', False) + \
              vu.create_daughter_aliases(hadron_var, [0, 0, 0], 'B_Kst_Ks_pi1', False) + \
              vu.create_daughter_aliases(hadron_var, [0, 0, 1], 'B_Kst_Ks_pi2', False) + \
              vu.create_daughter_aliases(lepton_var, [1], 'B_l1', False) + \
              vu.create_daughter_aliases(lepton_var, [2], 'B_l2', False)
kstp_pi0ll_var = b_var + ['cosHelicityAngle(0,0)'] +\
              vu.create_daughter_aliases(conjugate_particle_var + ['cosHelicityAngle(1,0)'], [0], 'B_Kst', False) + \
              vu.create_daughter_aliases(hadron_var, [0, 0], 'B_Kst_K', False) + \
              vu.create_daughter_aliases(conjugate_particle_var + pi0_var, [0, 1], 'B_Kst_pi0', False) + \
              vu.create_daughter_aliases(conjugate_particle_var + gamma_var, [0, 1, 0], 'B_Kst_pi0_g1', False) + \
              vu.create_daughter_aliases(conjugate_particle_var + gamma_var, [0, 1, 1], 'B_Kst_pi0_g2', False) + \
              vu.create_daughter_aliases(lepton_var, [1], 'B_l1', False) + \
              vu.create_daughter_aliases(lepton_var, [2], 'B_l2', False)


def varKst0mumu(path,  output_file):
    variablesToNtuple(decayString='B0:kst0mumu',
                      variables=kst0ll_var + ['roeE_BKst0mumu'],
                      filename=output_file, treename='kst0mumu', path=path)


def varKst0ee(path,  output_file):
    variablesToNtuple(decayString='B0:kst0ee',
                      variables=kst0ll_var + Xee_var + ['roeE_BKst0ee'],
                      filename=output_file, treename='kst0ee', path=path)


def varKst0ee_brem(path,  output_file):
    variablesToNtuple(decayString='B0:kst0ee_brem',
                      variables=kst0ll_var,
                      filename=output_file, treename='kst0ee_brem', path=path)


def varKstmumu_ks(path,  output_file):
    variablesToNtuple(decayString='B+:kst_ksmumu',
                      variables=kstp_ksll_var + ['roeE_BKst1mumu'],
                      filename=output_file, treename='kst_ksmumu', path=path)


def varKstee_ks(path,  output_file):
    variablesToNtuple(decayString='B+:kst_ksee',
                      variables=kstp_ksll_var + ['roeE_BKst1ee'] + Xee_var,
                      filename=output_file, treename='kst_ksee', path=path)


def varKstee_ks_brem(path,  output_file):
    variablesToNtuple(decayString='B+:kst_ksee_brem',
                      variables=kstp_ksll_var,
                      filename=output_file, treename='kst_ksee_brem', path=path)


def varKstmumu_pi0(path,  output_file):
    variablesToNtuple(decayString='B+:kst_pi0mumu',
                      variables=kstp_pi0ll_var + ['roeE_BKst2mumu'],
                      filename=output_file, treename='kst_pi0mumu', path=path)


def varKstee_pi0(path,  output_file):
    variablesToNtuple(decayString='B+:kst_pi0ee',
                      variables=kstp_pi0ll_var + ['roeE_BKst2ee'] + Xee_var,
                      filename=output_file, treename='kst_pi0ee', path=path)


def varKstee_pi0_brem(path,  output_file):
    variablesToNtuple(decayString='B+:kst_pi0ee_brem',
                      variables=kstp_pi0ll_var,
                      filename=output_file, treename='kst_pi0ee_brem', path=path)
