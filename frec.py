#!/usr/bin/env python3

#####################################################################
#
#
# This is a template for the US Belle II Summer School 2020
# It is intended as a starting point for an analysis of
#
# B0 -> J/psi Ks
#        |    |
#        |    +-> pi+ pi-
#        |
#        +-> mu- mu+
#
#
#####################################################################

import sys
import basf2 as b2
import modularAnalysis as ma
import stdV0s
import flavorTagger as ft
import vertex
from variables import variables
import variables.collections as vc
import variables.utils as vu

# filenumber = sys.argv[1]

# set analysis global tag
b2.use_central_database("analysis_tools_release-04-02")

# create path
main = b2.Path()

# load input data from mdst/udst file
filedirectory = '/gpfs/group/belle2/users/seemac/Kstll/signal/BtoKstjpsi/kst0jpsi'
ma.inputMdstList(environmentType='default', filelist=[f'{filedirectory}/mdst_000001_prod00012871_task10020000001.root'], path=main)


# fill final state particle lists
impactcut = 'abs(d0)<2.0 and abs(z0)<5.0 '
ma.fillParticleList(decayString='mu+:all',
                     cut=impactcut + ' and muonID > 0.5', path=main)
stdV0s.stdKshorts(path=main)



# combine final state particles to form composite particles
ma.reconstructDecay('J/psi:mumu -> mu+:all mu-:all', cut='dM < 0.11', path=main)

# perform vertex fit of J/psi candidates
vertex.KFit('J/psi:mumu', conf_level=0.0, path=main)

# combine J/psi and KS candidates to form B0 candidates
ma.reconstructDecay('B0 -> J/psi:mumu K_S0:merged', cut='Mbc > 5.2 and abs(deltaE) < 0.3', path=main)

# match reconstructed with MC particles
ma.matchMCTruth('B0', path=main)

# build the rest of the event
ma.buildRestOfEvent('B0', fillWithMostLikely=True, path=main)

# call flavor tagging
ft.flavorTagger('B0', path=main)

# remove B0 candidates without a valid flavor information
ma.applyCuts('B0', 'qrOutput(FBDT) > -2', path=main)

# fit B vertex on the tag-side
vertex.TagV('B0', constraintType='tube', fitAlgorithm='Rave', path=main)

# perform best candidate selection
#b2.set_random_seed('USBelleIISummerSchool')
#ma.rankByHighest('B0', variable='random', numBest=1, path=main)

# create list of variables for output ntuple
standard_vars = vc.kinematics + vc.mc_kinematics + vc.mc_truth
fs_vars = vc.pid + vc.track + vc.track_hits + standard_vars
jpsi_ks_vars = vc.inv_mass + vc.vertex + vc.mc_vertex + standard_vars
b_vars = vc.deltae_mbc + vc.tag_vertex + vc.mc_tag_vertex + ft.flavor_tagging + standard_vars

#b_vars += vu.create_aliases_for_selected([*fs_vars, 'isBremsCorrected'], 'B0 -> [J/psi -> ^e+ ^e-] K_S0', prefix=['ep', 'em'])
b_vars += vu.create_aliases_for_selected(fs_vars, 'B0 -> J/psi [K_S0 -> ^pi+ ^pi-]', prefix=['pip', 'pim'])
b_vars += vu.create_aliases_for_selected(jpsi_ks_vars, 'B0 -> ^J/psi ^K_S0')

cmskinematics = vu.create_aliases(vc.kinematics, 'useCMSFrame({variable})', 'CMS')
b_vars += vu.create_aliases_for_selected(cmskinematics, '^B0 -> [^J/psi -> ^mu+ ^mu-] [^K_S0 -> ^pi+ ^pi-]')

#variables.addAlias('withBremsCorrection', 'passesCut(passesCut(ep_isBremsCorrected == 1) or passesCut(em_isBremsCorrected == 1))')
#b_vars += ['withBremsCorrection']

# save variables to an ntuple
ma.variablesToNtuple('B0', variables=b_vars, filename='Bd2JpsiKS.root', treename='tree', path=main)

# process the events
b2.process(main)

# print out the summary
print(b2.statistics)


