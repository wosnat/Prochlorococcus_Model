#!/usr/bin/env python
# coding: utf-8

# # N starvation 
# 
# ## Which N sources allow MED4 to grow

import pandas as pd
import matplotlib.pyplot as plt
import cobra
import numpy as np
import seaborn as sns
import itertools
from matplotlib.colors import LogNorm, Normalize
import os

from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import production_envelope
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

# values from Ofaim at el.

PARAMETER_VALUES = [#"Name",     "Reaction ID",          "Lower bound", "UpperBound"
                   ["HCO3",      "HCO3EXcar",            -8,            0],
                   #["Nitrogen",  "AmmoniaEX",            -0.56,         0],
                   #["Phosphate", "FAKEOrthophosphateEX", -0.1,          0],
                   ["Light",     "LightEX",              -150,          0]]
CO2MAX = 0.82

# Block fake reactions
FAKE_TRANSPORT = ["AminosugarsOUT", "FAKEAAOUT", "FAKEABPOUT", "FAKEacpTRANS", "FAKEApoacpTRANS", "FAKEThioredoxinTRANS", "FreefattyacidsOUT", "7NMeth7carbOUT", "ArtificialproteinOUT", "FADOUT", "LipoylproteinTRANS", "MenaquinoneOUT", "NicotinateOUT", "THFpolyglutOUT", "Thiamin_dpOUT"]


# # Import model and manipulate based on Ofaim at el
def load_model(remove_blocked = False):
    model_dpath = os.path.join('..', 'Model_files')
    model_fname = 'iSO595v7.xml'
    model_fpath = os.path.join(model_dpath, model_fname)
    model = cobra.io.read_sbml_model(model_fpath)


    # manipulations copied from Ofaim at el.

    # Block H2S
    model.reactions.H2SEX.lower_bound = 0

    # Block fake transports
    for rid in FAKE_TRANSPORT:
        model.reactions.get_by_id(rid).bounds = (0,0)

    # Remove blocked reactions
    if remove_blocked:
        blocked = cobra.flux_analysis.find_blocked_reactions(model, open_exchanges = True)
        print('blocked', len(blocked), blocked)
        model.remove_reactions([model.reactions.get_by_id(r_id) for r_id in blocked])

    # Block maximum CO2 production
    model.reactions.CO2EX.bounds = (0, CO2MAX)

    for i, row in enumerate(PARAMETER_VALUES):
        # Row: Name, Reaction ID, lower bound, upper bound
        key = row[0]
        reaction_id = row[1]
        lower_bound = row[2]
        upper_bound = row[3]
        r = model.reactions.get_by_id(reaction_id)
        # Fix flux
        r.bounds = (lower_bound, upper_bound)
    return model
    


#exchange_ids = [r.id for r in model.exchanges] # + ["R00024"]

# # Identify N sources where MED4 grows on

#double_gene_deletion()

def _gene_deletion_feeding(model, uptake1, sec1, maxflux):
    with model:
        medium = model.medium
        medium["AmmoniaEX"] = 0.0
        medium[uptake1] = 1000.0
        model.medium = medium
        model.reactions.get_by_id(uptake1).upper_bound = -1e-5
        model.reactions.get_by_id(sec1).upper_bound = maxflux
        model.reactions.get_by_id(sec1).lower_bound = maxflux
        #model.reactions.BIOMASS.upper_bound = 1e-2
        #model
        solution = model.optimize()
        print(uptake1, sec1, model.summary())
        df = single_gene_deletion(model, processes=10)
        df['uptake'] = uptake1
        df['secretion'] = sec1
        df = df.loc[df.status != 'optimal']
        print(df.head())
        return df




if __name__ == '__main__':
    import argparse

    # def generate_json_and_run_from_X(X,  json_dpath, out_dpath, out_fprefix, timeout=10*60):
    parser = argparse.ArgumentParser(description='run gene knockouts to find uptake/secretion pathways.')
    parser.add_argument('--uptake', help='uptake reaction', required=True)
    parser.add_argument('--secretion', help='secretion reaction', required=True)
    parser.add_argument('--bound', help='secretion bound', required=True, type=float)
    parser.add_argument("--out_dpath", help="output dir", default='.')
    
    args = parser.parse_args()
    dpath = args.out_dpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
    else:
        dpath = '.'
    out_fpath = os.path.join(dpath, f"secretion_knockout_1gene_{args.uptake}_{args.secretion}.csv")

    model = load_model()
    resdf = _gene_deletion_feeding(model, args.uptake, args.secretion, args.bound)
    resdf.to_csv(out_fpath)



