from dgap import m9
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob

#load models information as a dataframe
all_inf=pd.read_csv('/home/omidard/all_inf.csv')
#add model path to the dataframe
wild_type_growth_rate=[]
for i in all_inf['file']:
    i = i+'.json'
    model = load_json_model(i)
    m9(model)
    growth = model.optimize().fluxes.BIOMASS2
    wild_type_growth_rate.append(growth)
all_inf['wild_type_growth_rate']=wild_type_growth_rate
all_inf.to_csv('/home/omidard/all_inf.csv')
#%%
#define a list of exchange reactions in the model
exchanges = ['EX_ac_e','EX_actn_S_e','EX_ade_e','EX_adn_e','EX_akg_e','EX_ala_D_e',
'EX_ala_L_e','EX_amp_e','EX_arg_L_e','EX_ascb_L_e','EX_asp_L_e','EX_btn_e','EX_cit_e',
'EX_co2_e','EX_cys_L_e','EX_fol_e','EX_for_e','EX_fum_e','EX_gam_e','EX_glc_D_e','EX_glu_L_e',
'EX_gly_e','EX_glyc_e','EX_gua_e','EX_h2o_e','EX_h2s_e','EX_h_e','EX_hdca_e','EX_his_L_e',
'EX_hxan_e','EX_ile_L_e','EX_ins_e','EX_leu_L_e','EX_lys_L_e','EX_mal_L_e','EX_man_e','EX_met_L_e',
'EX_mnl_e','EX_nh4_e','EX_orn_e','EX_orot_e','EX_oxa_e','EX_phe_L_e','EX_pi_e','EX_pnto_R_e',
'EX_pro_L_e','EX_ptrc_e','EX_pydx_e','EX_pydxn_e','EX_pyr_e','EX_ribflv_e','EX_ser_L_e','EX_succ_e',
'EX_thm_e','EX_thr_L_e','EX_thymd_e','EX_trp_L_e','EX_tyr_L_e','EX_ura_e','EX_uri_e','EX_val_L_e','EX_xan_e',]
#add exchange reactions to data frame
for i in exchanges:
    all_inf[i]=[0 for i in range(len(all_inf['file']))]
#collect production consumption rates for each exchange reaction per model
for i in range(len(all_inf['file'])):
    model = load_json_model(all_inf['file'][i]+'.json')
    m9(model)
    flx = model.optimize().fluxes
    for j in exchanges:
        if j in flx.index:
            all_inf[j][i]=flx[j]
all_inf.to_csv('/home/omidard/all_inf2.csv')



