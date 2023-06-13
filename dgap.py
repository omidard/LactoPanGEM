#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:27:09 2022

@author: omidard
"""

#diverse template gapfilling: dgap

#imports


import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
import matplotlib.pyplot as plt
from cobra import Model, Reaction, Metabolite
import multiprocessing as mp
import numpy as np
import os




#formulation of a chemically defined media
def m9(model):
    model.solver = 'gurobi'
    for reaction in model.reactions:
        if 'EX_' in  reaction.id:
            reaction.lower_bound=0

    #amino acids
    model.reactions.EX_arg_L_e.lower_bound = -2
    model.reactions.EX_cys_L_e.lower_bound = -2
    model.reactions.EX_glu_L_e.lower_bound = -2
    model.reactions.EX_ile_L_e.lower_bound = -2
    model.reactions.EX_leu_L_e.lower_bound = -2
    model.reactions.EX_met_L_e.lower_bound = -2
    model.reactions.EX_tyr_L_e.lower_bound = -2
    model.reactions.EX_phe_L_e.lower_bound = -2
    model.reactions.EX_thr_L_e.lower_bound = -2
    model.reactions.EX_val_L_e.lower_bound = -2
    model.reactions.EX_gly_e.lower_bound = -2
    model.reactions.EX_ala_L_e.lower_bound = -2
    model.reactions.EX_asp_L_e.lower_bound = -2
    model.reactions.EX_his_L_e.lower_bound = -2
    model.reactions.EX_lys_L_e.lower_bound = -2
    model.reactions.EX_pro_L_e.lower_bound = -2
    model.reactions.EX_ser_L_e.lower_bound = -2
    model.reactions.EX_trp_L_e.lower_bound = -2


    #carbon
    model.reactions.EX_glc_D_e.lower_bound = -15
    #model.reactions.EX_glc_D_e.upper_bound = -5
    model.reactions.EX_ac_e.lower_bound = -1
    model.reactions.EX_cit_e.lower_bound = -1

    #nuclotides
    model.reactions.EX_thymd_e.lower_bound = -0.1
    model.reactions.EX_ura_e.lower_bound = -1
    model.reactions.EX_gua_e.lower_bound = -1
    model.reactions.EX_ins_e.lower_bound = -1
    model.reactions.EX_ade_e.lower_bound = -1
    model.reactions.EX_xan_e.lower_bound = -1
    model.reactions.EX_orot_e.lower_bound = -1
    
    #vitamins
    model.reactions.EX_btn_e.lower_bound = -0.1
    model.reactions.EX_pnto_R_e.lower_bound = -0.5
    model.reactions.EX_thm_e.lower_bound = -0.1
    model.reactions.EX_pydam_e.lower_bound = -0.1
    model.reactions.EX_pydxn_e.lower_bound = -0.1
    model.reactions.EX_ribflv_e.lower_bound = -0.1
    model.reactions.EX_fol_e.lower_bound = -0.1
    model.reactions.EX_ascb_L_e.lower_bound = -0.5
    model.reactions.EX_4abz_e.lower_bound = -1
    model.reactions.EX_nac_e.lower_bound = -1
    
    #minerals
    #model.reactions.EX_o2_e.lower_bound = -10
    model.reactions.EX_cl_e.lower_bound = -10
    model.reactions.EX_h_e.lower_bound = -1000
    model.reactions.EX_h2o_e.lower_bound = -10
    model.reactions.EX_nh4_e.lower_bound = -10
    model.reactions.EX_ca2_e.lower_bound = -10
    model.reactions.EX_co_e.lower_bound = -10
    model.reactions.EX_co2_e.lower_bound = -10
    model.reactions.EX_pi_e.lower_bound = -100
    model.reactions.EX_cobalt2_e.lower_bound = -10
    model.reactions.EX_cu2_e.lower_bound = -10
    model.reactions.EX_fe3_e.lower_bound = -10
    model.reactions.EX_k_e.lower_bound = -10
    model.reactions.EX_mn2_e.lower_bound = -10
    model.reactions.EX_so4_e.lower_bound = -10
    model.reactions.EX_na1_e.lower_bound = -10
    model.reactions.EX_mg2_e.lower_bound = -10
    model.reactions.EX_zn2_e.lower_bound = -10
   
    #correction bounds HCLTr
    for reaction in model.reactions:
        if 'HCLTr' in  reaction.id:
            reaction.lower_bound=-10
            reaction.upper_bound=10
        if 'PTRCt2' in  reaction.id:
            reaction.lower_bound=-10
    #model.reactions.MNLpts.lower_bound = -1
    model.reactions.ATPM.upper_bound = 1
    model.reactions.ATPM.lower_bound = 1 #0.38
    model.reactions.ATPM.upper_bound = 1 #0.38
    model.reactions.EX_ptrc_e.upper_bound = 0
    model.reactions.EX_ptrc_e.lower_bound = 0
    #model.reactions.PTRCORNt7.bounds=(-10,10)
    
    
    #model.reactions.BIOMASS2.upper_bound = 1000
    #model.reactions.BIOMASS2.lower_bound = 0
    model.objective = 'BIOMASS2'
    return model
    


#scanning models and retrive important information
def scan(directory):
    models =glob('%s/*.json'%directory)
    gr=[]
    mid=[]
    reactions_total = []
    reactions_gap = []
    failed=[] #failed models
    temps=[] #gapfilled mopdels
    for mod in models:
        model=load_json_model(mod)
        m9(model)
        fba=model.optimize()
        failed.append(model.id)
        if fba.fluxes.BIOMASS2 >= 0.01:
            temps.append(model.id)
        gr.append(fba.fluxes.BIOMASS2)
        mid.append(model.id)
        reactions_gap.append(len(model.genes.GAP.reactions))
        reactions_total.append(len(model.reactions))
    all_gems = pd.DataFrame()
    all_gems['id']=mid
    all_gems['growth']=gr
    all_gems['total_reactions']=reactions_total
    all_gems['gapfilled_reactions']=reactions_gap
    return all_gems,failed,temps
 
        
 
#find best template
def tempfind(all_gems,failed,temps):
    tr=[]
    gapr=[]
    for i in range(len(all_gems)):
        if all_gems.id[i] in temps:
            tr.append(all_gems.total_reactions[i])
            gapr.append(all_gems.gapfilled_reactions[i])
    temps_inf = pd.DataFrame()
    temps_inf['id'] = temps
    temps_inf['total_reactions'] = tr
    temps_inf['gapfilled_reactions']=gapr
    temps_inf.sort_values(by='total_reactions', axis=0, ascending=False, inplace=True, kind='quicksort', na_position='last', ignore_index=True, key=None)
    return temps_inf



#find candidate reactions
def gaps(failed,temp_name):
    temp_re_id=[] #list1
    dir2='/home/omidard/allgems/Limosilactobacillus_fermentum/biomassed'
    template = load_json_model(dir2+'/'+temp_name)
    for reaction in template.reactions:
        temp_re_id.append(reaction.id)
    allmissing=[]
    for i in failed:
        missing=pd.DataFrame()
        model = load_json_model(dir2+'/'+i)
        failed_re_id=[] #list two
        for reaction in model.reactions:
            failed_re_id.append(reaction.id)
        missing_reactions = list(set(temp_re_id).difference(failed_re_id))
        missing[model.id] = missing_reactions
        allmissing.append(missing)
    return allmissing



#alternate function development
def zx(allmissing):
    xlist=[] #potential gaps
    modelsz=[] #models
    for i in allmissing:
        for x in i:
            modelsz.append(x)
            zlist=[]
            for z in i[x]:
                zlist.append(z)
            xlist.append(zlist)
    return xlist,modelsz



#add candidate reactions and save models based on gapfilling status
def addgaps(allgapz,modelsz,temp_name):
    for i in range(len(modelsz)):
        dir2='/home/omidard/allgems/Limosilactobacillus_fermentum/biomassed'
        dir3='/home/omidard/allgems/Limosilactobacillus_fermentum/feasible/'
        dir4='/home/omidard/allgems/Limosilactobacillus_fermentum/failed2/'
        model = load_json_model(dir2+'/'+modelsz[i])
        template = load_json_model(dir2+'/'+temp_name)
        for x in allgapz[i]:
            if x != 'no':
                reaction = template.reactions.get_by_id(x)
                reaction2 = Reaction(reaction.id)
                reaction2.name = reaction.name
                reaction2.subsystem = reaction.subsystem
                reaction2.lower_bound = reaction.lower_bound
                reaction2.upper_bound = reaction.upper_bound
                reaction2.add_metabolites(reaction.metabolites)
                reaction2.gene_reaction_rule = '(GAP)'
                model.add_reactions([reaction2])
                model.repair()
                print(reaction2.id)
                print('-----done')
                m9(model)
                print(model.optimize().fluxes.BIOMASS2)
            if model.optimize().fluxes.BIOMASS2>0.01:
                cobra.io.json.save_json_model(model,dir3+model.id)
            if model.optimize().fluxes.BIOMASS2<=0.01:
                cobra.io.json.save_json_model(model,dir4+model.id)



def fluxanalyze (rea):
    dir2='/home/omidard/allgems/Limosilactobacillus_fermentum/biomassed'
    temp_name = 'GCF_009556455.1.json'
    template = load_json_model(dir2+'/'+temp_name)
    m9(template)
    template.reactions.get_by_id(rea).lower_bound = 0
    template.reactions.get_by_id(rea).upper_bound = 0
    if template.optimize().fluxes.BIOMASS2 < 0.1:
        tar = rea
    else:
        tar = 'no'
    return tar




def gap_extract(mod):
    model=load_json_model(mod)
    total_gaps =pd.DataFrame({model.id:[re.id for re in model.genes.GAP.reactions]})
    tg = total_gaps.T
    return tg

#find non essential gaps
def none_ess_gap(total_gaps,directory1,directory2):
    model = load_json_model(directory1+total_gaps.columns[0])
    print('in process',model.id)
    m9(model)
    for re in total_gaps[total_gaps.columns[0]]:
        try:
            m9(model)
            low = model.reactions.get_by_id(re).lower_bound
            up = model.reactions.get_by_id(re).upper_bound
            model.reactions.get_by_id(re).lower_bound = 0
            model.reactions.get_by_id(re).upper_bound = 0
            if model.optimize().fluxes.BIOMASS2 > 0.1:
                reaction = model.reactions.get_by_id(re)
                reaction.remove_from_model()
                print(re,'removed')
            else:
                model.reactions.get_by_id(re).lower_bound = low
                model.reactions.get_by_id(re).upper_bound = up
        except AttributeError:
            pass
    cobra.io.json.save_json_model(model,directory2+model.id)
    print(directory2,'finished')



def total_dir(rootdir):
    dirs = []
    for file in os.listdir(rootdir):
        local = os.path.join(rootdir, file)
        if os.path.isdir(local):
            dirs.append(local)
    return dirs




def exchange_reactions(directory):
    models =glob('%s/*.json'%directory)
    model = load_json_model(models[0])
    model = m9(model)
    exchanges = []
    for reaction in model.reactions:
        if 'EX_' in reaction.id and model.reactions.get_by_id(reaction.id).lower_bound < 0:
            exchanges.append(reaction.id)
    df = pd.DataFrame({'metabolites': exchanges})
    df2 = df.set_index(df['metabolites'])
    substrates = df2.drop('metabolites', axis=1)
    return substrates,exchanges



def essential_substrates(models,substrates):
    for mod in models:
        model = load_json_model(mod)
        growth_rates =[]
        for i in substrates[1]:
            m9(model)
            model.reactions.get_by_id(i).lower_bound = 0
            try:
                growth_rates.append(model.optimize().fluxes.BIOMASS2)
            except:
                growth_rates.append(0)
        col=model.id
        name=col.replace('.json','')
        substrates[0][name]=growth_rates
    return substrates



def exchanges_fluxes(model):
    ex_flux=[]
    ex_reaction=[]
    m9(model)
    all_flux = model.optimize().fluxes
    for i in range(len(all_flux)):
        if 'EX_' in all_flux.index[i]:
            ex_flux.append(all_flux[i])
            ex_reaction.append(all_flux.index[i])
    flux_frame = pd.DataFrame()
    flux_frame['exchanges']=ex_reaction
    flux_frame[model.id]=ex_flux
    flux_frame.sort_values(by=['exchanges'],inplace=True, ignore_index=True)
    flux_frame.set_index('exchanges',inplace=True)
    return flux_frame


def knockout_fluxes(model):
    ex_flux=[]
    ex_reaction=[]
    m9(model)
    all_flux = model.optimize().fluxes
    for i in range(len(all_flux)):
        if 'EX_' in all_flux.index[i]:
            ex_flux.append(all_flux[i])
            ex_reaction.append(all_flux.index[i])
    flux_frame = pd.DataFrame()
    flux_frame['exchanges']=ex_reaction
    flux_frame[model.id]=ex_flux
    flux_frame.sort_values(by=['exchanges'],inplace=True, ignore_index=True)
    flux_frame.set_index('exchanges',inplace=True)
    return flux_frame


def ess_re(mod):
    ess_collector=[]
    model=load_json_model(mod)
    for re in model.reactions:
        if 'EX_' not in re.id:
            model = load_json_model(mod)
            m9(model)
            model.reactions.get_by_id(re.id).lower_bound = 0
            model.reactions.get_by_id(re.id).upper_bound = 0
            try:
                if model.optimize().fluxes.BIOMASS2 <= 0:
                    ess_collector.append(re.id)
            except AttributeError:
                ess_collector.append(re.id)
    ess_collector2=pd.DataFrame()
    ess_collector2[model.id]=ess_collector
    return ess_collector2



def allinf (modelid):
    aux= ['EX_orot_e', 'EX_co_e', 'EX_pydam_e', 'EX_pydxn_e', 'EX_h_e', 'EX_h2o_e', 'EX_4abz_e', 'EX_ac_e', 'EX_ade_e', 'EX_ala_L_e', 'EX_arg_L_e', 'EX_asp_L_e', 'EX_btn_e', 'EX_ca2_e', 'EX_cit_e', 'EX_cl_e', 'EX_co2_e', 'EX_cys_L_e', 'EX_fol_e', 'EX_glc_D_e', 'EX_glu_L_e', 'EX_gly_e', 'EX_gua_e', 'EX_his_L_e', 'EX_ile_L_e', 'EX_ascb_L_e', 'EX_ins_e', 'EX_k_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_met_L_e', 'EX_nh4_e', 'EX_phe_L_e', 'EX_pi_e', 'EX_pnto_R_e', 'EX_pro_L_e', 'EX_ribflv_e', 'EX_ser_L_e', 'EX_thm_e', 'EX_thr_L_e', 'EX_thymd_e', 'EX_trp_L_e', 'EX_tyr_L_e', 'EX_ura_e', 'EX_val_L_e', 'EX_xan_e']
    carbs=['EX_acgam_e', 'EX_gal_e', 'EX_2ddglcn_e', 'EX_cellb_e', 'EX_fru_e', 'EX_gam_e', 'EX_glc_D_e', 'EX_glcn_e', 'EX_lcts_e', 'EX_malt_e', 'EX_man_e', 'EX_raffin_e', 'EX_rib_D_e', 'EX_sbt_D_e', 'EX_sucr_e']
    exchanges = ['EX_ac_e', 'EX_actn_S_e', 'EX_ade_e', 'EX_adn_e', 'EX_akg_e', 'EX_ala_D_e', 'EX_ala_L_e', 'EX_amp_e', 'EX_arg_L_e', 'EX_ascb_L_e', 'EX_asp_L_e', 'EX_btn_e', 'EX_cit_e', 'EX_co2_e', 'EX_cys_L_e', 'EX_fol_e', 'EX_for_e', 'EX_fum_e', 'EX_gam_e', 'EX_glc_D_e', 'EX_glu_L_e', 'EX_gly_e', 'EX_glyc_e', 'EX_gua_e', 'EX_h2o_e', 'EX_h2s_e', 'EX_h_e', 'EX_hdca_e', 'EX_his_L_e', 'EX_hxan_e', 'EX_ile_L_e', 'EX_ins_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_mal_L_e', 'EX_man_e', 'EX_met_L_e', 'EX_mnl_e', 'EX_nh4_e', 'EX_orn_e', 'EX_orot_e', 'EX_oxa_e', 'EX_phe_L_e', 'EX_pi_e', 'EX_pnto_R_e', 'EX_pro_L_e', 'EX_ptrc_e', 'EX_pydx_e', 'EX_pydxn_e', 'EX_pyr_e', 'EX_ribflv_e', 'EX_ser_L_e', 'EX_succ_e', 'EX_thm_e', 'EX_thr_L_e', 'EX_thymd_e', 'EX_trp_L_e', 'EX_tyr_L_e', 'EX_ura_e', 'EX_uri_e', 'EX_val_L_e', 'EX_xan_e', 'EX_lipoate_e', 'EX_co_e', 'EX_dad_5_e', 'EX_pydam_e', 'EX_12ppd_S_e', 'EX_acgam_e', 'EX_acnam_e', 'EX_fuc_L_e', 'EX_gal_e', 'EX_lac_L_e', 'EX_lald_L_e', 'EX_12dgr180_e', 'EX_13ppd_e', 'EX_2ddglcn_e', 'EX_2dmmq8_e', 'EX_2obut_e', 'EX_34dhpha_e', 'EX_34dhphe_e', 'EX_3mop_e', 'EX_4abut_e', 'EX_4abz_e', 'EX_4hbz_e', 'EX_5htrp_e', 'EX_5mthf_e', 'EX_acald_e', 'EX_adocbl_e', 'EX_ahcys_e', 'EX_alltn_e', 'EX_arab_L_e', 'EX_arbt_e', 'EX_arsenb_e', 'EX_asn_L_e', 'EX_btd_RR_e', 'EX_butso3_e', 'EX_C02528_e', 'EX_ca2_e', 'EX_cd2_e', 'EX_cellb_e', 'EX_cgly_e', 'EX_chol_e', 'EX_cholate_e', 'EX_chols_e', 'EX_cl_e', 'EX_cobalt2_e', 'EX_csn_e', 'EX_ctbt_e', 'EX_cu2_e', 'EX_cytd_e', 'EX_dad_2_e', 'EX_dcyt_e', 'EX_ddca_e', 'EX_dextrin_e', 'EX_dgsn_e', 'EX_diact_e', 'EX_din_e', 'EX_dopa_e', 'EX_dpcoa_e', 'EX_drib_e', 'EX_etha_e', 'EX_ethso3_e', 'EX_etoh_e', 'EX_fe3_e', 'EX_fecrm_e', 'EX_fru_e', 'EX_galt_e', 'EX_galur_e', 'EX_gcald_e', 'EX_gchola_e', 'EX_glcn_e', 'EX_glcur_e', 'EX_gln_L_e', 'EX_glyb_e', 'EX_glyclt_e', 'EX_gsn_e', 'EX_gthox_e', 'EX_gthrd_e', 'EX_h2_e', 'EX_hexs_e', 'EX_hg2_e', 'EX_hista_e', 'EX_ind3ac_e', 'EX_indole_e', 'EX_inost_e', 'EX_k_e', 'EX_lac_D_e', 'EX_lcts_e', 'EX_Lcyst_e', 'EX_malt_e', 'EX_malthx_e', 'EX_malttr_e', 'EX_melib_e', 'EX_met_D_e', 'EX_metsox_S_L_e', 'EX_mn2_e', 'EX_mops_e', 'EX_mqn8_e', 'EX_mso3_e', 'EX_n2o_e', 'EX_nac_e', 'EX_ncam_e', 'EX_ni2_e', 'EX_nmn_e', 'EX_no_e', 'EX_no2_e', 'EX_no3_e', 'EX_o2_e', 'EX_ocdcea_e', 'EX_pb_e', 'EX_pheme_e', 'EX_ppa_e', 'EX_ppi_e', 'EX_q8_e', 'EX_raffin_e', 'EX_rib_D_e', 'EX_rmn_e', 'EX_salcn_e', 'EX_sbt_D_e', 'EX_ser_D_e', 'EX_sheme_e', 'EX_so4_e', 'EX_spmd_e', 'EX_srtn_e', 'EX_sucr_e', 'EX_sulfac_e', 'EX_taur_e', 'EX_tchola_e', 'EX_tdchola_e', 'EX_thf_e', 'EX_tre_e', 'EX_trypta_e', 'EX_tsul_e', 'EX_ttdca_e', 'EX_tym_e', 'EX_urea_e', 'EX_xyl_D_e', 'EX_na1_e', 'EX_alaasp_e', 'EX_alagln_e', 'EX_alaglu_e', 'EX_alagly_e', 'EX_alahis_e', 'EX_alaleu_e', 'EX_alathr_e', 'EX_crn_e', 'EX_glyasn_e', 'EX_glyasp_e', 'EX_glycys_e', 'EX_glygln_e', 'EX_glyglu_e', 'EX_glyleu_e', 'EX_glymet_e', 'EX_glyphe_e', 'EX_glypro_e', 'EX_glytyr_e', 'EX_isetac_e', 'EX_mantr_e', 'EX_metsox_R_L_e', 'EX_mg2_e', 'EX_Ser_Thr_e', 'EX_stys_e', 'EX_zn2_e', 'EX_dha_e', 'EX_12ppd_R_e', 'EX_4ahmmp_e']
    auxf = []
    carbsf=[]
    carbsp = []
    exf =[]
    wgr = []
    names = []
    ids = modelid.replace('/home/omidard/gems/allgems/','')
    model = load_json_model(modelid)
    m9(model)
    fx = model.optimize().fluxes
    wgr.append(model.optimize().fluxes.BIOMASS2)
    for e in exchanges:
        exf.append(fx[e])
        
    for c in carbs:
        model = load_json_model(modelid)
        m9(model)
        model.reactions.EX_glc_D_e.lower_bound = 0
        model.reactions.get_by_id(c).lower_bound = -15
        fxc = model.optimize().fluxes
        carbsf.append(model.optimize().fluxes.BIOMASS2)
        for e in exchanges:
            carbsp.append(fxc[e])
            name = e+'_'+'growth_on_'+c
            names.append(name)
            
    for a in aux:
        model = load_json_model(modelid)
        m9(model)
        model.reactions.get_by_id(a).lower_bound = 0
        auxf.append(model.optimize().fluxes.BIOMASS2)   
    
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    df3 = pd.DataFrame()
    df4 = pd.DataFrame()
    #df1 setup
    df1.index = ['growth_on_'+c for c in carbs]
    df1[ids.replace('.json','')] = carbsf
    df1 = df1.T
    #df2
    df2.index = ['eliminated_'+a for a in aux]
    df2[ids.replace('.json','')] = auxf
    df2 = df2.T
    #df3
    df3.index = [n for n in names]
    df3[ids.replace('.json','')] = carbsp
    df3 = df3.T
    #df4
    df4.index = [e for e in exchanges]
    df4[ids.replace('.json','')] = exf
    df4 = df4.T
    #sumup
    df = pd.concat([df1,df2,df3,df4],axis = 1)
    df['wild_type_growth_rate'] = wgr
    return df


def all_reactions(modelid):
    reactions_list=[]
    model = load_json_model(modelid)
    for re in model.reactions:
        reactions_list.append(re.id)
    df = pd.DataFrame()
    df['all_reactions'] = reactions_list
    return df




def gene_associated_reactions (modelid):
    model = load_json_model(modelid)
    nonmetgenes=[]
    metgenes=[]
    model = load_json_model(modelid)
    for re in model.genes.GAP.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.EXCHANGE.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.ORPHAN.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.DEMAND.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.BIOMASS.reactions:
        nonmetgenes.append(re.id)
    for re in model.genes.spontaneous.reactions:
        nonmetgenes.append(re.id)
    for re in model.reactions:
        if re.id not in nonmetgenes:
            if 'SINK' not in re.id:
                metgenes.append(re.id)
    
    df = pd.DataFrame()
    df[model.id]=metgenes
    df.set_index(model.id,inplace=True)
    df[model.id]=[1 for n in range(len(df.index))]
    
    df2 = pd.DataFrame()
    df2[model.id]=nonmetgenes
    df2.set_index(model.id,inplace=True)
    df2[model.id]=[0 for n in range(len(df2.index))]
    
    df3=pd.concat([df,df2],axis=0)
    
    return df3

    

def basic_counts(modelid):
    reactions=[]
    genes=[]
    gaps=[]
    model = load_json_model(modelid)
    reactions.append(len(model.reactions))
    genes.append(len(model.genes))
    gaps.append(len(model.genes.GAP.reactions))
    df = pd.DataFrame()
    df['reactions']=reactions
    df['genes']=genes
    df['gaps']=gaps
    df['id']=[model.id]
    df.set_index('id',inplace=True)
    return df


def all_fluxes(modelid):
    model = load_json_model(modelid)
    m9(model)
    reaction = []
    fx = model.optimize().fluxes
    for i in range(len(fx)):
        if fx[i] != 0:
            reaction.append(fx.index[i])
            
    dfx = pd.DataFrame()
    dfx['reactions']=reaction
    dfx[model.id] = [1 for v in range(len(dfx))]
    dfx.set_index('reactions',inplace=True)
    return dfx





def product_associated_genes(targetmet):
    allflux=pd.read_csv('/home/omidard/allflux.csv')
    allflux.set_index('reactions',inplace=True)
    producers=[]
    for c in allflux.columns:
        if allflux[c][targetmet] == 1:
            producers.append(allflux[c])
    allflux_producers_df = pd.concat(producers,axis=1)
    allflux_producers_df.fillna(0, inplace=True)
    return allflux_producers_df
            


def product_essential_genes(mod):
    allflux_producers_df=pd.read_csv('/home/omidard/allflux_producers_df.csv')
    allflux_producers_df.set_index('reactions',inplace=True)
    inv_genes_collect = []
    for i in allflux_producers_df.index:
        if allflux_producers_df[mod][i] == 1:
            df = pd.DataFrame()
            model = load_json_model('/home/omidard/gems/allgems/'+mod)
            m9(model)
            df['id'] = [model.id]
            df.set_index('id',inplace=True)
            df['gene']=[i]
            df['presence'] = [model.optimize().fluxes.EX_mnl_e]
            model.reactions.get_by_id(i).lower_bound = 0
            model.reactions.get_by_id(i).upper_bound = 0
            if model.optimize().fluxes.BIOMASS2 != 0:
                df['absence'] = [model.optimize().fluxes.EX_mnl_e]
            if model.optimize().fluxes.BIOMASS2 == 0:
                df['absence'] = ['essential gene']
            inv_genes_collect.append(df)
    gene_product_effect = pd.concat(inv_genes_collect,axis = 0)
    return gene_product_effect



            
def product_based_essential_genes(mod):
    targetmet = 'EX_4abut_e'  #change target met
    model = load_json_model(mod)
    m9(model)
    fluxed_reactions = []
    fx = model.optimize().fluxes
    if fx['EX_4abut_e'] !=0: #change target met
        for i in range(len(fx)):
            if fx[i] != 0:
                fluxed_reactions.append(fx.index[i])
    coll=[]
    if len(fluxed_reactions) > 0:
        for i in fluxed_reactions:
            model = load_json_model(mod)
            m9(model)
            df = pd.DataFrame()
            df['gene'] = [i]
            df.set_index('gene',inplace=True)
            df['wildtype'+model.id] = [model.optimize().fluxes.EX_4abut_e] #change target met
            model.reactions.get_by_id(i).lower_bound = 0
            model.reactions.get_by_id(i).upper_bound = 0
            if model.optimize().fluxes.BIOMASS2 !=0:
                df['knocked_out'+model.id] = [model.optimize().fluxes.EX_4abut_e] #change target met
            if model.optimize().fluxes.BIOMASS2 ==0:
                df['knocked_out'+model.id] = ['essential']
            coll.append(df)
    if len(fluxed_reactions) == 0:
        df = pd.DataFrame()
        coll.append(df)
    genes_effect = pd.concat(coll,axis=0)
    
    return genes_effect
                
                

                
def coreflux(mod):
    model = load_json_model(mod)
    m9(model)
    reaction=[]
    flux=[]
    flx = model.optimize().fluxes
    for i in flx.index:
        if flx[i] !=0:
            reaction.append(i)
            flux.append(flx[i])
    df = pd.DataFrame()
    df['reactions'] = reaction
    df[model.id] = flux
    df.set_index('reactions',inplace=True)
    return df       
    
        
        
def biomass2(model):
    model = load_json_model(model)
    reaction = Reaction('BIOMASS2')
    reaction.name = 'Biomass production'
    reaction.subsystem = 'Biomass'
    reaction.lower_bound = 0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default

    #metabolites
    thmpp_c = model.metabolites.get_by_id('thmpp_c')
    CPS_LBR_c = model.metabolites.get_by_id('CPS_LBR_c')
    DNA_LBR_c = model.metabolites.get_by_id('DNA_LBR_c')
    LIP_LBR_c = model.metabolites.get_by_id('LIP_LBR_c')
    LTAtotal_LBR_c = model.metabolites.get_by_id('LTAtotal_LBR_c')
    PGlac2_c = model.metabolites.get_by_id('PGlac2_c')
    PROT_LBR_c = model.metabolites.get_by_id('PROT_LBR_c')
    RNA_LBR_c = model.metabolites.get_by_id('RNA_LBR_c')
    atp_c = model.metabolites.get_by_id('atp_c')
    btn_c = model.metabolites.get_by_id('btn_c')
    coa_c = model.metabolites.get_by_id('coa_c')
    h2o_c = model.metabolites.get_by_id('h2o_c')
    nad_c = model.metabolites.get_by_id('nad_c')
    pydx5p_c = model.metabolites.get_by_id('pydx5p_c')
    thf_c = model.metabolites.get_by_id('thf_c')
    udcpdp_c = model.metabolites.get_by_id('udcpdp_c')
    adp_c = model.metabolites.get_by_id('adp_c')
    h_c = model.metabolites.get_by_id('h_c')
    pi_c = model.metabolites.get_by_id('pi_c')
    
    reaction.add_metabolites({
        thmpp_c: -0.0001,
        CPS_LBR_c: -0.078,
        DNA_LBR_c: -0.205,
        atp_c: -41.2,
        LIP_LBR_c: -0.106,
        LTAtotal_LBR_c: -0.006,
        PGlac2_c: -0.009,
        PROT_LBR_c: -3.311,
        RNA_LBR_c: -0.926,
        btn_c: -0.00001,
        coa_c: -0.002,
        h2o_c: -41.2,
        nad_c: -0.002,
        pydx5p_c: -0.000001,
        thf_c: -0.00001,
        udcpdp_c: -0.002,
        adp_c: 41.2,
        h_c: 41.2,
        pi_c: 41.2
    })
    reaction.gene_reaction_rule = '(BIOMASS)'
    reaction.reaction
    model.add_reactions([reaction])
    model.repair()
    
    cobra.io.json.save_json_model(model,'/home/omidard/biomassed/'+model.id)        
        
        
        
        
        
def producers(mod):
    producers=[]
    model = load_json_model(mod)
    m9(model)
    if model.optimize().fluxes.EX_mnl_e != 0: #targetmet
        producers.append(model.id)
    return producers

def pro_act_ge(mod):
    reactions=[]
    model = load_json_model('/home/omidard/gems/allgems/'+mod)
    m9(model)
    fx = model.optimize().fluxes
    for i in fx.index:
        if fx[i] !=0:
            reactions.append(i)
    df = pd.DataFrame()
    df[model.id] = reactions
    return df

def pro_ess_ge(df):
    non_ess = []
    for i in df[df.columns[0]]:
        model = load_json_model('/home/omidard/gems/allgems/'+df.columns[0])
        m9(model)
        model.reactions.get_by_id(i).lower_bound=0
        model.reactions.get_by_id(i).upper_bound=0
        if model.optimize().fluxes.BIOMASS2 != 0:
            non_ess.append(i)
    noess = pd.DataFrame()
    noess['index'] = non_ess
    noess.set_index('index',inplace=True)
    noess[df.columns[0]] = [1000 for v in range(len(noess))]
    return noess

def pro_ess_eff(df):
    flxc=[]
    for i in df.index:
        model = load_json_model('/home/omidard/gems/allgems/'+df.columns[0])
        m9(model)
        model.reactions.get_by_id(i).lower_bound=0
        model.reactions.get_by_id(i).upper_bound=0
        flx = model.optimize().fluxes.EX_mnl_e #targetmet
        flxc.append(flx)
    df[df.columns[0]] = flxc
    return df
        
        
        
    
        
            
    
    
    
            
    
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    







