"""
@author Ziyi Huang

@brief Functions for writing initial condition of MASCARET or REZO.
"""

import os
from argparse import ArgumentParser
import json
from datetime import datetime
import numpy as np
from shutil import copy2
from decimal import Decimal
import sys

from config import CFGS
from telapy.tools.study_masc_driven import MascaretStudy

def get_param():
    """
    Get configuration from JSON file
    @return param(dict) loaded configuration
    """

    parser = ArgumentParser(description = 'name of input json file')
    parser.add_argument('name_json', 
                        type=str, 
                        help='absolute directory of input json file')
    args = parser.parse_args()
    try:
        with open(args.name_json) as file:
            masc_ic_cnfg = json.load(file)
    except(ValueError, IOError):
        raise ValueError('incorrect input json file.\nSTOP.')

    param = {}

    if "system_config" in masc_ic_cnfg.keys():
        if "config_file" in masc_ic_cnfg['system_config'].keys():
            param['config_file'] = masc_ic_cnfg['system_config']['config_file']
        if "config_name" in masc_ic_cnfg['system_config'].keys():
            param['config_name'] = masc_ic_cnfg['system_config']['config_name']

    param['nom'] = list(masc_ic_cnfg['models_configs'].keys())[0]

    param['run_to_stdy'] = masc_ic_cnfg['run_configs']['run_ic']['run_to_stdy']
    if param['run_to_stdy']:
        param['prd_run_to_stdy'] = \
        masc_ic_cnfg['run_configs']['run_ic']['run_prd']
        param['tlrc_rltv_dpth_run_to_stdy'] = \
        masc_ic_cnfg['run_configs']['run_ic']['run_tlrc']['rltv_dpth']
        param['tlrc_abslt_dpth_run_to_stdy'] = \
        masc_ic_cnfg['run_configs']['run_ic']['run_tlrc']['abslt_dpth']
        param['tlrc_rltv_flrt_run_to_stdy'] = \
        masc_ic_cnfg['run_configs']['run_ic']['run_tlrc']['rltv_flrt']
        param['tlrc_abslt_flrt_run_to_stdy'] = \
        masc_ic_cnfg['run_configs']['run_ic']['run_tlrc']['abslt_flrt']

    param['config'] = masc_ic_cnfg['models_configs'][param['nom']]
    if "lmtd_dpth_fll_intl_lvl" in param['config'].keys():
        param['lmtd_dpth_fll_intl_lvl'] = \
        param['config']['lmtd_dpth_fll_intl_lvl']
    if "Froude_lmtd_bdry" in param['config'].keys():
        param['Froude_lmtd_bdry'] = param['config']['Froude_lmtd_bdry']
    if "tlrc_fll_intl_lvl" in param['config'].keys():
        if not "rltv" in param['config']['tlrc_fll_intl_lvl'].keys():
            raise NameError('model {}: relative tolerance for '
            'filling initial level missing.\nSTOP.'.format(param['nom']))
        if not "abslt" in param['config']['tlrc_fll_intl_lvl'].keys():
            raise NameError('model {}: absolute tolerance for '
            'filling initial level missing.\nSTOP.'.format(param['nom']))
        param['tlrc_rltv_fll_intl_lvl'] = \
        param['config']['tlrc_fll_intl_lvl']['rltv']
        param['tlrc_abslt_fll_intl_lvl'] = \
        param['config']['tlrc_fll_intl_lvl']['abslt']
    if "incrmt_fll_intl_lvl" in param['config'].keys():
        param['incrmt_fll_intl_lvl'] = param['config']['incrmt_fll_intl_lvl']

    param['iprint'] = masc_ic_cnfg['iprint']

    return param

def get_boundary_condition_cas_initial(chnl, tm_intl):
    """
    Get boundary condition at initial time stated in case file
    @param chnl(obj) instance of channel model study
    @param tm_intl(float) initial time of model
    @return bc1_intl, bc2_intl(dict) prescribed initial boundary condition
    """

    from telapy.tools.masc_boundary_condition import BoundaryConditionCas

    """
    only time-dependent boundary conditions
    """
    nb_bc_cas = chnl.masc.get_nb_cl()
    id_nom_bc = [chnl.masc.get_name_cl(i + 1) for i in range(nb_bc_cas)]
    id_bc = [id_nom[0] for id_nom in id_nom_bc]
    nom_bc = [id_nom[1] for id_nom in id_nom_bc]

    bc1_intl = {}
    bc2_intl = {}
    for i in range(nb_bc_cas):
        bc = BoundaryConditionCas(chnl, id_bc[i], nom_bc[i])
        if bc.type == 1 or bc.type == 2:
            bc1_intl[id_bc[i]] = bc.interpolate_boundary_condtion(tm_intl)
        elif bc.type == 3 or bc.type == 7:
            bc1_intl[id_bc[i]], bc2_intl[id_bc[i]] = \
            bc.interpolate_boundary_condtion(tm_intl)

    return bc1_intl, bc2_intl

def parse_initial_condition_prescribed(file_ic_prscrb, absc):
    """
    Parse prescribed initial condition file, 
    and checkes whether it includes complete contents
    @param file_ic_prscrb(str) path of prescribed initial condition file
    @param absc(list) relative abscissas of all cross-sections(same as input)
    @return ic_prscrb(array) water level and flowrate read from file
    """

    try:
        with open(file_ic_prscrb, 'r') as file:
            ic_prscrb = np.genfromtxt(file, dtype=float, delimiter=',', 
            skip_header=1, invalid_raise=True)
    except(ValueError, IOError):
        raise ValueError('incorrect input initial condition file.\nSTOP.')

    """
    raise error if there is NaN in abscissa
    """
    if not np.array_equal(
    np.asarray(absc, dtype=float), ic_prscrb[:, 0], equal_nan=False):
        raise ValueError('incorrect cross-sections '
        'in input initial condition file.\nSTOP.')
    """
    NaN is also interpreted as an element
    """
    if ic_prscrb[:, 1].size != len(absc):
        raise ValueError('incorrect number of level '
        'in input initial condition file.\nSTOP.')
    if ic_prscrb[:, 2].size != len(absc):
        raise ValueError('incorrect number of flowrate '
        'in input initial condition file.\nSTOP.')

    return ic_prscrb[:, 1], ic_prscrb[:, 2]

def fill_initial_level(chnl, param, bdry_tp, IDT, bc1_intl, xsID_fll=None):
    """
    Fill initial water level at cross-section when there is no input, 
    aiming to avoid super-critical flow at boundary, and no dry zone exists. 
    Regardless of input level.
    
    Initial condition could not be dependent on boundary condition. 
    If the boundary condition at initial time 
    is experienced over a sufficiently long period, 
    a run-to-steady should be executed using unsteady kernel.
    
    For a single reach with two prescribed level boundaries, 
    tests show that levels computed by steady kernel 
    could cause super-critical flow at boundary. 
    Because in steady kernel, 
    computed levels are dependent only on 
    prescribed level at start or end side, not both sides, 
    hence computed levels could be greatly different from 
    the prescribed level at the side 
    which is not assciated with steady-computation, 
    and so results in unlimited flowrate inducing super-critical flow.
    
    Execute steady kernel to get level 
    in reachs with prescribed level boundary, 
    defining limited flowrate by prescribed level and limited Froude number, 
    defining flow direction by bottom slope, 
    and setting boundary as downstream end when zero bottom slope.
    
    For level in other reaches: 
    if level of only one side has been filled, 
    equality of conveyance is set at all cross-sections in a reach; 
    if level of both sides have been filled, 
    use linear interpolation between two sides 
    in order to reduce distinguish of conveyance at two sides. 
    When dealing with complex network: 
    flow direction associated with real case cannot be defined, 
    since real flow could be opposite to direction of bottom slope. 
    Also, flowrate in other reaches is unknown, 
    since setting distribution of flowrate at junctions 
    based on magnitude of conveyance could not be the real case, 
    and also getting flowrate by limited Froude number 
    is appropriate for boundary but not for connection cross-section. 
    Moreover, there is numerical risk using steady kernel. 
    Hence, steady kernel is not used in other reaches.
    
    @param chnl(obj) instance of channel model study
    @param param(dict) loaded definitions of caller
    @param bdry_tp(list) model boundary condition type
    @param IDT(list) index of input model cross-section
    @param bc1_intl(dict) prescribed initial boundary condition
    @param xsID_fll(1d array) index of cross-section to be filled(from 0)
    @return z_fll(1d array) filled water level
    """

    def run_steady(z_fll, rchID_run, q_vct, z_prscrb, rvs_stt_ed):
        """
        execute steady-run to fill initial level of a reach
        @param z_fll(1d array) initial level to be filled
        @param rchID_run(int) index of reach to run
        @param q_vct(float) vetor flowrate as upstream boundary condition
        @param z_prscrb(float) prescribed level as downstream boundary condition
        @param rvs_stt_ed(bool) whether or not to reverse start and end
        @return z_fll(1d array) filled water level
        """

        from data_manip.formats.mascaretgeo_file import MascaretGeoFile
        from data_manip.formats.mascaret_file import Reach
        from execution.mascaret_cas import MascaretCas

        pth_prvs = os.getcwd()

        """
        working directory and model files
        """
        dir_rch = os.path.join('.', 
        '{}_reach{}'.format(param['nom'], (rchID_run + 1)))
        os.mkdir(dir_rch)
        dir_rch_files = os.path.join(dir_rch, 'files')
        os.mkdir(dir_rch_files)

        file_geo = os.path.join(dir_rch_files, 
        '{}_reach{}.geo'.format(param['nom'], (rchID_run + 1)))
        geo = MascaretGeoFile(
        os.path.join(chnl.paths['data'], param['config']['files']['geo']))

        xsID_stt_rch_run_geo_1 = 1
        for rchID in range(rchID_run):
            xsID_stt_rch_run_geo_1 += geo.reaches[(rchID + 1)].nsections
        xsID_ed_rch_run_geo_1 = \
        xsID_stt_rch_run_geo_1 + geo.reaches[(rchID_run + 1)].nsections - 1
        abscss_ed_rch_run_geo = geo.reaches[(rchID_run + 1)][-1].pk

        for rchID in list(geo.reaches.keys()):
            if rchID != (rchID_run + 1): geo.reaches.pop(rchID)
        if rvs_stt_ed:
            rch_std = Reach(1, 'reach_{}_steady'.format((rchID_run + 1)))
            j = geo.reaches[(rchID_run + 1)].nsections - 1
            for i in range(geo.reaches[(rchID_run + 1)].nsections):
                rch_std.add_section(geo.reaches[(rchID_run + 1)][j])
                rch_std[i].pk = \
                abscss_ed_rch_run_geo - geo.reaches[(rchID_run + 1)][j].pk
                rch_std[i].name = "p{}".format((i + 1))
                j -= 1
            geo.reaches[(rchID_run + 1)] = rch_std
            del rch_std
        geo.save(file_geo)

        dir_rch_files_bc = os.path.join(dir_rch_files, 'bc')
        os.mkdir(dir_rch_files_bc)
        file_bc_up = os.path.join(dir_rch_files_bc, 'up.loi')
        file_bc_down = os.path.join(dir_rch_files_bc, 'down.loi')
        with open(file_bc_up, 'w') as f:
            np.savetxt(f, np.array([[0, 0], [32e+06, 0]]), 
            fmt=('%14.3f'), delimiter=' ', header='S', comments='')
        with open(file_bc_down, 'w') as f:
            np.savetxt(f, np.array([[0, 0], [32e+06, 0]]), 
            fmt=('%14.3f'), delimiter=' ', header='S', comments='')

        file_cas = os.path.join(dir_rch_files, 
        '{}_reach{}.xcas'.format(param['nom'], (rchID_run + 1)))
        copy2(
        os.path.join(chnl.paths['data'], param['config']['files']['xcas']), 
        file_cas)
        cas = MascaretCas(file_cas, check_files=False, convert_xcas=True)

        cas.case.set('CALCULATION KERNEL', 1)
        cas.case.set('STORAGE AREAS', False)
        cas.case.set('COURLIS OPTION', False)
        cas.case.set('FLOOD WAVE CALCULATION', False)

        cas.case.set('FICHIER DE GEOMETRIE', os.path.basename(file_geo))
        cas.case.set('NOMBRE DE BRANCHES', 1)
        cas.case.set('BRANCHE NUMERO', [1])
        cas.case.set('ABSCISSE DEBUT', [geo.reaches[(rchID_run + 1)][0].pk])
        cas.case.set('ABSCISSE FIN', [geo.reaches[(rchID_run + 1)][-1].pk])
        cas.case.set("NUM DE L'EXTREMITE DE DEBUT", [1])
        cas.case.set("NUM DE L'EXTREMITE DE FIN", [2])
        cas.case.set('NODES NUMBER', 0)
        cas.case.set('FREE BOUNDARY NUMBER', 2)
        cas.case.set('EXTREMITE NUMERO', [1, 2])
        cas.case.set('NOM EXTREMITE', ['up', 'down'])
        cas.case.set('TYPE DE CONDITION', [1, 2])
        cas.case.set('NUMERO DE LA LOI', [1, 2])
        cas.case.set('NOMBRE DE CONFLUENTS', 0)

        xsID_stt_zn_vtcl_dscrtz_1 = cas.case.get('NUMERO DU PREMIER PROFIL')
        if type(xsID_stt_zn_vtcl_dscrtz_1) != list:
            xsID_stt_zn_vtcl_dscrtz_1 = [xsID_stt_zn_vtcl_dscrtz_1]
        xsID_ed_zn_vtcl_dscrtz_1 = cas.case.get('NUMERO DU DERNIER PROFIL')
        if type(xsID_ed_zn_vtcl_dscrtz_1) != list:
            xsID_ed_zn_vtcl_dscrtz_1 = [xsID_ed_zn_vtcl_dscrtz_1]
        znID_stt_vtcl_dscrtz = np.intersect1d(
        np.where(
        np.asarray(xsID_stt_zn_vtcl_dscrtz_1) <= xsID_stt_rch_run_geo_1), 
        np.where(np.asarray(xsID_ed_zn_vtcl_dscrtz_1) >= xsID_stt_rch_run_geo_1)
        )[0]
        znID_ed_vtcl_dscrtz = np.intersect1d(
        np.where(
        np.asarray(xsID_stt_zn_vtcl_dscrtz_1) <= xsID_ed_rch_run_geo_1), 
        np.where(np.asarray(xsID_ed_zn_vtcl_dscrtz_1) >= xsID_ed_rch_run_geo_1)
        )[0]
        xsID_stt_zn_vtcl_dscrtz_1 = list(map(int, list(
        np.asarray([xsID_stt_rch_run_geo_1] + xsID_stt_zn_vtcl_dscrtz_1[
        (znID_stt_vtcl_dscrtz + 1):(znID_ed_vtcl_dscrtz + 1)]) - 
        xsID_stt_rch_run_geo_1 + 1)))
        xsID_ed_zn_vtcl_dscrtz_1 = list(map(int, list(
        np.asarray(xsID_ed_zn_vtcl_dscrtz_1[
        znID_stt_vtcl_dscrtz:znID_ed_vtcl_dscrtz] + [xsID_ed_rch_run_geo_1]) - 
        xsID_stt_rch_run_geo_1 + 1)))
        if rvs_stt_ed:
            cas.case.set('NUMERO DU PREMIER PROFIL', list(map(int, 
            list(- np.flip(np.asarray(xsID_ed_zn_vtcl_dscrtz_1)) + 
            geo.reaches[(rchID_run + 1)].nsections + 1))))
            cas.case.set('NUMERO DU DERNIER PROFIL', list(map(int, 
            list(- np.flip(np.asarray(xsID_stt_zn_vtcl_dscrtz_1)) + 
            geo.reaches[(rchID_run + 1)].nsections + 1))))
        else:
            cas.case.set('NUMERO DU PREMIER PROFIL', xsID_stt_zn_vtcl_dscrtz_1)
            cas.case.set('NUMERO DU DERNIER PROFIL', xsID_ed_zn_vtcl_dscrtz_1)
        sz_vtcl_dscrtz = cas.case.get('VALEUR DU PAS')
        if type(sz_vtcl_dscrtz) != list:
            sz_vtcl_dscrtz = [sz_vtcl_dscrtz]
        if rvs_stt_ed:
            cas.case.set('VALEUR DU PAS', 
            list(np.flip(np.asarray(
            sz_vtcl_dscrtz[znID_stt_vtcl_dscrtz:(znID_ed_vtcl_dscrtz + 1)]))))
        else:
            cas.case.set('VALEUR DU PAS', 
            sz_vtcl_dscrtz[znID_stt_vtcl_dscrtz:(znID_ed_vtcl_dscrtz + 1)])
        cas.case.set('NOMBRE DE ZONES DE PLANIMETRAGE', 
        int(znID_ed_vtcl_dscrtz - znID_stt_vtcl_dscrtz + 1))

        """
        longtitudinal discretization approach by input file 
        is not supported by the kernel
        """
        if cas.case.get('MODE DE SAISIE DU MAILLAGE') != 2:
            raise ValueError('model {}: longtitudinal discretization approach '
            'by input file not supported.\nSTOP.'.format(param['nom']))
        tp_lgtdl_dscrtz = cas.case.get('METHODE DE CALCUL DU MAILLAGE')
        if tp_lgtdl_dscrtz == 2:
            rchID_zn_lgtdl_dscrtz = cas.case.get('NUMERO DE BRANCHE DE ZONE')
            if type(rchID_zn_lgtdl_dscrtz) != list:
                rchID_zn_lgtdl_dscrtz = [rchID_zn_lgtdl_dscrtz]
            znID_run = np.sort(
            np.where(np.asarray(rchID_zn_lgtdl_dscrtz) == (rchID_run + 1))[0])
            cas.case.set('NUMERO DE BRANCHE DE ZONE', 
            list(map(int, list(np.ones(znID_run.size)))))
            abscss_stt_zn_lgtdl_dscrtz = \
            cas.case.get('ABSCISSE DE DEBUT DE ZONE')
            if type(abscss_stt_zn_lgtdl_dscrtz) != list:
                abscss_stt_zn_lgtdl_dscrtz = [abscss_stt_zn_lgtdl_dscrtz]
            abscss_ed_zn_lgtdl_dscrtz = cas.case.get('ABSCISSE DE FIN DE ZONE')
            if type(abscss_ed_zn_lgtdl_dscrtz) != list:
                abscss_ed_zn_lgtdl_dscrtz = [abscss_ed_zn_lgtdl_dscrtz]
            nb_nd_zn_lgtdl_dscrtz = \
            cas.case.get('NOMBRE DE SECTIONS DE LA ZONE')
            if type(nb_nd_zn_lgtdl_dscrtz) != list:
                nb_nd_zn_lgtdl_dscrtz = [nb_nd_zn_lgtdl_dscrtz]
            if rvs_stt_ed:
                abscss_stt_zn_lgtdl_dscrtz_dcml = []
                for abscss in list(map(str, (
                abscss_ed_rch_run_geo - 
                np.asarray(abscss_ed_zn_lgtdl_dscrtz)[znID_run]))):
                    abscss_stt_zn_lgtdl_dscrtz_dcml.append(
                    float(Decimal(abscss).quantize(Decimal('0.000'))))
                cas.case.set('ABSCISSE DE DEBUT DE ZONE', 
                list(np.flip(np.asarray(abscss_stt_zn_lgtdl_dscrtz_dcml))))
                abscss_ed_zn_lgtdl_dscrtz_dcml = []
                for abscss in list(map(str, (
                abscss_ed_rch_run_geo - 
                np.asarray(abscss_stt_zn_lgtdl_dscrtz)[znID_run]))):
                    abscss_ed_zn_lgtdl_dscrtz_dcml.append(
                    float(Decimal(abscss).quantize(Decimal('0.000'))))
                cas.case.set('ABSCISSE DE FIN DE ZONE', 
                list(np.flip(np.asarray(abscss_ed_zn_lgtdl_dscrtz_dcml))))
                cas.case.set('NOMBRE DE SECTIONS DE LA ZONE', list(map(int, 
                list(np.flip(np.asarray(nb_nd_zn_lgtdl_dscrtz)[znID_run])))))
            else:
                cas.case.set('ABSCISSE DE DEBUT DE ZONE', 
                list(np.asarray(abscss_stt_zn_lgtdl_dscrtz)[znID_run]))
                cas.case.set('ABSCISSE DE FIN DE ZONE', 
                list(np.asarray(abscss_ed_zn_lgtdl_dscrtz)[znID_run]))
                cas.case.set('NOMBRE DE SECTIONS DE LA ZONE', list(map(int, 
                list(np.asarray(nb_nd_zn_lgtdl_dscrtz)[znID_run]))))
            cas.case.set('NOMBRE DE ZONES DE DISCRETISATION', znID_run.size)
        elif tp_lgtdl_dscrtz == 3:
            rchID_sct_lgtdl_dscrtz = \
            cas.case.get('BRANCHES DES SECTIONS DE CALCUL')
            if type(rchID_sct_lgtdl_dscrtz) != list:
                rchID_sct_lgtdl_dscrtz = [rchID_sct_lgtdl_dscrtz]
            sctID_run = np.sort(
            np.where(np.asarray(rchID_sct_lgtdl_dscrtz) == (rchID_run + 1))[0])
            cas.case.set('BRANCHES DES SECTIONS DE CALCUL', 
            list(map(int, list(np.ones(sctID_run.size)))))
            abscss_nd_lgtdl_dscrtz = \
            cas.case.get('ABSCISSES DES SECTIONS DE CALCUL')
            if type(abscss_nd_lgtdl_dscrtz) != list:
                abscss_nd_lgtdl_dscrtz = [abscss_nd_lgtdl_dscrtz]
            if rvs_stt_ed:
                abscss_nd_lgtdl_dscrtz_dcml = []
                for abscss in list(map(str, (abscss_ed_rch_run_geo - 
                np.asarray(abscss_nd_lgtdl_dscrtz)[sctID_run]))):
                    abscss_nd_lgtdl_dscrtz_dcml.append(
                    float(Decimal(abscss).quantize(Decimal('0.000'))))
                cas.case.set('ABSCISSES DES SECTIONS DE CALCUL', 
                list(np.flip(np.asarray(abscss_nd_lgtdl_dscrtz_dcml))))
            else:
                cas.case.set('ABSCISSES DES SECTIONS DE CALCUL', 
                list(np.asarray(abscss_nd_lgtdl_dscrtz)[sctID_run]))
            cas.case.set('NOMBRE DE SECTIONS DE CALCUL', sctID_run.size)
        elif tp_lgtdl_dscrtz == 5:
            xsID_stt_rg_lgtdl_dscrtz_1 = \
            cas.case.get('NUMERO DU PREMIER PROFIL DE LA SERIE')
            if type(xsID_stt_rg_lgtdl_dscrtz_1) != list:
                xsID_stt_rg_lgtdl_dscrtz_1 = [xsID_stt_rg_lgtdl_dscrtz_1]
            xsID_ed_rg_lgtdl_dscrtz_1 = \
            cas.case.get('NUMERO DU DERNIER PROFIL DE LA SERIE')
            if type(xsID_ed_rg_lgtdl_dscrtz_1) != list:
                xsID_ed_rg_lgtdl_dscrtz_1 = [xsID_ed_rg_lgtdl_dscrtz_1]
            rgID_stt_lgtdl_dscrtz = \
            np.where(
            np.asarray(xsID_stt_rg_lgtdl_dscrtz_1) == xsID_stt_rch_run_geo_1
            )[0][0]
            rgID_ed_lgtdl_dscrtz = \
            np.where(
            np.asarray(xsID_ed_rg_lgtdl_dscrtz_1) == xsID_ed_rch_run_geo_1
            )[0][0]
            xsID_stt_rg_lgtdl_dscrtz_1 = list(map(int, list(
            np.asarray(xsID_stt_rg_lgtdl_dscrtz_1[
            rgID_stt_lgtdl_dscrtz:(rgID_ed_lgtdl_dscrtz + 1)]) - 
            xsID_stt_rch_run_geo_1 + 1)))
            xsID_ed_rg_lgtdl_dscrtz_1 = list(map(int, list(
            np.asarray(xsID_ed_rg_lgtdl_dscrtz_1[
            rgID_stt_lgtdl_dscrtz:(rgID_ed_lgtdl_dscrtz + 1)]) - 
            xsID_stt_rch_run_geo_1 + 1)))
            sz_rg_lgtdl_dscrtz = cas.case.get("PAS D'ESPACE DE LA SERIE")
            if type(sz_rg_lgtdl_dscrtz) != list:
                sz_rg_lgtdl_dscrtz = [sz_rg_lgtdl_dscrtz]
            if rvs_stt_ed:
                cas.case.set('NUMERO DU PREMIER PROFIL DE LA SERIE', 
                list(map(int, 
                list(- np.flip(np.asarray(xsID_ed_rg_lgtdl_dscrtz_1)) + 
                geo.reaches[(rchID_run + 1)].nsections + 1))))
                cas.case.set('NUMERO DU DERNIER PROFIL DE LA SERIE', 
                list(map(int, 
                list(- np.flip(np.asarray(xsID_stt_rg_lgtdl_dscrtz_1)) + 
                geo.reaches[(rchID_run + 1)].nsections + 1))))
                cas.case.set("PAS D'ESPACE DE LA SERIE", 
                list(np.flip(np.asarray(sz_rg_lgtdl_dscrtz[
                rgID_stt_lgtdl_dscrtz:(rgID_ed_lgtdl_dscrtz + 1)]))))
            else:
                cas.case.set('NUMERO DU PREMIER PROFIL DE LA SERIE', 
                xsID_stt_rg_lgtdl_dscrtz_1)
                cas.case.set('NUMERO DU DERNIER PROFIL DE LA SERIE', 
                xsID_ed_rg_lgtdl_dscrtz_1)
                cas.case.set("PAS D'ESPACE DE LA SERIE", 
                sz_rg_lgtdl_dscrtz[
                rgID_stt_lgtdl_dscrtz:(rgID_ed_lgtdl_dscrtz + 1)])
            cas.case.set('NOMBRE DE PLAGES DE DISCRETISATION', 
            int(rgID_ed_lgtdl_dscrtz - rgID_stt_lgtdl_dscrtz + 1))

        """
        Initial level should be dependent on 
        prescribed level at end and limited Fround number, 
        not weir, regulation, energy dissipation, 
        lateral flow, lateral weir, or storage, 
        because these objects could cause high level 
        which results in supercritical or low level which results in dry.
        """
        cas.case.remove('NUM BRANCHE DU BARRAGE PRINCIPAL')
        cas.case.remove('ABSCISSE DU BARRAGE PRINCIPAL')
        cas.case.remove('TYPE DE RUPTURE DU BARRAGE PRINCIPAL')
        cas.case.remove('COTE DE CRETE DU BARRAGE PRINCIPAL')
        cas.case.set('WEIRS NUMBER', 0)
        cas.case.set('NOMBRE DE PERTES DE CHARGE SINGULIERES', 0)

        cas.case.set('NOMBRE DE CASIERS', 0)

        cas.case.set('LATERAL INFLOW DISCHARGES NUMBER', 0)
        cas.case.set('NOMBRE DE DEVERSOIRS', 0)
        cas.case.set("NOMBRE D'APPORTS DE PLUIE", 0)

        cas.case.set('TRACERS PRESENCE', False)
        cas.case.set('CALCULATION WITH SAND', False)

        cas.case.set('NUMBER OF HYDRAULIC LAWS', 2)
        cas.case.set('LAW NAME', ['up', 'down'])
        cas.case.set('LOI TYPE', [1, 2])
        cas.case.set("LOIS MODE D'ENTREE", [1, 1])
        cas.case.set('LOIS FICHIER', ['up.loi', 'down.loi'])

        cas.case.set('VARIABLE TIME STEP WITH COURANT NUMBER', False)

        cas.case.set('STORAGE OPTION', 1)

        """
        For case of start and end is reversed, 
        friction coefficient is set after steady model imported.
        """
        cas.case.set('OPTION AUTO CALIBRATION', False)
        if rvs_stt_ed:
            cas.case.set('NUMBER OF FRICTION ZONES', 1)
            cas.case.set('REACH NUMBER FOR THE FRICTION ZONE', [1])
            cas.case.set('FRICTION ZONE UPSTREAM ABSCISSA', 
            [geo.reaches[(rchID_run + 1)][0].pk])
            cas.case.set('FRICTION ZONE DOWNSTREAM ABSCISSA', 
            [geo.reaches[(rchID_run + 1)][-1].pk])
            cas.case.set('MAIN CHANNEL COEFFICIENT', [0.03])
            cas.case.set('FLOODPLAIN COEFFICIENT', [0.03])
        else:
            rchID_zn_frct_cffct = \
            cas.case.get('REACH NUMBER FOR THE FRICTION ZONE')
            if type(rchID_zn_frct_cffct) != list:
                rchID_zn_frct_cffct = [rchID_zn_frct_cffct]
            abscss_stt_zn_frct_cffct = \
            cas.case.get('FRICTION ZONE UPSTREAM ABSCISSA')
            if type(abscss_stt_zn_frct_cffct) != list:
                abscss_stt_zn_frct_cffct = [abscss_stt_zn_frct_cffct]
            abscss_ed_zn_frct_cffct = \
            cas.case.get('FRICTION ZONE DOWNSTREAM ABSCISSA')
            if type(abscss_ed_zn_frct_cffct) != list:
                abscss_ed_zn_frct_cffct = [abscss_ed_zn_frct_cffct]
            cffct_mc_zn_frct_cffct = cas.case.get('MAIN CHANNEL COEFFICIENT')
            if type(cffct_mc_zn_frct_cffct) != list:
                cffct_mc_zn_frct_cffct = [cffct_mc_zn_frct_cffct]
            cffct_fp_zn_frct_cffct = cas.case.get('FLOODPLAIN COEFFICIENT')
            if type(cffct_fp_zn_frct_cffct) != list:
                cffct_fp_zn_frct_cffct = [cffct_fp_zn_frct_cffct]
            znID_run = np.sort(
            np.where(np.asarray(rchID_zn_frct_cffct) == (rchID_run + 1))[0])
            cas.case.set('NUMBER OF FRICTION ZONES', znID_run.size)
            cas.case.set('REACH NUMBER FOR THE FRICTION ZONE', 
            list(map(int, list(np.ones(znID_run.size)))))
            cas.case.set('FRICTION ZONE UPSTREAM ABSCISSA', 
            list(np.asarray(abscss_stt_zn_frct_cffct)[znID_run]))
            cas.case.set('FRICTION ZONE DOWNSTREAM ABSCISSA', 
            list(np.asarray(abscss_ed_zn_frct_cffct)[znID_run]))
            cas.case.set('MAIN CHANNEL COEFFICIENT', 
            list(np.asarray(cffct_mc_zn_frct_cffct)[znID_run]))
            cas.case.set('FLOODPLAIN COEFFICIENT', 
            list(np.asarray(cffct_fp_zn_frct_cffct)[znID_run]))

        cas.write_xcas_file()
        del cas

        cnfg_rch = {}
        cnfg_rch['files'] = {}
        cnfg_rch['files']['xcas'] = file_cas
        cnfg_rch['files']['geo'] = file_geo
        cnfg_rch['files']['res'] = \
        "result_{}_reach{}.csv".format(param['nom'], (rchID_run + 1))
        cnfg_rch['files']['listing'] = \
        "listing_{}_reach{}.lis".format(param['nom'], (rchID_run + 1))
        cnfg_rch['files']['damocle'] = \
        "pre_{}_reach{}.damoc".format(param['nom'], (rchID_run + 1))
        cnfg_rch['files']['loi'] = [file_bc_up, file_bc_down]

        """
        import steady model, revise coefficient, and compute
        """
        chnl_rch = MascaretStudy(cnfg_rch, log_lvl='CRITICAL', iprint=0, 
        working_directory='{}_reach{}'.format(param['nom'], (rchID_run + 1)))
        IDT_mdl_rch = np.array([chnl_rch.masc.get('Model.IDT', i) 
        for i in range(chnl_rch.model_size)], dtype=np.int_)
        IDT_mdl_rch[-1] = xsID_ed_rch_run_geo_1 - xsID_stt_rch_run_geo_1 + 1
        ndID_all_mdl_rch = [i for i in range(chnl_rch.model_size)]

        """
        If start and end is reversed, 
        original friction coefficient of computational nodes could be changed. 
        By the kernel, the coefficient of upstream is set for 
        the computational nodes between gap of zones, and after reverse, 
        the coefficient of downstream will be set for 
        those computational nodes. 
        Then the coefficient of original model and 
        this steady model will be different, which induces incoherence. 
        Here the coefficient of original model is set for 
        each computational node. 
        The relative location of all computational nodes is not changed by 
        any type of longitudinal discretization, 
        although abscissa could be changed.
        """
        if rvs_stt_ed:
            ndID_stt = chnl.masc.get('Model.Connect.FirstNdNum', rchID_run) - 1
            ndID_ed = chnl.masc.get('Model.Connect.LastNdNum', rchID_run) - 1
            if chnl_rch.model_size != (ndID_ed - ndID_stt + 1):
                raise Exception('model {}: '
                'number of computational nodes of reversed steady-run reach '
                'different from that of original model.\nSTOP.'
                .format(param['nom']))
            cffct_mc = chnl.get_friction_minor()
            cffct_fp = chnl.get_friction_major()
            for i in range(chnl_rch.model_size):
                chnl_rch.set_friction_minor({
                'value': cffct_mc[ndID_stt:(ndID_ed + 1)]
                [(chnl_rch.model_size - i - 1)], 
                'index': i
                })
                chnl_rch.set_friction_major({
                'value': cffct_fp[ndID_stt:(ndID_ed + 1)]
                [(chnl_rch.model_size - i - 1)], 
                'index': i
                })

        """
        Computation of steady-run does not require initial conditions, 
        but in the kernel, 
        all variables regarding model state need to be initialized, 
        otherwise these variables are not allocated or false randoms.
        """
        os.chdir(chnl_rch.paths['output'])
        chnl_rch.initialize_model(
        z_init=np.array([min([
        chnl_rch.masc.get('Model.LevRightBk', i), 
        chnl_rch.masc.get('Model.LevLeftBk', i)]) 
        for i in range(chnl_rch.model_size)]), 
        q_init=np.zeros(chnl_rch.model_size, dtype=np.double))
        os.chdir(chnl_rch.paths['study'])

        chnl_rch.masc.compute_bc(0.0, 1.0, 1.0, np.array([0.0, 1.0]), 2, 2, 
        np.array([[q_vct, z_prscrb], [q_vct, z_prscrb]]), 
        np.array([[0.0, 0.0], [0.0, 0.0]]))

        for xsID_mdl_rch in np.where(
        np.isnan(z_fll[(xsID_stt_rch_run_geo_1 - 1):xsID_ed_rch_run_geo_1]))[0]:
            xsID_mdl = xsID_mdl_rch + xsID_stt_rch_run_geo_1 - 1
            if rvs_stt_ed:
                xsID_mdl_rch = \
                geo.reaches[(rchID_run + 1)].nsections - xsID_mdl_rch - 1
            z_fll[xsID_mdl] = \
            chnl_rch.masc.get('State.Z', ndID_all_mdl_rch[
            np.searchsorted(IDT_mdl_rch, (xsID_mdl_rch + 1), side='left')])

        del geo
        del chnl_rch

        os.chdir(pth_prvs)

        return z_fll

    if not "Froude_lmtd_bdry" in param.keys():
        raise NameError('model {}: subcritical kernel and '
        'boundary condition type of prescribed level but '
        'boundary maximum Frounde number missing.\nSTOP.'.format(param['nom']))
    if param['Froude_lmtd_bdry'] < 0.0:
        raise ValueError(
        'model {}: negative limited Froude number.\nSTOP.'.format(param['nom']))
    if not "incrmt_fll_intl_lvl" in param.keys():
        raise NameError('model {}: subcritical kernel and '
        'boundary condition type of prescribed level but '
        'iteration level increment missing.\nSTOP.'.format(param['nom']))
    if param['incrmt_fll_intl_lvl'] < 0.0:
        raise ValueError(
        'model {}: negative interation increment.\nSTOP.'.format(param['nom']))
    if not "tlrc_rltv_fll_intl_lvl" in param.keys() or \
    not "tlrc_abslt_fll_intl_lvl" in param.keys():
        raise NameError('model {}: subcritical kernel and '
        'boundary condition type of prescribed level but '
        'iteration tolerance missing.\nSTOP.'.format(param['nom']))
    if param['tlrc_rltv_fll_intl_lvl'] < 0.0:
        raise ValueError(
        'model {}: negative relative tolerance.\nSTOP.'.format(param['nom']))
    if param['tlrc_abslt_fll_intl_lvl'] < 0.0:
        raise ValueError(
        'model {}: negative absolute tolerance.\nSTOP.'.format(param['nom']))

    """
    General model property.
    """
    nb_xs = len(chnl.cross_sections['bottom']['abscissa'])
    IDT = np.array(IDT)
    nb_cnnct = chnl.masc.get_var_size('Model.Connect.NumReachJunction')[0]
    z_bttm = [chnl.masc.get('Model.CrossSection.Zbot', i) for i in range(nb_xs)]

    z_fll = np.zeros(nb_xs, dtype=np.double)
    z_fll[:] = np.nan
    xsID_all = np.array([i for i in range(nb_xs)])
    if xsID_fll is None:
        xsID_fll = np.array([i for i in range(nb_xs)])
    ndID_all = np.array([i for i in range(chnl.model_size)])

    """
    single reach model, and prescribed level at both start and end boundaries
    """
    if nb_cnnct == 0 and np.unique(np.asarray(bdry_tp)).size == 1:
        if not "lmtd_dpth_fll_intl_lvl" in param.keys():
            raise NameError('model {}: subcritical kernel, single reach and '
            'boundary condition type of prescribed level at both boundaries, '
            'but limited water depth missing.\nSTOP.'.format(param['nom']))

        for bdryID in range(2):
            z_fll[(
            IDT[(chnl.masc.get('Model.Connect.NodeNumFreeOutflow', bdryID) - 1)]
             - 1)] = bc1_intl[chnl.masc.get('Model.Boundary.GraphNum', bdryID)]
            if np.all(np.isin(xsID_fll, 
            np.delete(xsID_all, np.where(np.isnan(z_fll))[0]), 
            assume_unique=True)):
                return z_fll

        z_fll[1:(nb_xs - 1)] = np.interp(
        np.asarray(chnl.cross_sections['bottom']['abscissa'][1:(nb_xs - 1)]), 
        np.asarray([chnl.cross_sections['bottom']['abscissa'][0], 
        chnl.cross_sections['bottom']['abscissa'][(nb_xs - 1)]]), 
        z_fll[[0, (nb_xs - 1)]])
        for xsID_slc in np.where(
        (z_fll[1:(nb_xs - 1)] - np.asarray(z_bttm[1:(nb_xs - 1)])) < 
        param['lmtd_dpth_fll_intl_lvl'])[0]:
            z_fll[(xsID_slc + 1)] = \
            z_bttm[(xsID_slc + 1)] + param['lmtd_dpth_fll_intl_lvl']

        if np.all(np.isin(xsID_fll, 
        np.delete(xsID_all, np.where(np.isnan(z_fll))[0]), 
        assume_unique=True)):
            return z_fll

        raise Exception(
        'model {}: filling initial level for single reach failed.\nSTOP.'
        .format(param['nom']))

    """
    Reach with prescribed level boundary.
    
    In steady-flow kernel, end cross-section must be prescribed-level.
    """
    for bdryID in np.where(np.asarray(bdry_tp) == 2)[0]:
        ndID_bdry = \
        chnl.masc.get('Model.Connect.NodeNumFreeOutflow', bdryID) - 1
        xsID_bdry = IDT[ndID_bdry] - 1
        rchID = chnl.masc.get('Model.Connect.ReachNumFreeOutflow', bdryID) - 1
        ndID_stt = chnl.masc.get('Model.Connect.FirstNdNum', rchID) - 1
        ndID_ed = chnl.masc.get('Model.Connect.LastNdNum', rchID) - 1
        bcID = chnl.masc.get('Model.Boundary.GraphNum', bdryID) - 1

        z_fll[xsID_bdry] = bc1_intl[(bcID + 1)]
        if np.all(np.isin(xsID_fll, 
        np.delete(xsID_all, np.where(np.isnan(z_fll))[0]), 
        assume_unique=True)):
            return z_fll

        dpth_prscrb = bc1_intl[(bcID + 1)] - z_bttm[xsID_bdry]
        if (dpth_prscrb > 
        (chnl.masc.get('Model.CrossSection.Step', xsID_bdry) * 
        (chnl.masc.get('Model.CrossSection.NumStep', xsID_bdry) - 1))) \
        or (dpth_prscrb <= 0.):
            raise ValueError('model {}: initial level of '
            'input boundary condition No.{} exceeds vertical discretisation '
            'range.\nSTOP.'.format(param['nom'], (bcID + 1)))
        _, _, _, _, _, area_mc, area_fp, _, _, _, beta = \
        chnl.masc.get_hydraulic_variable(ndID_bdry, bc1_intl[(bcID + 1)])

        """
        Exactly same algorithm as FORTRAN subroutine "FROUDE_S".
        Compound channel coefficient is always 1.0 for "bottom/bank" method.
        """
        if chnl.masc.get('Model.CSectionLayout') == 2:
            vlct_mn_lmtd = \
            (((9.81 * dpth_prscrb) ** 0.5) * param['Froude_lmtd_bdry'])
        else:
            dnmnt = \
            ((beta / param['Froude_lmtd_bdry']) ** 2) - ((beta - 1.) * beta)
            if dnmnt <= 0.:
                raise ValueError('model {}: input boundary maximum Froude '
                'too large.\nSTOP.'.format(param['nom']))
            vlct_mn_lmtd = ((9.81 * dpth_prscrb) / dnmnt) ** 0.5

        if ndID_bdry == ndID_ed:
            if chnl.masc.get('Model.Zbot', ndID_bdry) <= \
            chnl.masc.get('Model.Zbot', ndID_stt):
                q_vct = vlct_mn_lmtd * (area_mc + area_fp)
            else:
                q_vct = - (vlct_mn_lmtd * (area_mc + area_fp))
            rvs_stt_ed = False
        elif ndID_bdry == ndID_stt:
            if chnl.masc.get('Model.Zbot', ndID_bdry) <= \
            chnl.masc.get('Model.Zbot', ndID_ed):
                q_vct = vlct_mn_lmtd * (area_mc + area_fp)
            else:
                q_vct = - (vlct_mn_lmtd * (area_mc + area_fp))
            rvs_stt_ed = True
        z_fll = run_steady(
        z_fll, rchID, q_vct, bc1_intl[(bcID + 1)], rvs_stt_ed)
        if np.all(np.isin(xsID_fll, 
        np.delete(xsID_all, np.where(np.isnan(z_fll))[0]), 
        assume_unique=True)):
            return z_fll

    if nb_cnnct == 0:
        raise Exception(
        'model {}: filling initial level for single reach failed.\nSTOP.'
        .format(param['nom']))

    """
    other reaches if there is connection
    """

    """
    several general array variables
    """
    nb_rch_cnnct = [chnl.masc.get('Model.Connect.NumReachJunction', 
    i) for i in range(nb_cnnct)]
    rchID_cnnct = np.zeros((nb_cnnct, max(nb_rch_cnnct)), dtype=np.int_)
    ndID_cnnct = np.zeros((nb_cnnct, max(nb_rch_cnnct)), dtype=np.int_)
    for i in range(nb_cnnct):
        for j in range(nb_rch_cnnct[i]):
            rchID_cnnct[i, j] = \
            chnl.masc.get('Model.Connect.ReachNum', i, j)
            ndID_cnnct[i, j] = chnl.masc.get('Model.Connect.NodeNum', i, j)
    rchID_cnnct -= 1
    ndID_cnnct -= 1
    cnnctID_intnlID_rch_nan = np.asarray(np.where(rchID_cnnct == -1))

    """
    loop to fill level
    """
    z_cnnct = np.zeros((nb_cnnct, max(nb_rch_cnnct)), dtype=np.double)
    z_cnnct[:, :] = np.nan
    for i in range(nb_cnnct):
        for j in range(nb_rch_cnnct[i]):
            z_cnnct[i, j] = z_fll[(IDT[ndID_cnnct[i, j]] - 1)]
    cnnctID_intnlID_flld = \
    np.asarray(np.where(np.logical_not(np.isnan(z_cnnct))))
    cnnctID_intnlID_to_fll = np.asarray(np.where(np.isnan(z_cnnct)))
    for cnnctID_intnlID in cnnctID_intnlID_rch_nan.T:
        cnnctID_intnlID_to_fll = np.delete(cnnctID_intnlID_to_fll, 
        np.intersect1d(
        np.where(cnnctID_intnlID_to_fll[0] == cnnctID_intnlID[0])[0], 
        np.where(cnnctID_intnlID_to_fll[1] == cnnctID_intnlID[1])[0]), 
        axis=1)
    cnnctID_to_fll = \
    np.intersect1d(cnnctID_intnlID_flld[0], cnnctID_intnlID_to_fll[0])

    while cnnctID_to_fll.size:
        intnlID_to_fll = {}
        rchID_spcl = []

        """
        first, cross-sections consisting of connection
        """
        for cnnctID in cnnctID_to_fll:
            intnlID_to_fll[cnnctID] = cnnctID_intnlID_to_fll[1]\
            [np.where(cnnctID_intnlID_to_fll[0] == cnnctID)]

            z_fll[(IDT[ndID_cnnct[cnnctID, intnlID_to_fll[cnnctID]]] - 1)] = \
            max(z_cnnct[cnnctID, 
            cnnctID_intnlID_flld[1]
            [np.where(cnnctID_intnlID_flld[0] == cnnctID)]])
            if np.all(np.isin(xsID_fll, 
            np.delete(xsID_all, np.where(np.isnan(z_fll))[0]), 
            assume_unique=True)):
                return z_fll

        """
        then, reach by reach, except for that with both sides filled
        """
        for cnnctID in cnnctID_to_fll:
            for intnlID in intnlID_to_fll[cnnctID]:
                ndID_stt = chnl.masc.get(
                'Model.Connect.FirstNdNum', rchID_cnnct[cnnctID, intnlID]) - 1
                ndID_ed = chnl.masc.get(
                'Model.Connect.LastNdNum', rchID_cnnct[cnnctID, intnlID]) - 1

                if np.all(
                np.logical_not(np.isnan(z_fll[
                [(IDT[ndID_stt] - 1), (IDT[ndID_ed] - 1)]]))):
                    rchID_spcl.append(rchID_cnnct[cnnctID, intnlID])
                else:
                    dpth_fll = z_fll[(IDT[ndID_cnnct[cnnctID, intnlID]] - 1)] \
                    - z_bttm[(IDT[ndID_cnnct[cnnctID, intnlID]] - 1)]
                    if (dpth_fll > 
                    (chnl.masc.get('Model.CrossSection.Step', 
                    (IDT[ndID_cnnct[cnnctID, intnlID]] - 1)) * 
                    (chnl.masc.get('Model.CrossSection.NumStep', 
                    (IDT[ndID_cnnct[cnnctID, intnlID]] - 1)) - 1))) \
                    or (dpth_fll <= 0.):
                        raise ValueError(
                        'model {}: filling level at connection No.{} '
                        'exceeds vertical discretisation range of '
                        'cross-section No.{} of reach No.{}. '
                        'Revise cross-section geometry of the connection.'
                        '\nSTOP.'.format(param['nom'], (cnnctID + 1), 
                        IDT[ndID_cnnct[cnnctID, intnlID]], 
                        (rchID_cnnct[cnnctID, intnlID] + 1)))
                    _, _, _, _, _, _, _, _, _, cnvc_rch, _ = \
                    chnl.masc.get_hydraulic_variable(
                    ndID_cnnct[cnnctID, intnlID], 
                    z_fll[(IDT[ndID_cnnct[cnnctID, intnlID]] - 1)])
                    if ndID_cnnct[cnnctID, intnlID] == ndID_stt:
                        xsID_to_fll = [i for i in range(
                        (IDT[ndID_stt]), IDT[ndID_ed])]
                    elif ndID_cnnct[cnnctID, intnlID] == ndID_ed:
                        xsID_to_fll = [i for i in range(
                        (IDT[ndID_stt] - 1), (IDT[ndID_ed] - 1))]

                    for xsID in xsID_to_fll:
                        ndID = \
                        ndID_all[np.searchsorted(IDT, (xsID + 1), side='left')]
                        z_lmtd = min(
                        (z_bttm[xsID] + 
                        (chnl.masc.get('Model.CrossSection.Step', xsID) * 
                        (chnl.masc.get('Model.CrossSection.NumStep', xsID) - 
                        1))), 
                        min([chnl.masc.get('Model.CrossSection.Zbank', xsID, i) 
                        for i in range(2)])
                        )

                        _, _, _, _, _, _, _, _, _, cnvc, _ = \
                        chnl.masc.get_hydraulic_variable(ndID, z_lmtd)
                        if cnvc < cnvc_rch or \
                        np.isclose(cnvc, cnvc_rch, 
                        rtol=param['tlrc_rltv_fll_intl_lvl'], 
                        atol=param['tlrc_abslt_fll_intl_lvl'], 
                        equal_nan=False):
                            z_fll[xsID] = z_lmtd
                        else:
                            z_prdct = z_lmtd - param['incrmt_fll_intl_lvl']
                            while True:
                                if z_prdct <= z_bttm[xsID]:
                                    raise Exception('model {}: '
                                    'filling initial level for '
                                    'cross-section No.{} failed.\nSTOP.'
                                    .format(param['nom'], (xsID + 1)))
                                _, _, _, _, _, _, _, _, _, cnvc, _ = \
                                chnl.masc.get_hydraulic_variable(
                                ndID, z_prdct)
                                if np.isclose(cnvc, cnvc_rch, 
                                rtol=param['tlrc_rltv_fll_intl_lvl'], 
                                atol=param['tlrc_abslt_fll_intl_lvl'], 
                                equal_nan=False):
                                    z_fll[xsID] = z_prdct
                                    break
                                else:
                                    z_prdct -= param['incrmt_fll_intl_lvl']

                        if np.all(np.isin(xsID_fll, 
                        np.delete(xsID_all, np.where(np.isnan(z_fll))[0]), 
                        assume_unique=True)):
                            return z_fll

        """
        reach with both sides filled
        """
        if rchID_spcl and not "lmtd_dpth_fll_intl_lvl" in param.keys():
            raise NameError('model {}: subcritical kernel, single reach and '
            'boundary condition type of prescribed level at both boundaries, '
            'but limited water depth missing.\nSTOP.'.format(param['nom']))

        for rchID in set(rchID_spcl):
            ndID_stt = chnl.masc.get('Model.Connect.FirstNdNum', rchID) - 1
            ndID_ed = chnl.masc.get('Model.Connect.LastNdNum', rchID) - 1

            z_fll[IDT[ndID_stt]:(IDT[ndID_ed] - 1)] = np.interp(
            np.asarray(chnl.cross_sections['bottom']['abscissa']
            [IDT[ndID_stt]:(IDT[ndID_ed] - 1)]), 
            np.asarray(
            [chnl.cross_sections['bottom']['abscissa'][(IDT[ndID_stt] - 1)], 
            chnl.cross_sections['bottom']['abscissa'][(IDT[ndID_ed] - 1)]]), 
            z_fll[[(IDT[ndID_stt] - 1), (IDT[ndID_ed] - 1)]])
            for xsID_slc in np.where(
            (z_fll[IDT[ndID_stt]:(IDT[ndID_ed] - 1)] - 
            np.asarray(z_bttm[IDT[ndID_stt]:(IDT[ndID_ed] - 1)])) < 
            param['lmtd_dpth_fll_intl_lvl'])[0]:
                z_fll[(xsID_slc + IDT[ndID_stt])] = \
                z_bttm[(xsID_slc + IDT[ndID_stt])] + \
                param['lmtd_dpth_fll_intl_lvl']
            if np.all(np.isin(xsID_fll, 
            np.delete(xsID_all, np.where(np.isnan(z_fll))[0]), 
            assume_unique=True)):
                return z_fll

        for i in range(nb_cnnct):
            for j in range(nb_rch_cnnct[i]):
                z_cnnct[i, j] = z_fll[(IDT[ndID_cnnct[i, j]] - 1)]
        cnnctID_intnlID_flld = \
        np.asarray(np.where(np.logical_not(np.isnan(z_cnnct))))
        cnnctID_intnlID_to_fll = np.asarray(np.where(np.isnan(z_cnnct)))
        for cnnctID_intnlID in cnnctID_intnlID_rch_nan.T:
            cnnctID_intnlID_to_fll = np.delete(cnnctID_intnlID_to_fll, 
            np.intersect1d(
            np.where(cnnctID_intnlID_to_fll[0] == cnnctID_intnlID[0])[0], 
            np.where(cnnctID_intnlID_to_fll[1] == cnnctID_intnlID[1])[0]), 
            axis=1)
        cnnctID_to_fll = \
        np.intersect1d(cnnctID_intnlID_flld[0], cnnctID_intnlID_to_fll[0])

    raise Exception(
    'model {}: filling initial level failed.\nSTOP.'.format(param['nom']))

def write_file_initial_condition(file_ic, absc, z, q, nb_rch):
    """
    Write initial condition file read by solver
    @param file_ic(str) full path of initial water line file
    @param absc(list) absolute abscissas
    @param z(array) water level
    @param q(array) flowrate
    @param nb_rch(int) number of reach in model
    """

    absc = np.array(absc, dtype=np.double)
    nb_fmt = int(absc.size / 5) * 5
    nb_rw = int(nb_fmt / 5)
    nb_mod = absc.size % 5

    with open(file_ic, 'w') as file:
        file.write('RESULTATS CALCUL, DATE :  '
        '{}\n'.format(datetime.now().strftime('%d/%m/%Y %H:%M')))
        file.write('FICHIER RESULTAT MASCARET\n')
        file.write('-------------------------------- \n')
        file.write(' IMAX =%6i ' % absc.size)
        file.write('NBBIEF=%5i\n' % nb_rch)
        file.write(((int((nb_rch - 1) / 5) + 1) * ' ENTETE NON RELUE\n'))

        np.savetxt(file, absc[0:nb_fmt].reshape(nb_rw, 5), 
        fmt=(5 * '%15.3f'), header='X', comments=' ')
        if nb_mod != 0:
            np.savetxt(file, absc[nb_fmt:].reshape(1, nb_mod), 
            fmt=(nb_mod * '%15.3f'))
        np.savetxt(file, z[0:nb_fmt].reshape(nb_rw, 5), 
        fmt=(5 * '%15.3f'), header='Z', comments=' ')
        if nb_mod != 0:
            np.savetxt(file, z[nb_fmt:].reshape(1, nb_mod), 
            fmt=(nb_mod * '%15.3f'))
        np.savetxt(file, q[0:nb_fmt].reshape(nb_rw, 5), 
        fmt=(5 * '%15.3f'), header='Q', comments=' ')
        if nb_mod != 0:
            np.savetxt(file, q[nb_fmt:].reshape(1, nb_mod), 
            fmt=(nb_mod * '%15.3f'))

        file.write(' FIN')

def run_to_steady(chnl, file_ic_tmp, param, bc1, bc2, 
                  z_in=None, q_in=None, IDT=None):
    """
    Run for initial condition depending on boundary condition
    @param chnl(obj) instance of channel model class
    @param file_ic_tmp(str) full path of temporary initial condition file
    @param param(dict) loaded definitions of caller
    @param bc1, bc2(2d array) prescribed initial boundary condition
    @param z_in(1d array) input initial water level
    @param q_in(1d array) input initial flowrate
    @param IDT(list) index of input model cross-section
    """

    def set_input_initial_condition(setID1_z, setID2_z, setID1_q, setID2_q):
        """
        Set input initial condition into model
        @param setID1_z(array) cross-section ID to set water level
        @param setID2_z(array) computational node ID to set water level
        @param setID1_q(array) cross-section ID to set flowrate
        @param setID2_q(array) computational node ID to set flowrate
        """

        for i in range(setID1_z.size):
            chnl.masc.set('State.Z', z_in[setID1_z[i]], setID2_z[i])

        for i in range(setID1_q.size):
            chnl.masc.set('State.Q', q_in[setID1_q[i]], setID2_q[i])

    def run_compute(time, dt):
        """
        Run model computation
        @param time(float) model time
        @param dt(float) timestep
        @return time(float) updated model time
        @return rslt_updt(array) updated model result
        """

        chnl.masc.compute_bc(time, time + dt, dt, 
        np.array([time, (time + dt)]), 2, bc1[0, :].size, bc1, bc2)
        time = chnl.masc.get('State.PreviousTime')
        if time > param['prd_run_to_stdy']:
            raise Exception('model {}: '
            'configured run period exceeded.\nSTOP.'.format(param['nom']))
        rslt_updt = np.array(
        [[chnl.masc.get('State.Z', i) for i in range(chnl.model_size)], 
        [chnl.masc.get('State.Q', i) for i in range(chnl.model_size)]])

        return time, rslt_updt

    os.chdir(chnl.paths['output'])
    chnl.initialize_model(lig=file_ic_tmp)
    os.chdir(chnl.paths['study'])

    zbot = np.array(
    [chnl.masc.get('Model.Zbot', i) for i in range(chnl.model_size)])
    chnl.masc.set('Model.InitTime', 0.0)
    time = 0.0
    rslt = np.array(
    [[chnl.masc.get('State.Z', i) for i in range(chnl.model_size)], 
    [chnl.masc.get('State.Q', i) for i in range(chnl.model_size)]])
    stdy = False

    if not z_in is None and not q_in is None:
        xsID_bdrys = []
        for ndID_bdry in [chnl.masc.get('Model.Connect.NodeNumFreeOutflow', i) 
        for i in range(
        chnl.masc.get_var_size('Model.Connect.NodeNumFreeOutflow')[0])]:
            xsID_bdrys.append(IDT[(ndID_bdry - 1)] - 1)

        setID1_z = np.delete(np.array([i for i in range(z_in.size)]), 
        np.append(np.where(np.isnan(z_in))[0], (np.array(xsID_bdrys))))
        setID2_z = np.array([i for i in range(chnl.model_size)])\
        [np.searchsorted(np.array(IDT), (setID1_z + 1), side='left')]

        setID1_q = np.delete(np.array([i for i in range(q_in.size)]), 
        np.append(np.where(np.isnan(q_in))[0], (np.array(xsID_bdrys))))
        setID2_q = np.array([i for i in range(chnl.model_size)])\
        [np.searchsorted(np.array(IDT), (setID1_q + 1), side='left')]

    """
    If variable timestep, timestep is gotten each step. 
    Otherwise, the timestep stated in case file is used 
    throughout overall computation.
    """
    if chnl.masc.get('Model.VarTimeStep'):
        """
        API term "State.DT" links to 
        varialbe "DT" inside type "ETAT_MASCARET_T".
        Initially, varialbe "DT" is set to 0.0, 
        but will be adjusted to 
        timestep stated in case file before computation if it is 0.0. 
        Then varialbe "DT" for next run is always 
        updated during current run after 
        computation of hydraulic variables.
        """
        dt = chnl.masc.get('Model.DT')

        if not z_in is None and not q_in is None:
            set_input_initial_condition(
            setID1_z, setID2_z, setID1_q, setID2_q)

        time, rslt_updt = run_compute(time, dt)

        if np.all(
        np.isclose((rslt_updt[0, :] - zbot), (rslt[0, :] - zbot), 
        rtol=param['tlrc_rltv_dpth_run_to_stdy'], 
        atol=param['tlrc_abslt_dpth_run_to_stdy'], equal_nan=False)) \
        and np.all(np.isclose(rslt_updt[1, :], rslt[1, :], 
        rtol=param['tlrc_rltv_flrt_run_to_stdy'], 
        atol=param['tlrc_abslt_flrt_run_to_stdy'], equal_nan=False)):
            stdy = True
        rslt = rslt_updt

        while not stdy:
            dt = chnl.masc.get('State.DT')

            if not z_in is None and not q_in is None:
                set_input_initial_condition(
                setID1_z, setID2_z, setID1_q, setID2_q)

            time, rslt_updt = run_compute(time, dt)

            if np.all(
            np.isclose((rslt_updt[0, :] - zbot), (rslt[0, :] - zbot), 
            rtol=param['tlrc_rltv_dpth_run_to_stdy'], 
            atol=param['tlrc_abslt_dpth_run_to_stdy'], equal_nan=False)) \
            and np.all(np.isclose(rslt_updt[1, :], rslt[1, :], 
            rtol=param['tlrc_rltv_flrt_run_to_stdy'], 
            atol=param['tlrc_abslt_flrt_run_to_stdy'], equal_nan=False)):
                stdy = True
            rslt = rslt_updt
        """
        could be steady but different from input
        """
        if not z_in is None and not q_in is None:
            set_input_initial_condition(
            setID1_z, setID2_z, setID1_q, setID2_q)
    else:
        dt = chnl.masc.get('Model.DT')
        while not stdy:
            if not z_in is None and not q_in is None:
                set_input_initial_condition(
                setID1_z, setID2_z, setID1_q, setID2_q)

            time, rslt_updt = run_compute(time, dt)

            if np.all(
            np.isclose((rslt_updt[0, :] - zbot), (rslt[0, :] - zbot), 
            rtol=param['tlrc_rltv_dpth_run_to_stdy'], 
            atol=param['tlrc_abslt_dpth_run_to_stdy'], equal_nan=False)) \
            and np.all(np.isclose(rslt_updt[1, :], rslt[1, :], 
            rtol=param['tlrc_rltv_flrt_run_to_stdy'], 
            atol=param['tlrc_abslt_flrt_run_to_stdy'], equal_nan=False)):
                stdy = True
            rslt = rslt_updt
        """
        could be steady but different from input
        """
        if not z_in is None and not q_in is None:
            set_input_initial_condition(
            setID1_z, setID2_z, setID1_q, setID2_q)

def main():
    """
    Main function
    """

    param = get_param()

    if "config_file" in param.keys():
        cfg_file = param['config_file']
    else:
        cfg_file = os.environ.get('SYSTELCFG')
        if cfg_file is None: cfg_file = ""
    if "config_name" in param.keys():
        cfg_name = param['config_name']
    else:
        cfg_name = os.environ.get('USETELCFG')
        if cfg_name is None: cfg_name = ""
    root_dir = os.environ.get('HOMETEL')
    if root_dir is None:
        root_dir = ""
        python_dir = ""
    else:
        python_dir = os.path.join(root_dir, 'scripts', 'python3')
    CFGS.parse_cfg_file(cfg_file, cfg_name, root_dir, python_dir)

    drct_ic = os.path.dirname(param['config']['files']['lig'])

    iprint = 1 if param['iprint'] else 0
    chnl = MascaretStudy(param['config'], 
    log_lvl='CRITICAL', iprint=iprint, working_directory='.')
    nb_rch = chnl.masc.get_var_size('Model.Connect.FirstNdNum')[0]
    nb_xs = len(chnl.cross_sections['bottom']['abscissa'])
    """
    The solver only accepts sorted IDT in order of natural numbers from 1.
    The last value of "Model.IDT" is false
    """
    IDT = [chnl.masc.get('Model.IDT', i) for i in range(chnl.model_size)]
    IDT[-1] = len(chnl.cross_sections['bottom']['abscissa'])
    param['tm_stt'] = chnl.masc.get('Model.InitTime')

    file_ic = os.path.join(drct_ic, 
    'WaterLine_{}_{}.lig'.format(param['nom'], param['tm_stt']))

    if param['run_to_stdy']:
        file_ic_tmp = os.path.join(chnl.paths['data'], 'WaterLine_tmp.lig')

        bc1_intl, bc2_intl = \
        get_boundary_condition_cas_initial(chnl, param['tm_stt'])
        """
        only integrate 1 time step
        """
        bc1_cmpt = np.zeros([2, len(bc1_intl.keys())], dtype=np.double)
        bc2_cmpt = np.zeros([2, len(bc1_intl.keys())], dtype=np.double)
        for bcID_1 in bc1_intl.keys():
            bc1_cmpt[:, (bcID_1 - 1)] = bc1_intl[bcID_1]
            if bcID_1 in bc2_intl.keys():
                bc2_cmpt[:, (bcID_1 - 1)] = bc2_intl[bcID_1]

    if "ic_cold" in param['config']['files'].keys():
        z_psd, q_psd = parse_initial_condition_prescribed(
        os.path.join(chnl.paths['data'], 
        os.path.basename(param['config']['files']['ic_cold'])), 
        chnl.cross_sections['bottom']['abscissa'])

        z_rvsd = np.zeros(nb_xs, dtype=np.double)
        q_rvsd = np.zeros(nb_xs, dtype=np.double)
        z_rvsd[:] = z_psd[:]
        q_rvsd[:] = q_psd[:]

        zID_nan = np.where(np.isnan(z_rvsd))[0]
        if zID_nan.size:
            """
            different approach of revising level for 
            kernels capable of trans-critical flow or not
            """
            if chnl.masc.get('Model.Kernel') == 2:
                bdry_tp = [chnl.masc.get('Model.Boundary.Type', i) for i in 
                range(chnl.masc.get_var_size('Model.Boundary.Name')[0])]
                """
                For case of prescribed level at boundary, 
                need to avoid super-critical at boundary.
                """
                if 2 in bdry_tp:
                    if not param['run_to_stdy']:
                        bc1_intl, _ = get_boundary_condition_cas_initial(
                        chnl, param['tm_stt'])
                    z_fll = fill_initial_level(
                    chnl, param, bdry_tp, IDT, bc1_intl, zID_nan)
                    z_rvsd[zID_nan] = z_fll[zID_nan]
                else:
                    for xsID in zID_nan:
                        z_rvsd[xsID] = min(
                        [chnl.masc.get('Model.CrossSection.Zbank', xsID, i) 
                        for i in range(2)])
            elif chnl.masc.get('Model.Kernel') == 3:
                for xsID in zID_nan:
                    z_rvsd[xsID] = min(
                    [chnl.masc.get('Model.CrossSection.Zbank', xsID, i) 
                    for i in range(2)])
        for xsID in np.where(np.isnan(q_rvsd))[0]:
            q_rvsd[xsID] = 0.0
        if np.any(np.isnan(z_rvsd)) or np.any(np.isnan(q_rvsd)):
            raise ValueError(
            'model {}: supplementing initial conditions failed.\nSTOP.'
            .format(param['nom']))

        if param['run_to_stdy']:
            write_file_initial_condition(file_ic_tmp, 
            chnl.cross_sections['bottom']['abscissa'], z_rvsd, q_rvsd, nb_rch)
            run_to_steady(
            chnl, file_ic_tmp, param, bc1_cmpt, bc2_cmpt, z_psd, q_psd, IDT)
            write_file_initial_condition(file_ic, chnl.xcoord, 
            np.array([chnl.masc.get('State.Z', i) 
            for i in range(chnl.model_size)]), 
            np.array([chnl.masc.get('State.Q', i) 
            for i in range(chnl.model_size)]), 
            nb_rch)
        else:
            write_file_initial_condition(file_ic, 
            chnl.cross_sections['bottom']['abscissa'], z_rvsd, q_rvsd, nb_rch)
    else:
        z = np.zeros(nb_xs, dtype=np.double)
        """
        different approach of revising level for 
        kernels capable of trans-critical flow or not
        """
        if chnl.masc.get('Model.Kernel') == 2:
            bdry_tp = [chnl.masc.get('Model.Boundary.Type', i) for i in 
            range(chnl.masc.get_var_size('Model.Boundary.Name')[0])]
            """
            For case of prescribed level at boundary, 
            need to avoid super-critical at boundary.
            """
            if 2 in bdry_tp:
                if not param['run_to_stdy']:
                    bc1_intl, _ = \
                    get_boundary_condition_cas_initial(chnl, param['tm_stt'])
                z = fill_initial_level(chnl, param, bdry_tp, IDT, bc1_intl)
            else:
                for xsID in range(nb_xs):
                    z[xsID] = min(
                    [chnl.masc.get('Model.CrossSection.Zbank', xsID, i) 
                    for i in range(2)])
        elif chnl.masc.get('Model.Kernel') == 3:
            for xsID in range(nb_xs):
                z[xsID] = min(
                [chnl.masc.get('Model.CrossSection.Zbank', xsID, i) 
                for i in range(2)])
        q = np.zeros(nb_xs, dtype=np.double)

        if param['run_to_stdy']:
            write_file_initial_condition(file_ic_tmp, 
            chnl.cross_sections['bottom']['abscissa'], z, q, nb_rch)
            run_to_steady(chnl, file_ic_tmp, param, bc1_cmpt, bc2_cmpt)
            write_file_initial_condition(file_ic, chnl.xcoord, 
            np.array([chnl.masc.get('State.Z', i) 
            for i in range(chnl.model_size)]), 
            np.array([chnl.masc.get('State.Q', i) 
            for i in range(chnl.model_size)]), 
            nb_rch)
        else:
            write_file_initial_condition(file_ic, 
            chnl.cross_sections['bottom']['abscissa'], z, q, nb_rch)

    del chnl

    sys.exit(0)

if __name__ == "__main__":
    main()