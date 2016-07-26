# CDmetaPOP.py
# Author: Erin L Landguth
# Created: February 2008
# v 1.0 Release: MARCH 2014
# ----------------------------------------------------------------------------
# General CDmetaPOP information
appName = "CDmetaPOP"
appVers = "version 1.03"
appRele = "2016.04.22-10:23:01MDT"
authorNames = "Erin L Landguth"

# ---------------
# Global symbols
#----------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False

# ------------------------------------------
# Import Modules with Except/Try statements
# ------------------------------------------
# Python specific functions
import datetime,time,pdb,os,sys,shutil,gc,multiprocessing

# Numpy functions
try:
    import numpy as np
except ImportError as eMsg:
    print("ImportError (%s) Numpy required."%(eMsg))
    sys.exit(-1)

# CDmetaPOP functions
from CDmetaPOP.modules import *
from CDmetaPOP.postprocess import *
from CDmetaPOP.preprocess import *


def doBatch(ibatch, params):
    # Timing events: start
    start_timeB = datetime.datetime.now()

    # Store all information and the type of each, also do some error checks
    xyfilename = datadir+params['xyfilename']
    #agefilename = datadir+params['agefilename']

    # Grab the nthfile list range specific to user input, list or sequence
    if not isinstance(params['output_years'], (list,tuple)):
        params['output_years'] = int(params['output_years'])
        if params['output_years'] != 0:
            nthfile = range(0,params['runtime']+params['output_years'],params['output_years'])
            del(nthfile[-1]) # Delete the last value 0, looptime - 1
        else:
            nthfile = [0]
    # If specified years with |
    else:
        nthfile = []
        # Split up list, removing space values, and appending to nthfile
        for inum in xrange(len(params['output_years'])):
            # Error check here if | at the end
            if len(params['output_years'][inum]) != 0:
                nthfile.append(int(params['output_years'][inum]))

    # Error check on nthfile, must be 1 less than looptime for indexing
    if max(nthfile) >= params['runtime']:
        print('nthfile selection maximum value must be less than to looptime.')
        sys.exit(-1)

    # Store cdmat file information - header file (loadFile()) passes tuple or string if only 1
    if not isinstance(params['cdclimgentime'], (list,tuple)):
        cdclimgentime = [params['cdclimgentime']]
    else:
        cdclimgentime = params['cdclimgentime']


    # ---------------------------------------------
    # Begin Monte-Carlo Looping
    # ---------------------------------------------

    # xrange(mcruns) is typically 10 - 50...and it takes a long time.
    for ithmcrun in xrange(params['mcruns']):

        # Timing events: start
        start_timeMC = datetime.datetime.now()

        # -----------------------------------------
        # Create storage variables
        # ------------------------------------------
        # These variables will be stored in output.csv at the end of the simulation

        # GetMetrics()
        Track_p1 = []
        Track_p2 = []
        Track_q1 = []
        Track_q2 = []
        Track_Alleles = []
        Track_He = []
        Track_Ho = []
        Track_N_Init_pop = []
        Track_N_Init_age = []
        Track_N_Init_class = []
        Track_K = []
        Track_CaptureCount_Out = []
        Track_CaptureCount_ClassOut = []
        Track_CaptureCount_Back = []
        Track_CaptureCount_ClassBack = []
        Track_MatureCount = []
        Track_ImmatureCount = []

        # DoMate()
        Track_FAvgMate = []
        Track_MAvgMate = []
        Track_FSDMate = []
        Track_MSDMate = []
        Track_MateDistCD = []
        Track_MateDistCDstd = []
        Track_ToTFemales = []
        Track_ToTMales = []
        Track_BreedFemales = []
        Track_BreedMales = []
        Track_BreedEvents = []

        # DoOffspring
        Track_Births = []
        Track_EggDeaths = []

        # DoUpdate
        Track_N_back_age = []
        Track_N_out_age = []

        # Emigration()
        N_Emigration_pop = []
        N_Emigration_age = []
        F_EmiDist = []
        M_EmiDist = []
        F_EmiDist_sd = []
        M_EmiDist_sd = []
        SelectionDeathsEmi = []
        DisperseDeathsEmi = []
        PackingDeathsEmi = []
        PackingDeathsEmiAge = []
        MgSuccess = []
        AdultNoMg = []

        # Mortlity after Emigration
        N_EmiMortality = []
        PopDeathsOUT = []
        AgeDeathsOUT = []
        SizeDeathsOUT = []

        # Immigration
        N_beforePack_Immi_pop = []
        N_beforePack_Immi_age = []
        N_Immigration_pop = []
        N_Immigration_age = []
        F_HomeDist = []
        M_HomeDist = []
        F_HomeDist_sd = []
        M_HomeDist_sd = []
        F_StrayDist = []
        M_StrayDist = []
        F_StrayDist_sd = []
        M_StrayDist_sd = []
        F_ZtrayDist = []
        M_ZtrayDist = []
        F_ZtrayDist_sd = []
        M_ZtrayDist_sd = []
        SelectionDeathsImm = []
        SelectionDeathsImm_Age0s = []
        DisperseDeathsImm = []
        PackingDeathsImmAge = []
        PackingDeathsImm = []
        StrSuccess = []

        # Mortality after immigration
        N_ImmiMortality = []
        PopDeathsIN = []
        AgeDeathsIN = []
        SizeDeathsIN = []

        # DoOutput()
        Infected = []
        Residors = []
        Strayers1 = []
        Strayers2 = []
        Immigrators = []
        PopSizes_Mean = []
        PopSizes_Std = []
        AgeSizes_Mean = []
        AgeSizes_Std = []
        ClassSizes_Mean = []
        ClassSizes_Std = []

        # ------------------------------------
        # Call DoPreProcess()
        # ------------------------------------

        # Timing events: start
        start_time1 = datetime.datetime.now()

        # Call function
        preprocess_params = DoPreProcess(outdir, datadir, ibatch, ithmcrun, xyfilename, params)

        # Grab first one only
        K = preprocess_params["K_mu"] # Initialize K with mu

        # Print to log
        stringout = 'DoPreProcess(): '+str(datetime.datetime.now() -start_time1) + ''
        logMsg(logfHndl,stringout)
        print(stringout)

        # ---------------------------------
        # Call GetMetrics()
        # ---------------------------------

        # Timing events: start
        start_time1 = datetime.datetime.now()

        GetMetrics(preprocess_params["SubpopIN_init"], K, Track_N_Init_pop, Track_K, params['loci'], params['alleles'],
                   0, Track_Ho,
                   Track_Alleles,
                   Track_He,
                   Track_p1, Track_p2, Track_q1, Track_q2, Infected, Residors, Strayers1, Strayers2, Immigrators,
                   PopSizes_Mean, PopSizes_Std, AgeSizes_Mean, AgeSizes_Std, Track_ToTMales, Track_ToTFemales,
                   Track_BreedMales, Track_BreedFemales, Track_N_Init_age, Track_MatureCount, Track_ImmatureCount,
                   params['sizecontrol'], preprocess_params["age_size_mean"], ClassSizes_Mean, ClassSizes_Std,
                   Track_N_Init_class,
                   params['sexans'], params['SNPanswer'])

        # Print to log
        stringout = 'GetMetrics() Initial: '+str(datetime.datetime.now() -start_time1) + ''
        logMsg(logfHndl,stringout)
        print(stringout)

        # ---------------------------------
        # Error statements
        # ---------------------------------
        # Error statement here in case no females or males, then break
        if Track_ToTFemales[0][0]==0 or Track_ToTMales[0][0]==0:
            print('There are no females or males to begin time loop.\n')
            break

        # ------------------------------------------
        # Call DoUpdate() - output initial file here ind-1.csv
        # ------------------------------------------

        # Timing events: start
        start_time1 = datetime.datetime.now()

        DoUpdate(preprocess_params["SubpopIN_init"], K, preprocess_params["xgridpop"], preprocess_params["ygridpop"],
                 -1, nthfile,
                 preprocess_params["ithmcrundir"], params['loci'],
                 params['alleles'],
                 logfHndl,
                 'Initial', 'N', 'N', [], params['startSelection'], [], [], [], [], [], [])

        # Print to log
        stringout = 'DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
        logMsg(logfHndl,stringout)
        print 'First DoUpdate(): ',str(datetime.datetime.now() -start_time1),''

        # -------------------------------------------
        # Start Generation Looping
        # -------------------------------------------
        # Begin generation loop
        for gen in xrange(params['runtime']):

            # Timing events: start
            start_timeGen = datetime.datetime.now()

            # If initial generation - update with initial populations
            if gen == 0:
                SubpopIN = preprocess_params["SubpopIN_init"]
                del preprocess_params["SubpopIN_init"]
                # Use NatalPop
                sourcePop = 'NatalPop'
            else:
                sourcePop = 'ImmiPop'

            # Exit the system if population is 0 or 1
            if sum(preprocess_params["N0"]) <= 1:
                stringout = 'Population went extinct, program ended.'
                logMsg(logfHndl, stringout)
                print('Population went extinct after generation ' + str(gen - 1) + '.\n')
                break

            # ---------------------------------
            # Call CDClimate()
            # ---------------------------------

            # Timing events: start
            start_time1 = datetime.datetime.now()

            # Check gen time equal to cdclimgentime
            for icdtime in xrange(len(cdclimgentime)):
                if gen == int(cdclimgentime[icdtime]):
                    climate_params = DoCDClimate(datadir, icdtime, cdclimgentime, params['mate_cdmat'],
                                             params['dispout_cdmat'],
                                             params['dispback_cdmat'], params['stray_cdmat'], params['matemoveno'],
                                             params['FdispmoveOutno'],
                                             params['MdispmoveOutno'], params['FdispmoveBackno'],
                                             params['MdispmoveBackno'],
                                             params['StrayBackno'],
                                             params['matemovethresh'], params['FdispmoveOutthresh'],
                                             params['MdispmoveOutthresh'],
                                             params['FdispmoveBackthresh'], params['MdispmoveBackthresh'],
                                             params['StrayBackthresh'],
                                             params['matemoveparA'], params['matemoveparB'], params['matemoveparC'],
                                             params['FdispmoveOutparA'],
                                             params['FdispmoveOutparB'], params['FdispmoveOutparC'],
                                             params['MdispmoveOutparA'],
                                             params['MdispmoveOutparB'],
                                             params['MdispmoveOutparC'],
                                             params['FdispmoveBackparA'],
                                             params['FdispmoveBackparB'], params['FdispmoveBackparC'],
                                             params['MdispmoveBackparA'],
                                             params['MdispmoveBackparB'],
                                             params['MdispmoveBackparC'], params['StrayBackparA'],
                                             params['StrayBackparB'],
                                             params['StrayBackparC'],
                                             preprocess_params["Mg_pass"],
                                             preprocess_params["Str_pass"], preprocess_params["Kmu_pass"],
                                             preprocess_params["outsizevals_pass"],
                                             preprocess_params["backsizevals_pass"],
                                             preprocess_params["outgrowdays_pass"],
                                             preprocess_params["backgrowdays_pass"], preprocess_params["fitvals_pass"],
                                             preprocess_params["popmort_back_pass"],
                                             preprocess_params["popmort_out_pass"],
                                             preprocess_params["eggmort_pass"], preprocess_params["Kstd_pass"],
                                             preprocess_params["popmort_back_sd_pass"],
                                             preprocess_params["popmort_out_sd_pass"],
                                             preprocess_params["eggmort_sd_pass"],
                                             preprocess_params["outsizevals_sd_pass"],
                                             preprocess_params["backsizevals_sd_pass"],
                                             preprocess_params["outgrowdays_sd_pass"],
                                             preprocess_params["backgrowdays_sd_pass"],
                                             preprocess_params["pop_capture_back_pass"],
                                             preprocess_params["pop_capture_out_pass"], params['cdevolveans'],
                                             preprocess_params["N0_pass"],
                                             preprocess_params["allefreqfiles_pass"],
                                             preprocess_params["classvarsfiles_pass"])

                    # ----------------------------------------
                    # Introduce new individuals
                    # ----------------------------------------
                    if (gen != 0 and len(preprocess_params["N0_pass"][0].split('|')) > 1):
                        SubpopIN = AddIndividuals(SubpopIN, climate_params["tempN0"], climate_params["tempAllelefile"],
                                                  climate_params["tempClassVarsfile"],
                                                  datadir,
                                                  params['loci'],
                                                  params['alleles'], params['sizecontrol'], params['cdinfect'],
                                                  params['SNPanswer'],
                                                  params['cdevolveans'], params['startSelection'],
                                                  climate_params["fitvals"],
                                                  params['eggFrequency'], params['mature_length_female'],
                                                  params['mature_length_male'],
                                                  params['mature_int_female'], params['mature_slope_female'],
                                                  params['mature_int_male'], params['mature_slope_male'],
                                                  preprocess_params["dtype"],
                                                  preprocess_params["N0"],
                                                  preprocess_params["natal"], gen)
            # -------------------------------------------
            # Update stochastic parameters each year here
            # -------------------------------------------
            tupStoch = DoStochasticUpdate(climate_params["K_mu"], climate_params["K_std"], climate_params["popmort_back_mu"],
                                          climate_params["popmort_back_sd"],
                                          climate_params["popmort_out_mu"],
                                          climate_params["popmort_out_sd"],
                                          climate_params["eggmort_mu"], climate_params["eggmort_sd"],
                                          climate_params["outsizevals_mu"],
                                          climate_params["outsizevals_sd"],
                                          climate_params["backsizevals_mu"],
                                          climate_params["backsizevals_sd"], climate_params["outgrowdays_mu"],
                                          climate_params["outgrowdays_sd"],
                                          climate_params["backgrowdays_mu"],
                                          climate_params["backgrowdays_sd"], preprocess_params["age_percmort_out_mu"],
                                          preprocess_params["age_percmort_out_sd"],
                                          preprocess_params["age_percmort_back_mu"],
                                          preprocess_params["age_percmort_back_sd"],
                                          preprocess_params["size_percmort_out_mu"],
                                          preprocess_params["size_percmort_out_sd"],
                                          preprocess_params["size_percmort_back_mu"],
                                          preprocess_params["size_percmort_back_sd"],
                                          params['Egg_Mortality'], params['Egg_Mortality_StDev'],
                                          preprocess_params["cor_mat"])
            K = tupStoch[0]
            popmort_back = tupStoch[1]
            popmort_out = tupStoch[2]
            eggmort_patch = tupStoch[3]
            outsizevals = tupStoch[4]
            backsizevals = tupStoch[5]
            outgrowdays = tupStoch[6]
            backgrowdays = tupStoch[7]
            age_percmort_out = tupStoch[8]
            age_percmort_back = tupStoch[9]
            size_percmort_out = tupStoch[10]
            size_percmort_back = tupStoch[11]
            eggmort_pop = tupStoch[12]

            # Print to log
            stringout = 'DoCDClimate(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'DoCDClimate(): ',str(datetime.datetime.now() -start_time1),''
            # ---------------------------------------
            # Call DoMate()
            # ---------------------------------------

            # Timing events: start
            start_time1 = datetime.datetime.now()

            Bearpairs = DoMate(SubpopIN, K, \
                               params['Freplace'], params['Mreplace'], climate_params["mateno"], climate_params["thresh_mate"], \
                               climate_params["cdmatrix_mate"], Track_MateDistCD, preprocess_params["xgridpop"], \
                               preprocess_params["ygridpop"], Track_MateDistCDstd, Track_FAvgMate, Track_MAvgMate,
                               Track_FSDMate,
                               Track_MSDMate, Track_BreedEvents, gen, sourcePop, preprocess_params["dtype"],
                               climate_params["mate_ScaleMax"], climate_params["mate_ScaleMin"],
                               params['matemoveparA'], params['matemoveparB'], params['matemoveparC'],
                               params['Egg_FemalePercent'], params['eggFrequency'],
                               params['sexans'], params['selfans'])

            # Print to log
            stringout = 'DoMate(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'DoMate(): ',str(datetime.datetime.now() -start_time1),''

            # ---------------------------------------
            # Call DoOffspring()
            # ---------------------------------------

            # Timing events: start
            start_time1 = datetime.datetime.now()

            noOffspring, Bearpairs = DoOffspring(params['offno'], Bearpairs, \
                                                 Track_Births, params['transmissionprob'], gen, K, sourcePop, \
                                                 preprocess_params["age_mu"], preprocess_params["age_sigma"],
                                                 params['sizecontrol'], \
                                                 params['Egg_Mean_par1'], params['Egg_Mean_par2'],
                                                 params['Egg_Mean_ans'],
                                                 params['equalClutchSize'], preprocess_params["dtype"],
                                                 params['mature_length_male'], params['mature_length_female'],
                                                 eggmort_patch, Track_EggDeaths,
                                                 eggmort_pop)

            # Print to log
            stringout = 'DoOffspring(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'DoOffspring(): ',str(datetime.datetime.now() -start_time1),''

            # ----------------------------------------------------------------
            # Call 2nd DoUpdate() - grow, age/mature (selection option),egglay,capture, output ind.csv file;no Age0s; ind.csv
            # ----------------------------------------------------------------

            # Timing events: start
            start_time1 = datetime.datetime.now()

            SubpopIN = DoUpdate(SubpopIN, K, preprocess_params["xgridpop"], preprocess_params["ygridpop"], gen, nthfile,
                                preprocess_params["ithmcrundir"],
                                params['loci'],
                                params['alleles'],
                                logfHndl,
                                'Middle', params['growth_option'], params['cdevolveans'], climate_params["fitvals"],
                                params['startSelection'],
                                preprocess_params["age_capture_back"],
                                climate_params["pop_capture_back"],
                                Track_CaptureCount_Back, Track_CaptureCount_ClassBack, params['sizecontrol'],
                                preprocess_params["age_size_mean"],
                                Track_N_back_age,
                                params['eggFrequency'], backsizevals, params['growth_Loo'], params['growth_R0'],
                                params['growth_temp_max'],
                                params['growth_temp_CV'],
                                params['growth_temp_t0'], backgrowdays, sourcePop, params['sizecontrol'],
                                preprocess_params["M_mature"],
                                preprocess_params["F_mature"],
                                params['mature_slope_male'], params['mature_int_male'],
                                params['mature_slope_female'], params['mature_int_female'],
                                params['mature_length_male'],
                                params['mature_length_female'])

            # Print to log
            stringout = 'Second DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'Second DoUpdate(): ',str(datetime.datetime.now() -start_time1),''

            # ------------------------------------------
            # Call DoEmigration()
            # ------------------------------------------

            # Timing events: start
            start_time1 = datetime.datetime.now()

            SubpopIN = DoEmigration(SubpopIN, K, climate_params["FdispOutno"], \
                                    climate_params["MdispOutno"], climate_params["cdmatrix_FOut"], climate_params["cdmatrix_MOut"],
                                    gen,
                                    preprocess_params["xgridpop"],
                                    preprocess_params["ygridpop"], F_EmiDist,
                                    M_EmiDist, params['cdevolveans'], climate_params["fitvals"], F_EmiDist_sd, M_EmiDist_sd,
                                    preprocess_params["subpopemigration"], \
                                    SelectionDeathsEmi, DisperseDeathsEmi, \
                                    params['startSelection'], climate_params["Mg"], \
                                    MgSuccess, AdultNoMg, sum(params['alleles']), preprocess_params["age_Mg"],
                                    climate_params["thresh_FOut"], climate_params["thresh_MOut"],
                                    N_Emigration_pop, sourcePop, preprocess_params["dtype"],
                                    preprocess_params["setmigrate"],
                                    params['sizecontrol'],
                                    preprocess_params["age_size_mean"],
                                    PackingDeathsEmi, N_Emigration_age, params['loci'], params['muterate'],
                                    params['mtdna'],
                                    params['mutationtype'],
                                    climate_params["FdispOut_ScaleMax"], climate_params["FdispOut_ScaleMin"],
                                    climate_params["MdispOut_ScaleMax"],
                                    climate_params["MdispOut_ScaleMin"],
                                    params['FdispmoveOutparA'], params['FdispmoveOutparB'], params['FdispmoveOutparC'],
                                    params['MdispmoveOutparA'], params['MdispmoveOutparB'], params['MdispmoveOutparC'],
                                    params['popmodel'],
                                    PackingDeathsEmiAge,
                                    preprocess_params["ithmcrundir"], params['popmodel_par1'],
                                    params['implementSelection'],
                                    age_percmort_out, preprocess_params["migrate"])

            # Print to log
            stringout = 'DoEmigration(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'DoEmigration(): ',str(datetime.datetime.now() -start_time1),''

            # ----------------------------------------
            # Call DoMortality()
            # ----------------------------------------
            start_time1 = datetime.datetime.now() # Timing events: start

            SubpopIN = DoMortality(SubpopIN, K, PopDeathsOUT, \
                                   popmort_out, age_percmort_out, \
                                   gen, N_EmiMortality, AgeDeathsOUT, params['sizecontrol'],
                                   preprocess_params["age_size_mean"],
                                   size_percmort_out, SizeDeathsOUT, params['constMortans'], params['popmodel'])

            # Print to log
            stringout = 'DoOutMortality(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'DoOutMortality(): ',str(datetime.datetime.now() -start_time1),''

            # ----------------------------------------------------
            # Call DoUpdate() - grow, capture, and optional output indSample.csv
            # ----------------------------------------------------
            start_time1 = datetime.datetime.now()  # Timing events: start

            SubpopIN = DoUpdate(SubpopIN, K, preprocess_params["xgridpop"], preprocess_params["ygridpop"], gen, nthfile,
                                preprocess_params["ithmcrundir"],
                                params['loci'],
                                params['alleles'],
                                logfHndl,
                                params['gridsampling'], params['growth_option'], 'N', [], params['startSelection'],
                                preprocess_params["age_capture_out"],
                                climate_params["pop_capture_out"],
                                Track_CaptureCount_Out, Track_CaptureCount_ClassOut, params['sizecontrol'],
                                preprocess_params["age_size_mean"],
                                Track_N_out_age,
                                params['eggFrequency'], outsizevals, params['growth_Loo'], params['growth_R0'],
                                params['growth_temp_max'],
                                params['growth_temp_CV'],
                                params['growth_temp_t0'], outgrowdays, 'EmiPop')

            # Print to log
            stringout = 'Third DoUpdate(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'Third DoUpdate(): ',str(datetime.datetime.now() -start_time1),''

            # ------------------------------------------
            # Call DoImmigration()
            # ------------------------------------------
            start_time1 = datetime.datetime.now()  # Timing events: start

            SubpopIN = DoImmigration(SubpopIN, K, preprocess_params["N0"], preprocess_params["natal"],
                                     climate_params["FdispBackno"], \
                                     climate_params["MdispBackno"], climate_params["cdmatrix_FBack"],
                                     climate_params["cdmatrix_MBack"], gen, \
                                     preprocess_params["xgridpop"], preprocess_params["ygridpop"],
                                     params['cdevolveans'], climate_params["fitvals"],
                                     preprocess_params["subpopimmigration"], \
                                     SelectionDeathsImm, DisperseDeathsImm, params['startSelection'], climate_params["Str"], \
                                     StrSuccess, \
                                     climate_params["Strno"], climate_params["cdmatrix_StrBack"], preprocess_params["age_S"],
                                     climate_params["thresh_FBack"],
                                     climate_params["thresh_MBack"],
                                     climate_params["thresh_Str"],
                                     N_Immigration_pop, preprocess_params["dtype"], params['sizecontrol'],
                                     preprocess_params["age_size_mean"], PackingDeathsImm,
                                     N_Immigration_age, climate_params["FdispBack_ScaleMax"],
                                     climate_params["FdispBack_ScaleMin"],
                                     climate_params["MdispBack_ScaleMax"],
                                     climate_params["MdispBack_ScaleMin"],
                                     params['FdispmoveBackparA'], params['FdispmoveBackparB'],
                                     params['FdispmoveBackparC'],
                                     params['MdispmoveBackparA'], params['MdispmoveBackparB'],
                                     params['MdispmoveBackparC'],
                                     climate_params["Str_ScaleMax"],
                                     climate_params["Str_ScaleMin"], params['StrayBackparA'], params['StrayBackparB'],
                                     params['StrayBackparC'],
                                     params['popmodel'],
                                     PackingDeathsImmAge,
                                     preprocess_params["ithmcrundir"], params['popmodel_par1'], noOffspring, Bearpairs,
                                     preprocess_params["age_size_std"],
                                     params['Egg_FemalePercent'],
                                     sourcePop, params['transmissionprob'], preprocess_params["M_mature"],
                                     preprocess_params["F_mature"],
                                     params['mature_slope_male'], params['mature_int_male'],
                                     params['mature_slope_female'],
                                     params['mature_int_female'], params['mature_length_male'],
                                     params['mature_length_female'],
                                     params['loci'],
                                     params['muterate'], params['mtdna'],
                                     params['mutationtype'],
                                     params['startGenes'],
                                     preprocess_params["allelst"], params['HomeAttempt'], params['implementSelection'],
                                     N_beforePack_Immi_pop,
                                     N_beforePack_Immi_age,
                                     SelectionDeathsImm_Age0s, F_StrayDist, M_StrayDist, F_StrayDist_sd, M_StrayDist_sd,
                                     F_ZtrayDist, M_ZtrayDist, F_ZtrayDist_sd, M_ZtrayDist_sd, F_HomeDist, M_HomeDist,
                                     F_HomeDist_sd, M_HomeDist_sd, params['SNPanswer'])
            del(Bearpairs)
            # Print to log
            stringout = 'DoImmigration(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'DoImmigration(): ',str(datetime.datetime.now() -start_time1),''

            # ------------------------------------------
            # Call DoMortality()
            # ------------------------------------------

            # Timing events: start
            start_time1 = datetime.datetime.now()

            SubpopIN = DoMortality(SubpopIN, K, PopDeathsIN, \
                                   popmort_back, age_percmort_back, \
                                   gen, N_ImmiMortality, AgeDeathsIN, params['sizecontrol'],
                                   preprocess_params["age_size_mean"],
                                   size_percmort_back, SizeDeathsIN,
                                   params['constMortans'], params['popmodel'])

            # Print to log
            stringout = 'DoInMortality(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'DoInMortality(): ',str(datetime.datetime.now() -start_time1),''

            # ---------------------------------
            # Call GetMetrics()
            # ---------------------------------

            # Timing events: start
            start_time1 = datetime.datetime.now()

            GetMetrics(SubpopIN, K, Track_N_Init_pop, Track_K, params['loci'], params['alleles'], gen + 1, Track_Ho,
                       Track_Alleles,
                       Track_He, Track_p1, Track_p2, Track_q1, Track_q2, Infected, Residors, Strayers1, Strayers2,
                       Immigrators, PopSizes_Mean, PopSizes_Std, AgeSizes_Mean, AgeSizes_Std, Track_ToTMales,
                       Track_ToTFemales, Track_BreedMales, Track_BreedFemales, Track_N_Init_age, Track_MatureCount,
                       Track_ImmatureCount, params['sizecontrol'], preprocess_params["age_size_mean"], ClassSizes_Mean, ClassSizes_Std, Track_N_Init_class,
                       params['sexans'], params['SNPanswer'])

            # Print to log
            stringout = 'GetMetrics(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'GetMetrics(): ',str(datetime.datetime.now() -start_time1),''

            # Error statement here in case no females or males, then break
            if Track_N_Init_pop[gen+1][0] == 0:
                print('Population went extinct after year'+str(gen)+'.\n')
                break
            if Track_ToTFemales[gen+1][0]==0 or Track_ToTMales[gen+1][0]==0:
                print('There are no more females or males left in population after year '+str(gen)+'.\n')
                break

            # Print to log
            stringout = 'End Generation/Year Loop'+str(gen)+': '+str(datetime.datetime.now() -start_timeGen) + '\n'
            logMsg(logfHndl,stringout)
            print 'End Generation/Year Loop',str(gen),': ',str(datetime.datetime.now() -start_timeGen),'\n'
        # End::generation loop

        # ------------------------------------------
        # Call DoPostProcess()
        # ------------------------------------------

        # Timing events: start
        start_time1 = datetime.datetime.now()

        DoPostProcess(preprocess_params["ithmcrundir"], params['gendmatans'], params['loci'], params['alleles'],
                      params['runtime'], \
                      Track_ToTFemales, Track_ToTMales, Track_BreedFemales, Track_BreedMales, Track_Births, PopDeathsIN, \
                      PopDeathsOUT, Track_Alleles, Track_He, Track_Ho, Track_MateDistCD, Track_MateDistCDstd, nthfile,
                      logfHndl, \
                      Track_p1, Track_p2, Track_q1, Track_q2, preprocess_params["subpopemigration"], \
                      preprocess_params["subpopimmigration"], Track_FAvgMate, Track_MAvgMate, Track_FSDMate,
                      Track_MSDMate, \
                      SelectionDeathsEmi, SelectionDeathsImm, \
                      DisperseDeathsEmi, DisperseDeathsImm, \
                      Track_BreedEvents, params['gridformat'], \
                      MgSuccess, AdultNoMg, StrSuccess, \
                      Track_EggDeaths, Track_K, Track_N_Init_pop, N_Emigration_pop, N_EmiMortality, N_Immigration_pop,
                      N_ImmiMortality, Infected, Residors, Strayers1, Strayers2, Immigrators, PopSizes_Mean,
                      PopSizes_Std, AgeSizes_Mean, AgeSizes_Std, PackingDeathsEmi, PackingDeathsImm, Track_N_Init_age,
                      N_Emigration_age, N_Immigration_age, AgeDeathsOUT, AgeDeathsIN, PackingDeathsEmiAge,
                      PackingDeathsImmAge, Track_MatureCount, Track_ImmatureCount, Track_N_back_age, Track_N_out_age,
                      params['summaryOutput'], gen, Track_CaptureCount_Back, Track_CaptureCount_ClassBack,
                      Track_CaptureCount_Out, Track_CaptureCount_ClassOut, preprocess_params["age_size_mean"],
                      params['sizecontrol'],
                      ClassSizes_Mean,
                      ClassSizes_Std, Track_N_Init_class, SizeDeathsOUT, SizeDeathsIN, N_beforePack_Immi_pop,
                      N_beforePack_Immi_age, SelectionDeathsImm_Age0s, F_StrayDist, M_StrayDist, F_StrayDist_sd,
                      M_StrayDist_sd, F_ZtrayDist, M_ZtrayDist, F_ZtrayDist_sd, M_ZtrayDist_sd, F_HomeDist, M_HomeDist,
                      F_HomeDist_sd, M_HomeDist_sd, F_EmiDist, M_EmiDist, F_EmiDist_sd, M_EmiDist_sd,
                      params['SNPanswer'])

        # Print to log
        stringout = 'DoPostProcess(): '+str(datetime.datetime.now() -start_time1) + ''
        logMsg(logfHndl,stringout)
        print 'DoPostProcess(): ',str(datetime.datetime.now() -start_time1),''

        # Print to log
        stringout = 'End Monte Carlo Loop'+str(ithmcrun)+': '+str(datetime.datetime.now() -start_timeMC) + '\n'
        logMsg(logfHndl,stringout)
        print 'End Monte Carlo Loop',str(ithmcrun),': ',str(datetime.datetime.now() -start_timeMC),'\n'


    stringout = 'End Batch Loop' + str(ibatch) + ': ' + str(datetime.datetime.now() - start_timeB) + '\n'
    logMsg(logfHndl, stringout)
    print 'End Batch Loop', str(ibatch), ': ', str(datetime.datetime.now() - start_timeB), '\n'
    # End::Monte Carlo Loop
    #End::Batch Loop

def initPoolProc(outdir):
    sys.stdout = open(outdir + "/" + multiprocessing.current_process().name + ".out", "a", buffering=0)
    sys.stderr = open(outdir + "/" + multiprocessing.current_process().name + ".err", "a", buffering=0)

#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------
if __name__ == '__main__':

    # ------------------------------------------------------
    # Start timer, get script arguments, create log writeout
    # ------------------------------------------------------
    # Timing events: start
    start_time = datetime.datetime.now()
    foldertime = int(time.time())

    if len(sys.argv) >= 4:
        datadir = sys.argv[1]+'/'
        fileans = datadir+sys.argv[2]
        outdir = datadir+sys.argv[3]+str(foldertime)+'/'
        if len(sys.argv) == 5:
            noproc = int(sys.argv[4])
        else: # assume one processor
            noproc = 1

    # If user did not specify .rip file
    else:
        print "User must specify data directory, input file name, and output file directory (e.g., at command line type CDmetaPOP.py ../CDmetaPOP_data/ inputvariables.csv exampleout_foldername)."
        sys.exit(-1)

    # If .ip file does not exist
    if not os.path.exists(fileans):
        print("Cannot find or open runtime inputs file(%s)"%(fileans))
        sys.exit(-1)

    # If entered more processors than machine
    if noproc > multiprocessing.cpu_count():
        print('Warning: Specified more CPUs than on local machine. Using one less than '+str(multiprocessing.cpu_count()))

    # Create output file directory - will automatically put in the data directory
    os.mkdir(outdir)

    # This properly names log file
    logSessionPath = outdir+"CDmetaPOP.log"
    logfHndl =open(logSessionPath,'w')

    msgVerbose = True
    logMsg(logfHndl,"\n%s Release %s Version %s\n"%(appName,appRele,appVers))
    logMsg(logfHndl,"Author(s): %s"%(authorNames)+'\n')
    logMsg(logfHndl,"Session runtime inputs from: %s"%(fileans)+'\n\n')
    msgVerbose = False

    # ------------------------------------
    # Call DoUserInput()
    # ------------------------------------
    # Timing events: start
    start_time1 = datetime.datetime.now()

    param_list = loadFile(fileans)
    convertTypes(param_list)
    warnings, errors = checkErrors(param_list)

    if errors > 0:
        sys.exit(-1)

    # Print to log
    stringout = 'DoUserInput(): '+str(datetime.datetime.now() -start_time1) + ''
    logMsg(logfHndl,stringout)
    print 'DoUserInput(): ',str(datetime.datetime.now() -start_time1),''

    # ----------------------------------------
    # Begin Batch Looping - assign processors
    # ----------------------------------------
    # This loop is defined by the number of rows in inputvariables.csv
    pool = multiprocessing.Pool(processes=noproc, initializer=(initPoolProc if noproc > 1 else None), initargs=(outdir,))
    results = []

    for i, params in enumerate(param_list):
        results.append(pool.apply_async(doBatch, (i, params)))

    for result in results:
        result.wait()

# End::Main Loop
# Print to log
stringout = 'Total CDmetaPOP Simulation Time: '+str(datetime.datetime.now() -start_time) + ''
logMsg(logfHndl,stringout)
logfHndl.close()
print(stringout)
