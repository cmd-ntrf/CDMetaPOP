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
        subpopemigration = []
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
        subpopimmigration = []
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
        tupPreProcess = DoPreProcess(outdir, datadir, ibatch, ithmcrun, xyfilename, params['loci'], params['alleles'],
                                     params['cdevolveans'],
                                     params['cdinfect'], subpopemigration, subpopimmigration, params['sizecontrol'],
                                     params['eggFrequency'],
                                     params['mature_length_female'],
                                     params['mature_length_male'], params['mature_int_female'],
                                     params['mature_slope_female'],
                                     params['mature_int_male'],
                                     params['mature_slope_male'],
                                     params['startSelection'],
                                     params['correlation_matrix'],
                                     params['SNPanswer'])

        ithmcrundir = tupPreProcess[0]
        fitvals_pass = tupPreProcess[1]
        allelst = tupPreProcess[2]
        subpopemigration = tupPreProcess[3]
        subpopimmigration = tupPreProcess[4]
        age_percmort_out_mu = tupPreProcess[5]
        age_percmort_back_mu = tupPreProcess[6]
        age_Mg = tupPreProcess[7]
        age_S = tupPreProcess[8]
        age_mu = tupPreProcess[9]
        age_size_mean = tupPreProcess[10]
        age_size_std = tupPreProcess[11]
        xgridpop = tupPreProcess[12]
        ygridpop = tupPreProcess[13]
        SubpopIN_init = tupPreProcess[14]
        N0 = tupPreProcess[15]
        K_mu = tupPreProcess[16]
        dtype = tupPreProcess[17]
        outsizevals_pass = tupPreProcess[18]
        backsizevals_pass = tupPreProcess[19]
        popmort_out_pass = tupPreProcess[20]
        popmort_back_pass = tupPreProcess[21]
        Mg_pass = tupPreProcess[22]
        Str_pass = tupPreProcess[23]
        eggmort_pass = tupPreProcess[24]
        setmigrate = tupPreProcess[25]
        M_mature = tupPreProcess[26]
        F_mature = tupPreProcess[27]
        age_sigma = tupPreProcess[28]
        outgrowdays_pass = tupPreProcess[29]
        backgrowdays_pass = tupPreProcess[30]
        Kmu_pass = tupPreProcess[31]
        age_capture_out = tupPreProcess[32]
        age_capture_back = tupPreProcess[33]
        Kstd_pass = tupPreProcess[34]
        K_std = tupPreProcess[35]
        popmort_out_sd_pass = tupPreProcess[36]
        popmort_back_sd_pass = tupPreProcess[37]
        eggmort_sd_pass = tupPreProcess[38]
        outsizevals_sd_pass = tupPreProcess[39]
        backsizevals_sd_pass = tupPreProcess[40]
        outgrowdays_sd_pass = tupPreProcess[41]
        backgrowdays_sd_pass = tupPreProcess[42]
        size_percmort_out_mu = tupPreProcess[43]
        size_percmort_back_mu = tupPreProcess[44]
        age_percmort_out_sd = tupPreProcess[45]
        age_percmort_back_sd = tupPreProcess[46]
        size_percmort_out_sd = tupPreProcess[47]
        size_percmort_back_sd = tupPreProcess[48]
        pop_capture_back_pass = tupPreProcess[49]
        pop_capture_out_pass = tupPreProcess[50]
        pop_capture_back = tupPreProcess[51]
        natal = tupPreProcess[52]
        cor_mat = tupPreProcess[53]
        migrate = tupPreProcess[54]
        N0_pass = tupPreProcess[55]
        allefreqfiles_pass = tupPreProcess[56]
        classvarsfiles_pass = tupPreProcess[57]

        # Grab first one only
        K = K_mu # Initialize K with mu
        #pdb.set_trace() # check SubPopIN for XY chromosomes
        # Print to log
        stringout = 'DoPreProcess(): '+str(datetime.datetime.now() -start_time1) + ''
        logMsg(logfHndl,stringout)
        print(stringout)

        # ---------------------------------
        # Call GetMetrics()
        # ---------------------------------

        # Timing events: start
        start_time1 = datetime.datetime.now()

        GetMetrics(SubpopIN_init, K, Track_N_Init_pop, Track_K, params['loci'], params['alleles'], 0, Track_Ho,
                   Track_Alleles,
                   Track_He,
                   Track_p1, Track_p2, Track_q1, Track_q2, Infected, Residors, Strayers1, Strayers2, Immigrators,
                   PopSizes_Mean, PopSizes_Std, AgeSizes_Mean, AgeSizes_Std, Track_ToTMales, Track_ToTFemales,
                   Track_BreedMales, Track_BreedFemales, Track_N_Init_age, Track_MatureCount, Track_ImmatureCount,
                   params['sizecontrol'], age_size_mean, ClassSizes_Mean, ClassSizes_Std, Track_N_Init_class,
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

        DoUpdate(SubpopIN_init, K, xgridpop, ygridpop, -1, nthfile, ithmcrundir, params['loci'], params['alleles'],
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
                SubpopIN = SubpopIN_init
                del SubpopIN_init
                # Use NatalPop
                sourcePop = 'NatalPop'
            else:
                sourcePop = 'ImmiPop'

            # Exit the system if population is 0 or 1
            if sum(N0) <= 1:
                stringout = 'Population went extinct, program ended.'
                logMsg(logfHndl,stringout)
                print('Population went extinct after generation '+str(gen-1)+'.\n')
                break

            # ---------------------------------
            # Call CDClimate()
            # ---------------------------------

            # Timing events: start
            start_time1 = datetime.datetime.now()

            # Check gen time equal to cdclimgentime
            for icdtime in xrange(len(cdclimgentime)):
                if gen == int(cdclimgentime[icdtime]):
                    tupClimate = DoCDClimate(datadir, icdtime, cdclimgentime, params['mate_cdmat'],
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
                                             Mg_pass,
                                             Str_pass, Kmu_pass, outsizevals_pass, backsizevals_pass, outgrowdays_pass,
                                             backgrowdays_pass, fitvals_pass, popmort_back_pass, popmort_out_pass,
                                             eggmort_pass, Kstd_pass, popmort_back_sd_pass, popmort_out_sd_pass,
                                             eggmort_sd_pass, outsizevals_sd_pass, backsizevals_sd_pass,
                                             outgrowdays_sd_pass, backgrowdays_sd_pass, pop_capture_back_pass,
                                             pop_capture_out_pass, params['cdevolveans'], N0_pass, allefreqfiles_pass,
                                             classvarsfiles_pass)

                    cdmatrix_mate = tupClimate[0]
                    cdmatrix_FOut = tupClimate[1]
                    cdmatrix_MOut = tupClimate[2]
                    cdmatrix_FBack = tupClimate[3]
                    cdmatrix_MBack = tupClimate[4]
                    cdmatrix_StrBack = tupClimate[5]
                    thresh_mate = tupClimate[6]
                    thresh_FOut = tupClimate[7]
                    thresh_MOut = tupClimate[8]
                    thresh_FBack = tupClimate[9]
                    thresh_MBack = tupClimate[10]
                    thresh_Str = tupClimate[11]
                    Mg = tupClimate[12]
                    Str = tupClimate[13]
                    Str_ScaleMin = tupClimate[14]
                    Str_ScaleMax = tupClimate[15]
                    FdispBack_ScaleMin = tupClimate[16]
                    FdispBack_ScaleMax = tupClimate[17]
                    MdispBack_ScaleMin = tupClimate[18]
                    MdispBack_ScaleMax = tupClimate[19]
                    FdispOut_ScaleMin = tupClimate[20]
                    FdispOut_ScaleMax = tupClimate[21]
                    MdispOut_ScaleMin = tupClimate[22]
                    MdispOut_ScaleMax = tupClimate[23]
                    mate_ScaleMin = tupClimate[24]
                    mate_ScaleMax = tupClimate[25]
                    outsizevals_mu = tupClimate[26]
                    backsizevals_mu = tupClimate[27]
                    outgrowdays_mu = tupClimate[28]
                    backgrowdays_mu = tupClimate[29]
                    fitvals = tupClimate[30]
                    K_mu = tupClimate[31]
                    popmort_back_mu = tupClimate[32]
                    popmort_out_mu = tupClimate[33]
                    eggmort_mu = tupClimate[34]
                    K_std = tupClimate[35]
                    popmort_back_sd = tupClimate[36]
                    popmort_out_sd = tupClimate[37]
                    eggmort_sd = tupClimate[38]
                    outsizevals_sd = tupClimate[39]
                    backsizevals_sd = tupClimate[40]
                    outgrowdays_sd = tupClimate[41]
                    backgrowdays_sd = tupClimate[42]
                    pop_capture_back = tupClimate[43]
                    pop_capture_out = tupClimate[44]
                    mateno = tupClimate[45]
                    FdispOutno = tupClimate[46]
                    MdispOutno = tupClimate[47]
                    FdispBackno = tupClimate[48]
                    MdispBackno = tupClimate[49]
                    Strno = tupClimate[50]
                    tempN0 = tupClimate[51]
                    tempAllelefile = tupClimate[52]
                    tempClassVarsfile = tupClimate[53]

                    # ----------------------------------------
                    # Introduce new individuals
                    # ----------------------------------------
                    if (gen != 0 and len(N0_pass[0].split('|')) > 1):
                        SubpopIN = AddIndividuals(SubpopIN, tempN0, tempAllelefile, tempClassVarsfile, datadir,
                                                  params['loci'],
                                                  params['alleles'], params['sizecontrol'], params['cdinfect'],
                                                  params['SNPanswer'],
                                                  params['cdevolveans'], params['startSelection'],
                                                  fitvals,
                                                  params['eggFrequency'], params['mature_length_female'],
                                                  params['mature_length_male'],
                                                  params['mature_int_female'], params['mature_slope_female'],
                                                  params['mature_int_male'], params['mature_slope_male'], dtype, N0,
                                                  natal, gen)
            # -------------------------------------------
            # Update stochastic parameters each year here
            # -------------------------------------------
            tupStoch = DoStochasticUpdate(K_mu, K_std, popmort_back_mu, popmort_back_sd, popmort_out_mu, popmort_out_sd,
                                          eggmort_mu, eggmort_sd, outsizevals_mu, outsizevals_sd, backsizevals_mu,
                                          backsizevals_sd, outgrowdays_mu, outgrowdays_sd, backgrowdays_mu,
                                          backgrowdays_sd, age_percmort_out_mu, age_percmort_out_sd,
                                          age_percmort_back_mu, age_percmort_back_sd, size_percmort_out_mu,
                                          size_percmort_out_sd, size_percmort_back_mu, size_percmort_back_sd,
                                          params['Egg_Mortality'], params['Egg_Mortality_StDev'], cor_mat)
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
                               params['Freplace'], params['Mreplace'], mateno, thresh_mate, \
                               cdmatrix_mate, Track_MateDistCD, xgridpop, \
                               ygridpop, Track_MateDistCDstd, Track_FAvgMate, Track_MAvgMate, Track_FSDMate,
                               Track_MSDMate, Track_BreedEvents, gen, sourcePop, dtype, mate_ScaleMax, mate_ScaleMin,
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
                                                 age_mu, age_sigma, params['sizecontrol'], \
                                                 params['Egg_Mean_par1'], params['Egg_Mean_par2'],
                                                 params['Egg_Mean_ans'],
                                                 params['equalClutchSize'], dtype,
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

            SubpopIN = DoUpdate(SubpopIN, K, xgridpop, ygridpop, gen, nthfile, ithmcrundir, params['loci'],
                                params['alleles'],
                                logfHndl,
                                'Middle', params['growth_option'], params['cdevolveans'], fitvals,
                                params['startSelection'],
                                age_capture_back,
                                pop_capture_back,
                                Track_CaptureCount_Back, Track_CaptureCount_ClassBack, params['sizecontrol'],
                                age_size_mean,
                                Track_N_back_age,
                                params['eggFrequency'], backsizevals, params['growth_Loo'], params['growth_R0'],
                                params['growth_temp_max'],
                                params['growth_temp_CV'],
                                params['growth_temp_t0'], backgrowdays, sourcePop, params['sizecontrol'], M_mature,
                                F_mature,
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

            SubpopIN = DoEmigration(SubpopIN, K, FdispOutno, \
                                    MdispOutno, cdmatrix_FOut, cdmatrix_MOut, gen, xgridpop, ygridpop, F_EmiDist,
                                    M_EmiDist, params['cdevolveans'], fitvals, F_EmiDist_sd, M_EmiDist_sd,
                                    subpopemigration, \
                                    SelectionDeathsEmi, DisperseDeathsEmi, \
                                    params['startSelection'], Mg, \
                                    MgSuccess, AdultNoMg, sum(params['alleles']), age_Mg, thresh_FOut, thresh_MOut,
                                    N_Emigration_pop, sourcePop, dtype, setmigrate, params['sizecontrol'],
                                    age_size_mean,
                                    PackingDeathsEmi, N_Emigration_age, params['loci'], params['muterate'],
                                    params['mtdna'],
                                    params['mutationtype'],
                                    FdispOut_ScaleMax, FdispOut_ScaleMin, MdispOut_ScaleMax, MdispOut_ScaleMin,
                                    params['FdispmoveOutparA'], params['FdispmoveOutparB'], params['FdispmoveOutparC'],
                                    params['MdispmoveOutparA'], params['MdispmoveOutparB'], params['MdispmoveOutparC'],
                                    params['popmodel'],
                                    PackingDeathsEmiAge,
                                    ithmcrundir, params['popmodel_par1'], params['implementSelection'],
                                    age_percmort_out, migrate)

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
                                   gen, N_EmiMortality, AgeDeathsOUT, params['sizecontrol'], age_size_mean,
                                   size_percmort_out, SizeDeathsOUT, params['constMortans'], params['popmodel'])

            # Print to log
            stringout = 'DoOutMortality(): '+str(datetime.datetime.now() -start_time1) + ''
            logMsg(logfHndl,stringout)
            print 'DoOutMortality(): ',str(datetime.datetime.now() -start_time1),''

            # ----------------------------------------------------
            # Call DoUpdate() - grow, capture, and optional output indSample.csv
            # ----------------------------------------------------
            start_time1 = datetime.datetime.now()  # Timing events: start

            SubpopIN = DoUpdate(SubpopIN, K, xgridpop, ygridpop, gen, nthfile, ithmcrundir, params['loci'],
                                params['alleles'],
                                logfHndl,
                                params['gridsampling'], params['growth_option'], 'N', [], params['startSelection'],
                                age_capture_out,
                                pop_capture_out,
                                Track_CaptureCount_Out, Track_CaptureCount_ClassOut, params['sizecontrol'],
                                age_size_mean,
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

            SubpopIN = DoImmigration(SubpopIN, K, N0, natal, FdispBackno, \
                                     MdispBackno, cdmatrix_FBack, cdmatrix_MBack, gen, \
                                     xgridpop, ygridpop, params['cdevolveans'], fitvals, subpopimmigration, \
                                     SelectionDeathsImm, DisperseDeathsImm, params['startSelection'], Str, \
                                     StrSuccess, \
                                     Strno, cdmatrix_StrBack, age_S, thresh_FBack, thresh_MBack, thresh_Str,
                                     N_Immigration_pop, dtype, params['sizecontrol'], age_size_mean, PackingDeathsImm,
                                     N_Immigration_age, FdispBack_ScaleMax, FdispBack_ScaleMin, MdispBack_ScaleMax,
                                     MdispBack_ScaleMin,
                                     params['FdispmoveBackparA'], params['FdispmoveBackparB'],
                                     params['FdispmoveBackparC'],
                                     params['MdispmoveBackparA'], params['MdispmoveBackparB'],
                                     params['MdispmoveBackparC'],
                                     Str_ScaleMax,
                                     Str_ScaleMin, params['StrayBackparA'], params['StrayBackparB'],
                                     params['StrayBackparC'],
                                     params['popmodel'],
                                     PackingDeathsImmAge,
                                     ithmcrundir, params['popmodel_par1'], noOffspring, Bearpairs, age_size_std,
                                     params['Egg_FemalePercent'],
                                     sourcePop, params['transmissionprob'], M_mature, F_mature,
                                     params['mature_slope_male'], params['mature_int_male'],
                                     params['mature_slope_female'],
                                     params['mature_int_female'], params['mature_length_male'],
                                     params['mature_length_female'],
                                     params['loci'],
                                     params['muterate'], params['mtdna'],
                                     params['mutationtype'],
                                     params['startGenes'],
                                     allelst, params['HomeAttempt'], params['implementSelection'],
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
                                   gen, N_ImmiMortality, AgeDeathsIN, params['sizecontrol'], age_size_mean,
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
                       Track_ImmatureCount, params['sizecontrol'], age_size_mean, ClassSizes_Mean, ClassSizes_Std, Track_N_Init_class,
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

        DoPostProcess(ithmcrundir, params['gendmatans'], params['loci'], params['alleles'], params['runtime'], \
                      Track_ToTFemales, Track_ToTMales, Track_BreedFemales, Track_BreedMales, Track_Births, PopDeathsIN, \
                      PopDeathsOUT, Track_Alleles, Track_He, Track_Ho, Track_MateDistCD, Track_MateDistCDstd, nthfile,
                      logfHndl, \
                      Track_p1, Track_p2, Track_q1, Track_q2, subpopemigration, \
                      subpopimmigration, Track_FAvgMate, Track_MAvgMate, Track_FSDMate, Track_MSDMate, \
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
                      Track_CaptureCount_Out, Track_CaptureCount_ClassOut, age_size_mean, params['sizecontrol'],
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
