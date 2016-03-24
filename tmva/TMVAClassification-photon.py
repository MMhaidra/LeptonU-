#!/usr/bin/env python
# @(#)root/tmva $Id$#
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser

# --------------------------------------------

# Default settings for command line arguments
DEFAULT_OUTFNAME = "TMVA.root"
#DEFAULT_INFNAMEBKG  = "/exp/LHCb/amhis/LeptonU-/tuples/data/real-data-electrons-jpsi-tight-23022016.root"
DEFAULT_INFNAMEBKG = "/exp/LHCb/amhis/LeptonU-/tuples/data/LeptonU-total-electrons-23022016.root"
#DEFAULT_INFNAMESIG  = "/exp/LHCb/amhis/LeptonU-/tuples/montecarlo/spring16/mc-15124011-leptonU.root"
DEFAULT_INFNAMESIG = "/exp/LHCb/amhis/LeptonU-/tuples/montecarlo/spring16/signalonly-mc-15124011-leptonU.root"
#DEFAULT_TREESIG  = "Tuple_Bu2LLK_eeLine2/DecayTree"
DEFAULT_TREESIG  = "DecayTree"
#DEFAULT_TREEBKG  = "DecayTree"
DEFAULT_TREEBKG  = "TupleFromData_Bu2LLK_eeLine2/DecayTree"
#DEFAULT_TREEBKG  = "DecayTree"
 
DEFAULT_METHODS  =  "BDTG"#,Cuts,CuBtsD,CutsPCA,CutsGA,CutsSA,Likelihood,LikelihoodD,LikelihoodPCA,LikelihoodKDE,LikelihoodMIX,PDERS,PDERSD,PDERSPCA,PDEFoam,PDEFoamBoost,KNN,LD,Fisher,FisherG,BoostedFisher,HMatrix,FDA_GA,FDA_SA,FDA_MC,FDA_MT,FDA_GAMT,FDA_MCMT,MLP,MLPBFGS,MLPBNN,CFMlpANN,TMlpANN,SVM,BDT,BDTD,BDTG,BDTB,RuleFit"

# Print usage help
def usage():
    print " "
    print "Usage: python %s [options]" % sys.argv[0]
    print "  -m | --methods    : gives methods to be run (default: all methods)"
    print "  -is | --inputfile signal  : name of input ROOT file (default: '%s')" %  DEFAULT_INFNAMESIG
    print "  -ib | --inputfile background : name of input ROOT file (default: '%s')" % DEFAULT_INFNAMEBKG
 
    print "  -o | --outputfile : name of output ROOT file containing results (default: '%s')" % DEFAULT_OUTFNAME
    print "  -t | --inputtrees signal : input ROOT Trees for signal and background (default: '%s %s')" % DEFAULT_TREESIG
   
    print "  -t | --inputtrees : input ROOT Trees for signal and background (default: '%s %s')" % DEFAULT_TREEBKG
    print "  -v | --verbose"
    print "  -? | --usage      : print this help message"
    print "  -h | --help       : print this help message"
    print " "

# Main routine
def main():

    try:
        # retrive command line options
        shortopts  = "m:is:ib:ts:tb:o:vh?"
        longopts   = ["methods=", "inputfile_signal=","inputfile_background=" , "inputtree_signal=", "inputtree_background=",  "outputfile=", "verbose", "help", "usage"]
        opts, args = getopt.getopt( sys.argv[1:], shortopts, longopts )

    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        usage()
        sys.exit(1)

    infnameSig     = DEFAULT_INFNAMESIG
    infnameBkg     = DEFAULT_INFNAMEBKG
    treeNameSig = DEFAULT_TREESIG
    treeNameBkg = DEFAULT_TREEBKG
    outfname    = DEFAULT_OUTFNAME
    methods     = DEFAULT_METHODS
    verbose     = False
    for o, a in opts:
        if o in ("-?", "-h", "--help", "--usage"):
            usage()
            sys.exit(0)
        elif o in ("-m", "--methods"):
            methods = a
        elif o in ("-is", "--inputfile_signal"):
            infnameSig = a
        elif o in ("-ib", "--inputfile_background"):
            infnameBkg = a
        elif o in ("-o", "--outputfile"):
            outfname = a
        elif o in ("-ts", "--inputree_signal"):
            treeNameSig = a 
        elif o in ("-tb", "--inputtree_background"):
            treeNameBkg = a
        elif o in ("-v", "--verbose"):
            verbose = True

    # Print methods
    mlist = methods.replace(' ',',').split(',')
    print "=== TMVAClassification: use method(s)..."
    for m in mlist:
        if m.strip() != '':
            print "=== - <%s>" % m.strip()

    # Import ROOT classes
    from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut
    
    # check ROOT version, give alarm if 5.18 
    if gROOT.GetVersionCode() >= 332288 and gROOT.GetVersionCode() < 332544:
        print "*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA"
        print "*** does not run properly (function calls with enums in the argument are ignored)."
        print "*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),"
        print "*** or use another ROOT version (e.g., ROOT 5.19)."
        sys.exit(1)
    
    # Logon not automatically loaded through PyROOT (logon loads TMVA library) load also GUI
    gROOT.SetMacroPath( "./" )
    gROOT.Macro       ( "./TMVAlogon.C" )    
    gROOT.LoadMacro   ( "./TMVAGui.C" )
    
    # Import TMVA classes from ROOT
    from ROOT import TMVA

    # Output file
    outputFile = TFile( outfname, 'RECREATE' )
    
    # Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for more factory options)
    # All TMVA output can be suppressed by removing the "!" (not) in 
    # front of the "Silent" argument in the option string
    factory = TMVA.Factory( "TMVAClassification", outputFile, 
                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

    # Set verbosity
    factory.SetVerbose( verbose )
    
    # If you wish to modify default settings 
    # (please check "src/Config.h" to see all available global options)
    #    gConfig().GetVariablePlotting()).fTimesRMS = 8.0
    #    gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory"

    # Define the input variables that shall be used for the classifier training
    # note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    # [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    #factory.AddVariable( "myvar1 := var1+var2", 'F' )
    #factory.AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' )
    #factory.AddVariable( "var3",                "Variable 3", "units", 'F' )
    #factory.AddVariable( "var4",                "Variable 4", "units", 'F' )
    #factory.AddVariable ("Polarity", "", "" , 'F')
    #factory.AddVariable ("GpsTime", "", "" , 'F')

    #Lb variables
    factory.AddVariable( "Lambdab_PT", "Lambdab_PT", "" , 'F' )
    factory.AddVariable( "Lambdab_DIRA_OWNPV",                "Lambdab_DIRAOWNPV", "", 'F' )
    factory.AddVariable( "Lambdab_FDCHI2_OWNPV",                "Lambdab_FDCHI2_OWNPV", "", 'F' )
    factory.AddVariable( "Lambdab_LOKI_DTF_CHI2NDOF",                "Lambdab_LOKI_DTF_CHI2NDOF", "", 'F' )
    factory.AddVariable( "Lambdab_ENDVERTEX_CHI2",                "Lambdab_ENDVERTEX_CHI2", "", 'F' )




    factory.AddVariable( "Jpsi_PT",   "Jpsi_PT", "" , 'F' )
    factory.AddVariable( "Jpsi_DIRA_OWNPV",                "Jpsi_DIRAOWNPV", "", 'F' )
    factory.AddVariable( "Jpsi_FDCHI2_OWNPV",                "Jpsi_FDCHI2_OWNPV", "", 'F' )
    factory.AddVariable( "Jpsi_LOKI_DTF_CHI2NDOF",                "Jpsi_LOKI_DTF_CHI2NDOF", "", 'F' )
    factory.AddVariable( "Jpsi_ENDVERTEX_CHI2",                "Jpsi_ENDVERTEX_CHI2", "", 'F' )

    
    factory.AddVariable( "Lambdastar_PT",   "Lambdastar", "" ,  'F' )
    factory.AddVariable( "Lambdastar_DIRA_OWNPV",                "Lambdastar_DIRAOWNPV", "", 'F' )
    factory.AddVariable( "Lambdastar_FDCHI2_OWNPV",                "Lambdastar_FDCHI2_OWNPV", "", 'F' )
    factory.AddVariable( "Lambdastar_LOKI_DTF_CHI2NDOF",                "Lambdastar_LOKI_DTF_CHI2NDOF", "", 'F' )
    factory.AddVariable( "Lambdastar_ENDVERTEX_CHI2",                "Lambdastar_ENDVERTEX_CHI2", "", 'F' )
        





    factory.AddVariable( "minHadron_PT:=min(Kaon_PT, Proton_PT)  ",   "minHadron_PT", "" , 'F' )
    factory.AddVariable( "maxHadron_PT:=max(Kaon_PT, Proton_PT)  ",   "maxHadron_PT", "" , 'F' )
    factory.AddVariable( "minHadron_IPCHI2:=min(Kaon_IPCHI2_OWNPV, Proton_IPCHI2_OWNPV)  ",   "minHadron_PTCHI2", "" , 'F' )
    factory.AddVariable( "maxHadron_IPCHI2:=max(Kaon_IPCHI2_OWNPV, Proton_IPCHI2_OWNPV)  ",   "maxHadron_IPCHI2", "" , 'F' )



     # You can add so-called "Spectator variables", which are not used in the MVA training, 
    # but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the 
    # input variables, the response values of all trained MVAs, and the spectator variables
    factory.AddSpectator( "Lambdab_MM",  "Lambdab_MM", "", 2000,8000 )
    factory.AddSpectator( "Jpsi_MM",  "Jpsi_MM", "", 2000,8000)

    # Read input data
    #if gSystem.AccessPathName( infname ) != 0: gSystem.Exec( "wget http://root.cern.ch/files/" + infname )
        
    #input = TFile.Open( infname )
    input_signal = TFile.Open( infnameSig )
    input_background = TFile.Open( infnameBkg )

    # Get the signal and background trees for training
    signal      = input_signal.Get( treeNameSig )
    background  = input_background.Get( treeNameBkg )
   


    #for hell in signal : 
        #print signal.Lambdab_ETA 


 
    # Global event weights (see below for setting event-wise weights)
    signalWeight     = 1.0
    backgroundWeight = 1.0

    print 'is this modern hell ? '

    print input_signal 
    print input_background
    print signal 
    print background 
  # ====== register trees ====================================================
    #
    # the following method is the prefered one:
    # you can add an arbitrary number of signal or background trees
    factory.AddSignalTree    ( signal,     signalWeight     )
    factory.AddBackgroundTree( background, backgroundWeight )

    # To give different trees for training and testing, do as follows:
    #    factory.AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" )
    #    factory.AddSignalTree( signalTestTree,     signalTestWeight,  "Test" )
    
    # Use the following code instead of the above two or four lines to add signal and background 
    # training and test events "by hand"
    # NOTE that in this case one should not give expressions (such as "var1+var2") in the input 
    #      variable definition, but simply compute the expression before adding the event
    #
    #    # --- begin ----------------------------------------------------------
    #    
    # ... *** please lookup code in TMVA/macros/TMVAClassification.C ***
    #    
    #    # --- end ------------------------------------------------------------
    #
    # ====== end of register trees ==============================================    
            
    # Set individual event weights (the variables must exist in the original TTree)
    #    for signal    : factory.SetSignalWeightExpression    ("weight1*weight2");
    #    for background: factory.SetBackgroundWeightExpression("weight1*weight2");
   # factory.SetBackgroundWeightExpression( "weight" )

    # Apply additional cuts on the signal and background sample. 
    # example for cut: mycut = TCut( "abs(var1)<0.5 && abs(var2-0.5)<1" )  
   #    from ROOT import * 
    #mycutSig = TCut( TString("Lambdab_BKGCAT==10") ) #signal MC 
    #mycutSig = TCut("Lambdab_BKGCAT<15" ) #signal MC 
   
 #   print mycutSig
    mycutSig = TCut( "" ) #signal MC 
        
    mycutBkg = TCut( "Lambdab_MM > 5700" ) #upper sideband
    
    # Here, the relevant variables are copied over in new, slim trees that are
    # used for TMVA training and testing
    # "SplitMode=Random" means that the input events are randomly shuffled before
    # splitting them into training and test samples
    factory.PrepareTrainingAndTestTree( mycutSig, mycutBkg,
                                        "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" )

    # --------------------------------------------------------------------------------------------------

    # ---- Book MVA methods
    #
    # please lookup the various method configuration options in the corresponding cxx files, eg:
    # src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    # it is possible to preset ranges in the option string in which the cut optimisation should be done:
    # "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

    # Cut optimisation
    if "Cuts" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "Cuts",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" )

    if "CutsD" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "CutsD",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" )

    if "CutsPCA" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "CutsPCA",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" )

    if "CutsGA" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "CutsGA",
                            "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" )

    if "CutsSA" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "CutsSA",
                            "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" )

    # Likelihood ("naive Bayes estimator")
    if "Likelihood" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "Likelihood",
                            "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )

    # Decorrelated likelihood
    if "LikelihoodD" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "LikelihoodD",
                            "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" )

    # PCA-transformed likelihood
    if "LikelihoodPCA" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "LikelihoodPCA",
                            "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ) 

    # Use a kernel density estimator to approximate the PDFs
    if "LikelihoodKDE" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "LikelihoodKDE",
                            "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ) 

    # Use a variable-dependent mix of splines and kernel density estimator
    if "LikelihoodMIX" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "LikelihoodMIX",
                            "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ) 

    # Test the multi-dimensional probability density estimator
    # here are the options strings for the MinMax and RMS methods, respectively:
    #      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
    #      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
    if "PDERS" in mlist:
        factory.BookMethod( TMVA.Types.kPDERS, "PDERS",
                            "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" )

    if "PDERSD" in mlist:
        factory.BookMethod( TMVA.Types.kPDERS, "PDERSD",
                            "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" )

    if "PDERSPCA" in mlist:
        factory.BookMethod( TMVA.Types.kPDERS, "PDERSPCA",
                             "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" )

   # Multi-dimensional likelihood estimator using self-adapting phase-space binning
    if "PDEFoam" in mlist:
        factory.BookMethod( TMVA.Types.kPDEFoam, "PDEFoam",
                            "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" )

    if "PDEFoamBoost" in mlist:
        factory.BookMethod( TMVA.Types.kPDEFoam, "PDEFoamBoost",
                            "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" )

    # K-Nearest Neighbour classifier (KNN)
    if "KNN" in mlist:
        factory.BookMethod( TMVA.Types.kKNN, "KNN",
                            "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" )

    # H-Matrix (chi2-squared) method
    if "HMatrix" in mlist:
        factory.BookMethod( TMVA.Types.kHMatrix, "HMatrix", "!H:!V" )

    # Linear discriminant (same as Fisher discriminant)
    if "LD" in mlist:
        factory.BookMethod( TMVA.Types.kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

    # Fisher discriminant (same as LD)
    if "Fisher" in mlist:
        factory.BookMethod( TMVA.Types.kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

    # Fisher with Gauss-transformed input variables
    if "FisherG" in mlist:
        factory.BookMethod( TMVA.Types.kFisher, "FisherG", "H:!V:VarTransform=Gauss" )

    # Composite classifier: ensemble (tree) of boosted Fisher classifiers
    if "BoostedFisher" in mlist:
        factory.BookMethod( TMVA.Types.kFisher, "BoostedFisher", 
                            "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" )

    # Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    if "FDA_MC" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_MC",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

    if "FDA_GA" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_GA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

    if "FDA_SA" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_SA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

    if "FDA_MT" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_MT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

    if "FDA_GAMT" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_GAMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

    if "FDA_MCMT" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_MCMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

    # TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
    if "MLP" in mlist:
        factory.BookMethod( TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" )

    if "MLPBFGS" in mlist:
        factory.BookMethod( TMVA.Types.kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" )

    if "MLPBNN" in mlist:
        factory.BookMethod( TMVA.Types.kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ) # BFGS training with bayesian regulators

    # CF(Clermont-Ferrand)ANN
    if "CFMlpANN" in mlist:
        factory.BookMethod( TMVA.Types.kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ) # n_cycles:#nodes:#nodes:...  

    # Tmlp(Root)ANN
    if "TMlpANN" in mlist:
        factory.BookMethod( TMVA.Types.kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ) # n_cycles:#nodes:#nodes:...

    # Support Vector Machine
    if "SVM" in mlist:
        factory.BookMethod( TMVA.Types.kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" )

    # Boosted Decision Trees
    if "BDTG" in mlist:
        factory.BookMethod( TMVA.Types.kBDT, "BDTG",
                            "!H:!V:NTrees=1000:MinNodeSize=1.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" )                        

    if "BDT" in mlist:
        factory.BookMethod( TMVA.Types.kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" )

    if "BDTB" in mlist:
        factory.BookMethod( TMVA.Types.kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" )

    if "BDTD" in mlist:
        factory.BookMethod( TMVA.Types.kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" )

    # RuleFit -- TMVA implementation of Friedman's method
    if "RuleFit" in mlist:
        factory.BookMethod( TMVA.Types.kRuleFit, "RuleFit",
                            "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" )

    # --------------------------------------------------------------------------------------------------
            
    # ---- Now you can tell the factory to train, test, and evaluate the MVAs. 

    # Train MVAs
    factory.TrainAllMethods()
    
    # Test MVAs
    factory.TestAllMethods()
    
    # Evaluate MVAs
    factory.EvaluateAllMethods()    
    
    # Save the output.
    outputFile.Close()
    
    print "=== wrote root file %s\n" % outfname
    print "=== TMVAClassification is done!\n"
    
    # open the GUI for the result macros    
    gROOT.ProcessLine( "TMVAGui(\"%s\")" % outfname )
    
    # keep the ROOT thread running
    gApplication.Run() 

# ----------------------------------------------------------

if __name__ == "__main__":
    main()
