
import os, sys, json
import numpy as np
import ROOT
import matplotlib.pyplot as plt
import collections
try:
    import matplotlib.patheffects as PathEffects
    PlotMatrices = True
except:
    print "Could not import PathEffects. Correlation matrices will not be made."
    PlotMatrices = False

plt.rcParams['text.usetex'] = True

from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
from matplotlib import cm

from runToleranceAnalysis import LoadData, InitialiseFittingFunction, CreateFluxMultiGraph, CreateFluxGraphArray, CreateDeviationGraphs

def PlotSaveFluxMultiGraph(fgs,leg,U1,U2,D1,D2,primary,outdir,emin,emax,fit_func,fname,extension='pdf'):

    canvas = ROOT.TCanvas("c2")
    # Draw the uncertainties first so they end up behind the data
    U2.Draw("AF")
    # U2.SetTitle("Experimental Cosmic Ray %s Flux"%primary)
    U2.GetYaxis().SetTitle('Flux [ (GeV/n m^{2} s sr)^{-1} ]')
    xmin = '1E%i'%(int(np.log10(emin)))
    xmax = '1E%i'%(int(np.log10(emax))+1)
    U2.GetXaxis().SetLimits(float(xmin),float(xmax))
    U2.GetXaxis().SetLabelOffset(0.001)
    U2.GetXaxis().SetTitleOffset(1.2)
    ymin = '1E%i'%(int(np.log10(fit_func.Eval(emax)))-1)
    ymax = '1E%i'%(int(np.log10(fit_func.Eval(emin)))+1)
    U2.GetYaxis().SetRangeUser(float(ymin),float(ymax))
    U2.GetYaxis().SetTitleOffset(1.1)
    # H3a is parameterised as total cosmic ray energy
    if fname == 'H3a':
        U2.GetXaxis().SetTitle('Primary Energy [GeV]')
    # GSHL and power law are in terms of energy per nucleon
    else:
        U2.GetXaxis().SetTitle('Primary Energy [GeV/n]')
    D2.Draw("Fsame")
    U1.Draw("Fsame")
    D1.Draw("Fsame")

    fit_func.SetLineWidth(1)
    fit_func.SetLineColor(1)
    fit_func.Draw("Lsame")

    # Now draw the data
    # This isn't done as a multigraph because ROOT
    for fg in fgs:
        fg.Draw("PEsame")
    canvas.SetLogy()
    canvas.SetLogx()
    leg.Draw()

    # Redraw the axes so the ticks aren't covered by the contours
    ROOT.gPad.RedrawAxis()

    if outdir == os.getcwd():
        canvas.SaveAs("%s/Global%sFit.%s"%(outdir,primary,extension))
    else:
        canvas.SaveAs("%s/FinalPlots/Global%sFit.%s"%(outdir,primary,extension))
    canvas.Close()

def CalculateSigmaF(fname,grad,covMat,energy,primary,Barr=False):

    TotalTerm = 0.0
    for i in range(0,len(grad)):
        for j in range(0,len(grad)):
            Term = grad[i]*grad[j] * covMat[i][j]
            if fname == 'GSHL':
                if Barr and energy > 200.0:
                    # Inflate proton uncertainties on d by factor 3
                    if primary == 'Proton':
                        if i == 3:
                            Term *= 3.0
                        if j == 3:
                            Term *= 3.0
                    # Inflate all other uncertainties on d by factor 2
                    else:
                        if i == 3:
                            Term *= 2.0
                        if j == 3:
                            Term *= 2.0
            TotalTerm += Term

    sigmaF = np.sqrt(TotalTerm)
    return sigmaF

def GenerateFitUncertaintyLines(emin,emax,fit_func,fname,primary,covMat,Barr=False):

    realup = ROOT.TGraph()
    realdown = ROOT.TGraph()
    realup2 = ROOT.TGraph()
    realdown2 = ROOT.TGraph()

    realup.SetPoint(0,emin,fit_func.Eval(emin))
    realdown.SetPoint(0,emin,fit_func.Eval(emin))
    realup2.SetPoint(0,emin,fit_func.Eval(emin))
    realdown2.SetPoint(0,emin,fit_func.Eval(emin))

    nPoint = 200
    step = (np.log10(emax)-np.log10(emin))/nPoint

    for i in range(0,nPoint+1):
        x = np.array([np.power(10,np.log10(emin)+i*step)])
        grad = np.empty(fit_func.GetNpar())
        fit_func.GradientPar(x,grad,0.001)
        sigmaF = CalculateSigmaF(fname,grad,covMat,x[0],primary,Barr)
        realup.SetPoint(i+1,x[0],sigmaF+fit_func.Eval(x[0]))
        realdown.SetPoint(i+1,x[0],-sigmaF+fit_func.Eval(x[0]))
        realup2.SetPoint(i+1,x[0],2*sigmaF+fit_func.Eval(x[0]))
        realdown2.SetPoint(i+1,x[0],-2*sigmaF+fit_func.Eval(x[0]))

    for i in range(0,nPoint+1):
        en = np.power(10,np.log10(emin)+(nPoint-i)*step)
        realup.SetPoint(realup.GetN(),en,fit_func.Eval(en))
        realdown.SetPoint(realdown.GetN(),en,fit_func.Eval(en))
        realup2.SetPoint(realup2.GetN(),en,fit_func.Eval(en))
        realdown2.SetPoint(realdown2.GetN(),en,fit_func.Eval(en))

    realup.SetFillColor(390);
    realup.SetLineColor(390);
    realdown.SetFillColor(390);
    realdown.SetLineColor(390);

    realup2.SetFillColor(406);
    realup2.SetLineColor(406);
    realdown2.SetFillColor(406);
    realdown2.SetLineColor(406);

    return realup, realdown, realup2, realdown2

def GenerateDeviationUncertaintyLines(emin,emax,fit_func,fname,outdir,primary,covMat,Barr=False):

    if outdir == os.getcwd():
        outfile1 = open('%s/%sOneSigmaContour.dat'%(outdir,fname),'w')
        outfile2 = open('%s/%sTwoSigmaContour.dat'%(outdir,fname),'w')
    else:
        outfile1 = open('%s/ContourData/%sOneSigmaContour.dat'%(outdir,fname),'w')
        outfile2 = open('%s/ContourData/%sTwoSigmaContour.dat'%(outdir,fname),'w')

    outfile1.write('# Energy PercetageUncertainty\n')
    outfile2.write('# Energy PercetageUncertainty\n')

    realup = ROOT.TGraph()
    realdown = ROOT.TGraph()
    realup2 = ROOT.TGraph()
    realdown2 = ROOT.TGraph()

    realup.SetPoint(0,emin,0)
    realdown.SetPoint(0,emin,0)
    realup2.SetPoint(0,emin,0)
    realdown2.SetPoint(0,emin,0)

    nPoint = 200
    step = (np.log10(emax)-np.log10(emin))/nPoint

    for i in range(0,nPoint+1):
        x = np.array([np.power(10,np.log10(emin)+i*step)])
        grad = np.empty(fit_func.GetNpar())
        fit_func.GradientPar(x,grad,0.001)
        sigmaF = CalculateSigmaF(fname,grad,covMat,x[0],primary,Barr)
        sigmaP = sigmaF/fit_func.Eval(x[0]) * 100
        realup.SetPoint(i+1,x[0],sigmaP)
        realdown.SetPoint(i+1,x[0],-sigmaP)
        realup2.SetPoint(i+1,x[0],2*sigmaP)
        realdown2.SetPoint(i+1,x[0],-2*sigmaP)
        outfile1.write('%.4f %.4f\n'%(x[0],sigmaP))
        outfile2.write('%.4f %.4f\n'%(x[0],2*sigmaP))

    realup.SetPoint(realup.GetN(),emax,0);
    realdown.SetPoint(realdown.GetN(),emax,0);
    realup.SetFillColor(390);
    realup.SetLineColor(390);
    realdown.SetFillColor(390);
    realdown.SetLineColor(390);

    realup2.SetPoint(realup2.GetN(),np.power(10,np.log10(emax)+step),0);
    realdown2.SetPoint(realdown2.GetN(),np.power(10,np.log10(emax)+step),0);
    realup2.SetFillColor(406);
    realup2.SetLineColor(406);
    realdown2.SetFillColor(406);
    realdown2.SetLineColor(406);

    outfile1.close()
    outfile2.close()

    return realup, realdown, realup2, realdown2

def DrawFinalDeviationPlot(dgs, dleg, U1, U2, D1, D2, emin, emax, primary, outdir, fname, extension='pdf'):

    # Set up style for canvas
    mystyle = ROOT.TStyle("Plain","Mystyle")
    mystyle.cd()
    mystyle.SetLabelSize(0.04,"XYZ")
    mystyle.SetLabelOffset(0.01,"XYZ")
    mystyle.SetTitleXSize(0.05)
    mystyle.SetTitleYSize(0.05)
    mystyle.SetTitleXOffset(1.1)
    mystyle.SetTitleYOffset(1.2)
    c1 = ROOT.TCanvas("c","c",800,800)
    c1.SetLogx()
    c1.SetTicks(0,1)
    c1.SetBottomMargin(0.12)
    c1.SetLeftMargin(0.13)
    c1.SetRightMargin(0.03)
    c1.SetTopMargin(0.03)
    c1.cd()

    # Draw the uncertainties first so they end up behind the data
    U2.Draw("AF")
    # U2.SetTitle("Experimental Cosmic Ray %s Flux"%primary)
    U2.GetYaxis().SetRangeUser(-100.0,100.0)
    U2.GetYaxis().SetTitle("Deviation from central fit value (%)")
    U2.GetXaxis().SetLimits(emin,emax)
    # H3a is parameterised as total cosmic ray energy
    if fname == 'H3a':
        U2.GetXaxis().SetTitle('Primary Energy [GeV]')
    # GSHL and power law are in terms of energy per nucleon
    else:
        U2.GetXaxis().SetTitle('Primary Energy [GeV/n]')
    U2.GetXaxis().CenterTitle()
    U2.GetYaxis().CenterTitle()
    D2.Draw("Fsame")
    U1.Draw("Fsame")
    D1.Draw("Fsame")

    # Draw all of the data
    # This isn't done as a multigraph because ROOT
    for dg in dgs:
        dg.Draw("PEsame")

    # Draw the legend
    dleg.Draw()

    # Draw a zero line to help guide the eye
    ZeroLine = ROOT.TLine(emin,0,emax,0)
    ZeroLine.Draw()

    # Redraw the axes so the ticks aren't covered by the contours
    ROOT.gPad.RedrawAxis()

    if outdir == os.getcwd():
        c1.SaveAs("%s/Global%sFitDeviation.%s"%(outdir,primary,extension))
    else:
        c1.SaveAs("%s/FinalPlots/Global%sFitDeviation.%s"%(outdir,primary,extension))

def PlotSaveCovCorMats(covMat,fit_func,outdir,primary,Plot=True,extension='pdf'):

    parNames = []
    for i in range(0,fit_func.GetNpar()):
        parNames.append(fit_func.GetParName(i))

    if outdir == os.getcwd():
        outfile = open('%s/CovCorMatrices.dat'%(outdir),'w')

    else:
        outfile = open('%s/FitReturnValues/CovCorMatrices.dat'%(outdir),'w')
        jsonOutfile1 = open('%s/JsonFormat/CovMatrix.json'%(outdir),'w')
        jsonOutfile2 = open('%s/JsonFormat/CorMatrix.json'%(outdir),'w')

    outfile.write('%s %s Fit\n'%(fit_func.GetName(),primary))
    outfile.write('\n')
    outfile.write('Covariance Matrix\n')
    outfile.write('\n')

    jsonMatrix1 = collections.OrderedDict()
    for i in range(0,fit_func.GetNpar()):
        jsonMatrix1[parNames[i]] = collections.OrderedDict()
        for j in range(0,fit_func.GetNpar()):
            outfile.write('%.4f '%(covMat[i][j]))
            jsonMatrix1[parNames[i]][parNames[j]] = covMat[i][j]
        outfile.write('\n')

    outfile.write('\n')
    outfile.write('Correlation Matrix\n')
    outfile.write('\n')
    json.dump(jsonMatrix1,jsonOutfile1)

    corMat = []
    jsonMatrix2 = collections.OrderedDict()
    for i in range(0,fit_func.GetNpar()):
        row = []
        jsonMatrix2[parNames[i]] = collections.OrderedDict()
        for j in range(0,fit_func.GetNpar()):
            corMatVal = covMat[i][j]/(np.sqrt(covMat[i][i]*covMat[j][j]))
            row.append(corMatVal)
            outfile.write('%.4f '%(corMatVal))
            jsonMatrix2[parNames[i]][parNames[j]] = corMatVal
        outfile.write('\n')
        corMat.append(row)
    json.dump(jsonMatrix2,jsonOutfile2)


    outfile.close()
    jsonOutfile1.close()
    jsonOutfile2.close()

    ParNames = []

    for i in range(0,fit_func.GetNpar()):
        ParNames.append(fit_func.GetParName(i))

    if Plot:

        plt.imshow(corMat,interpolation='none',cmap=cm.coolwarm,alpha=0.3,vmin=-1.0,vmax=1.0)
        plt.title('%s %s Correlation Matrix'%(fit_func.GetName(),primary))
        plt.colorbar()
        plt.xticks(np.arange(len(ParNames)),ParNames)
        plt.yticks(np.arange(len(ParNames)),ParNames)
        for i in range(0,len(corMat)):
            for j in range(0,len(corMat[0])):
                plt.text(i, j, '%.4f'%corMat[i][j],
                         verticalalignment='center',
                         horizontalalignment='center',
                         color='w',
                         path_effects=[PathEffects.withStroke(linewidth=3,
                                                              foreground='k')])

        if outdir == os.getcwd():
            plt.savefig('%s/%s%sCorrelationMatrix.%s'%(outdir,fit_func.GetName(),primary,extension))
        else:
            plt.savefig('%s/FitReturnValues/%s%sCorrelationMatrix.%s'%(outdir,fit_func.GetName(),primary,extension))

    plt.close()

def SaveGlobalFit(fit_func,outdir,primary):

    if outdir == os.getcwd():
        outfile = open('%s/FitParameters.dat'%(outdir),'w')

    else:
        outfile = open('%s/FitReturnValues/FitParameters.dat'%(outdir),'w')
        jsonOutfile1 = open('%s/JsonFormat/FitParameters.json'%(outdir),'w')
        jsonOutfile2 = open('%s/JsonFormat/FitResults.json'%(outdir),'w')


    outfile.write('%s %s Fit Parameters\n'%(fit_func.GetName(),primary))
    outfile.write('(Errors quoted as directly returned from ROOT)\n')
    outfile.write('\n')

    jsonData1 = collections.OrderedDict()
    jsonData2 = collections.OrderedDict()
    jsonData1['value'] = collections.OrderedDict()
    jsonData1['error'] = collections.OrderedDict()
    for i in range(0,fit_func.GetNpar()):
        jsonData1['value'][fit_func.GetParName(i)] = fit_func.GetParameter(i)
        jsonData1['error'][fit_func.GetParName(i)] = fit_func.GetParError(i)
        outfile.write('%s = %.4f +/- %.4f\n'%(fit_func.GetParName(i),
                                              fit_func.GetParameter(i),
                                              fit_func.GetParError(i)))

    outfile.write('\n')
    outfile.write('Chi2 = %.1f\n'%fit_func.GetChisquare())
    outfile.write('Number of Datapoints = %i\n'%fit_func.GetNumberFitPoints())
    outfile.write('Degrees of Freedom = %i\n'%fit_func.GetNDF())
    outfile.write('Reduced Chi2 = %.1f\n'%(fit_func.GetChisquare()/fit_func.GetNDF()))
    jsonData2['chi2'] = fit_func.GetChisquare()
    jsonData2['nDatapoints'] = fit_func.GetNumberFitPoints()
    jsonData2['dof'] = fit_func.GetNDF()
    jsonData2['redchi2'] = fit_func.GetChisquare()/float(fit_func.GetNDF())

    json.dump(jsonData1,jsonOutfile1)
    json.dump(jsonData2,jsonOutfile2)

    outfile.close()
    jsonOutfile1.close()
    jsonOutfile2.close()

if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tolerance',type=str,required=True,
                        help="""Either provide path to a DAT file containing
                        delta chi values from which to calculate the tolerance
                        or a number to use.""")
    parser.add_argument('--data',type=str,required=True,
                        help="""JSON file containing locations of all of the
                        data you want to include in the fit.""")
    parser.add_argument('--primary',type=str,default='Proton',
                        help="""Name of primary you are performing the global
                        fit on.""")
    parser.add_argument('--fname',type=str,default='GSHL',
                        help="""
                        Name of fitting function you want to use for the global
                        fits. The choices are:
                            * GSHL
                            * AMSGSHL (GSHL with AMS modification)
                            * Simple (Power Law)
                            * H3a
                        """)
    parser.add_argument('--Barr',action='store_true',default=False,
                        help="""Flag if performing paper a la Barr et al.""")
    parser.add_argument('--emin',type=float,default='0.5',
                        help="""Minimum energy value for fitter.""")
    parser.add_argument('--emax',type=float,default='250000',
                        help="""Maximum energy value for fitter.""")
    parser.add_argument('--smin',type=float,default='416',
                        help="""Minimum solar modulation value for which to
                        accept all data.""")
    parser.add_argument('--smax',type=float,default='555',
                        help="""Maximum solar modulation value for which to
                        accept all data.""")
    parser.add_argument('--outdir',type=str,default=None,
                        help="""Directory to save output. Script is designed to
                        create subdirectories and delete older plots if it needs
                        to, so set this carefully! If none provided, all plots
                        will be saved in CWD.""")
    parser.add_argument('--PNG',action='store_true',default=False,
                        help="""Flag if wanting to save plots as PNG rather
                        than the default of PDF.""")

    args = parser.parse_args()

    # Check requested function is valid. If not, exit.
    allowed_functions = ['GSHL','AMSGSHL','Simple','H3a']
    if args.fname not in allowed_functions:
        print "Requested fitting function is invalid. Please select from:"
        for allowed_function in allowed_functions:
            print "    %s"%allowed_function
        sys.exit()

    # Check file extension
    if args.PNG:
        extension = 'png'
    else:
        extension = 'pdf'

    #########################################################
    ###                                                   ###
    ### Make output directory structure if required       ###
    ###                                                   ###
    #########################################################

    if args.outdir is None:
        # Save everything in current working directory
        outdir = os.getcwd()
    else:
        # Create it if needs be
        if not os.path.isdir(args.outdir):
            os.makedirs(args.outdir)
        # Go on to create the subdirectories
        # Start with the parameterisation
        if not os.path.isdir(args.outdir+'/'+args.fname):
            os.makedirs(args.outdir+'/'+args.fname)
        # Then add the primary
        if not os.path.isdir(args.outdir+'/'+args.fname+'/'+args.primary):
            os.makedirs(args.outdir+'/'+args.fname+'/'+args.primary)
        # Then a subdirectory to store the resultant plots/data
        datadirs = ['FinalPlots','ContourData','FitReturnValues','JsonFormat']
        for datadir in datadirs:
            if not os.path.isdir(args.outdir+'/'+args.fname+'/'+args.primary+'/'+datadir):
                os.makedirs(args.outdir+'/'+args.fname+'/'+args.primary+'/'+datadir)
            # Delete the files which exist there, if any
            for outimg in os.listdir(args.outdir+'/'+args.fname+'/'+args.primary+'/'+datadir):
                os.unlink(os.path.join(args.outdir+'/'+args.fname+'/'+args.primary+'/'+datadir,outimg))
        # Set outdir to user's choice
        outdir = args.outdir+'/'+args.fname+'/'+args.primary

    #########################################################
    ###                                                   ###
    ### Calculate Tolerance                               ###
    ###                                                   ###
    #########################################################

    try:
        float(args.tolerance)
        tolerance = float(args.tolerance)

    except:
        chidata = np.genfromtxt(args.tolerance,
                                delimiter = " ",
                                names=True,
                                dtype=None)

        # Calculate delta chi sum
        totaldeltachi = np.sum(chidata['DeltaChi'])

        # Tolerance is then square root of this
        tolerance = np.sqrt(totaldeltachi)

    ROOT.TVirtualFitter.SetErrorDef(tolerance)

    #########################################################
    ###                                                   ###
    ### Load all Data in to a Dictionary                  ###
    ###                                                   ###
    #########################################################

    total_data_dict = {}
    fh = json.load(open(args.data))

    for data_key in fh.keys():
        total_data_dict[data_key], DataGood, DataExtracted = LoadData(fh[data_key],args.emin,args.emax,args.smin,args.smax,10.0,Barr=args.Barr)
        if not DataGood:
            sys.exit()
        if not DataExtracted:
            del total_data_dict[data_key]

    #########################################################
    ###                                                   ###
    ### Perform Regular Fit                               ###
    ###                                                   ###
    #########################################################

    # First initialise the desired fitting function
    fit_func, FuncGood = InitialiseFittingFunction(args.fname,args.emin,args.emax,args.primary)
    if not FuncGood:
        sys.exit()

    # Load all data from dictionary in to a multigraph object
    mg,leg = CreateFluxMultiGraph(total_data_dict)

    # Fit with user specified function
    # Save fitresult so we can access covariance matrix
    fitresult = mg.Fit(args.fname,"SR")
    covMat = fitresult.GetCovarianceMatrix()

    # Save covariance and correlation matrix to file
    # Correlation matrix will also be plotted if PathEffects could be imported
    PlotSaveCovCorMats(covMat,fit_func,outdir,args.primary,Plot=PlotMatrices,extension=extension)

    # Save these best fit parameters in a list
    globalFitPar = []
    for i in range(0,fit_func.GetNpar()):
        globalFitPar.append(fit_func.GetParameter(i))

    # Also save them to file
    SaveGlobalFit(fit_func,outdir,args.primary)

    # Calculate uncertainty contours for fit plot
    FitU1, FitD1, FitU2, FitD2 = GenerateFitUncertaintyLines(args.emin, args.emax,fit_func,args.fname,args.primary,covMat,args.Barr)
    # Create flux grah array with appropriate legend entries
    fgs, fleg = CreateFluxGraphArray(total_data_dict,U2=FitU2,U1=FitU1)

    # Plot and save this global fit
    PlotSaveFluxMultiGraph(fgs,fleg,FitU1,FitU2,FitD1,FitD2,args.primary,outdir,args.emin,args.emax,fit_func,args.fname,extension=extension)

    #########################################################
    ###                                                   ###
    ### Make Deviation Plot inc. 1 and 2 sigma bands      ###
    ###                                                   ###
    #########################################################

    # The 1 and 2 sigma bands are 4 different TGraphs
    DevU1, DevD1, DevU2, DevD2 = GenerateDeviationUncertaintyLines(args.emin,args.emax,fit_func,args.fname,outdir,args.primary,covMat,args.Barr)

    # Make all of the deviation graphs
    dgs, dleg = CreateDeviationGraphs(total_data_dict,fit_func,U2=DevU2,U1=DevU1)

    # Now draw it all
    DrawFinalDeviationPlot(dgs, dleg, DevU1, DevU2, DevD1, DevD2, args.emin, args.emax, args.primary, outdir, args.fname,extension=extension)
