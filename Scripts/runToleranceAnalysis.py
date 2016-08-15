
import os, sys, json
import numpy as np
import ROOT

from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter

def LoadData(filedict,emin,emax,smin,smax,lowcut,Barr=False):

    energycentres = []
    energylowedges = []
    energyhighedges = []
    fluxvalues = []
    fluxlowerrors = []
    fluxhigherrors = []

    DataGood = True

    try:
        filename = filedict['loc']
        sol_mod = filedict['sol_mod']
        try:
            open(filename)
        except:
            print "Something went wrong when opening %s"%filename
            DataGood = False
        for line in open(filename):
            if len(line.rstrip().split(' ')) == 10:
                Energy = float(line.rstrip().split(' ')[0])
                EnergyHighEnough = Energy >= emin
                EnergyLowEnough = Energy <= emax
                if sol_mod < smin or sol_mod > smax:
                    EnergyHighEnough = EnergyHighEnough and Energy >= lowcut
                if EnergyHighEnough and EnergyLowEnough:
                    energycentres.append(float(line.rstrip().split(' ')[0]))
                    energylowedges.append(float(line.rstrip().split(' ')[1]))
                    energyhighedges.append(float(line.rstrip().split(' ')[2]))
                    fluxvalues.append(float(line.rstrip().split(' ')[3]))
                    fluxlowerrors.append(float(line.rstrip().split(' ')[-2]))
                    fluxhigherrors.append(float(line.rstrip().split(' ')[-1]))

    except:
        for filenum in filedict.keys():
            try:
                int(filenum)
                filename = filedict[filenum]['loc']
                sol_mod = filedict[filenum]['sol_mod']
                if Barr:
                    if 'BESS-TeV' in filename:
                        lowcut = 120.0
                    else:
                        lowcut = 10.0
                try:
                    open(filename)
                except:
                    print "Something went wrong when opening %s"%filename
                    DataGood = False
                for line in open(filename):
                    if len(line.rstrip().split(' ')) == 10:
                        Energy = float(line.rstrip().split(' ')[0])
                        EnergyHighEnough = Energy >= emin
                        EnergyLowEnough = Energy <= emax
                        if sol_mod < smin or sol_mod > smax:
                            EnergyHighEnough = EnergyHighEnough and Energy > lowcut
                        if EnergyHighEnough and EnergyLowEnough:
                            energycentres.append(float(line.rstrip().split(' ')[0]))
                            energylowedges.append(float(line.rstrip().split(' ')[1]))
                            energyhighedges.append(float(line.rstrip().split(' ')[2]))
                            fluxvalues.append(float(line.rstrip().split(' ')[3]))
                            fluxlowerrors.append(float(line.rstrip().split(' ')[-2]))
                            fluxhigherrors.append(float(line.rstrip().split(' ')[-1]))
            except:
                continue

    energylowerrors = np.array(energycentres)-np.array(energylowedges)
    energyhigherrors = np.array(energyhighedges)-np.array(energycentres)

    data_dict = {}
    data_dict['en'] = energycentres
    data_dict['en_low'] = energylowerrors
    data_dict['en_high'] = energyhigherrors
    data_dict['flux'] = fluxvalues
    data_dict['flux_low'] = fluxlowerrors
    data_dict['flux_high'] = fluxhigherrors
    data_dict['expname'] = filedict['expname']
    data_dict['years'] = filedict['years']
    data_dict['line_colour'] = filedict['line_colour']
    data_dict['marker_colour'] = filedict['marker_colour']
    data_dict['marker_style'] = filedict['marker_style']

    if len(energycentres) == 0:
        DataExtracted = False
    else:
        DataExtracted = True

    return data_dict, DataGood, DataExtracted

def InitialiseFittingFunction(fname,emin,emax,primary):

    if fname == 'GSHL':
        func = ROOT.TF1("GSHL","[0]*((x)+[1]*exp([2]*sqrt(x)))^(-1*[3])",emin,emax)
        func.SetParName(0,"a")
        func.SetParName(1,"b")
        func.SetParName(2,"c")
        func.SetParName(3,"d")
        
        # Set initial values from Barr et al paper
        # [arXiv:astro-ph/0611266]
        if primary == 'Proton':
            func.SetParameter(0,14900)
            func.SetParameter(1,2.15)
            func.SetParameter(2,-0.21)
            func.SetParameter(3,2.74)
            FuncGood = True
        elif primary == 'Helium':
            func.SetParameter(0,460)
            func.SetParameter(1,1.16)
            func.SetParameter(2,-0.33)
            func.SetParameter(3,2.60)
            FuncGood = True
        # Values for higher energy primaries done on intuition
        elif primary == 'Nitrogen':
            func.SetParameter(0,10)
            func.SetParameter(1,1.16)
            func.SetParameter(2,-0.33)
            func.SetParameter(3,2.60)
            FuncGood = True
        else:
            print "%s currently not supported with this parameterisation. Please choose another"%primary
            FuncGood = False

    elif fname == 'Simple':
        func = ROOT.TF1("Simple","[0]*x^(-1*[1])",emin,emax)
        func.SetParName(0,"a")
        func.SetParName(1,"d")
        
        # Set initial values from Barr et al paper
        # [arXiv:astro-ph/0611266]
        # Ignore b and c values since we neglect them here
        if primary == 'Proton':
            func.SetParameter(0,14900)
            func.SetParameter(1,2.74)
            FuncGood = True
        elif primary == 'Helium':
            func.SetParameter(0,460)
            func.SetParameter(1,2.60)
            FuncGood = True
        # Values for higher energy primaries done on intuition
        elif primary == 'Nitrogen':
            func.SetParameter(0,10)
            func.SetParameter(3,2.60)
            FuncGood = True
        else:
            print "%s currently not supported with this parameterisation. Please choose another"%primary
            FuncGood = False

    elif fname == 'AMSGSHL':
        # Inspired by AMS-02 paper. Joins GSHL and their suggested
        # alternate parameterisation. Initial values from Davide's work.
        func = ROOT.TF1("AMSGSHL", "[0] * ((x)+[1]*exp([2]*sqrt(x)))^(-1*[3]) * (1+(x/[4])^([5]))^(([3]-[6])/[5])",emin,emax);
        func.SetParName(0,"a")
        func.SetParName(1,"b")
        func.SetParName(2,"c")
        func.SetParName(3,"d")
        func.SetParName(4,"k")
        func.SetParName(5,"s")
        func.SetParName(6,"e")
        if primary == 'Proton':
            func.SetParameter(0,14900)
            func.SetParameter(1,2.15)
            func.SetParameter(2,-0.21)
            func.SetParameter(3,2.74)
            func.SetParameter(4,500)
            func.SetParameter(5,1)
            func.SetParameter(6,2.68)
            FuncGood = True
        elif primary == 'Helium':
            func.SetParameter(0,460)
            func.SetParameter(1,1.16)
            func.SetParameter(2,-0.33)
            func.SetParameter(3,2.60)
            func.SetParameter(4,300)
            func.SetParameter(5,1)
            func.SetParameter(6,2.58)
            FuncGood = True
        else:
            print "%s currently not supported with this parameterisation. Please choose another"%primary
            FuncGood = False

    elif fname == 'H3a':
        # Inspired by Gaisser paper
        # Gaisser, T.K., Astroparticle Physics 35, 801 (2012)
        # Third component ignored, since it only comes in at very high energies

        if primary == 'Proton':
            func = ROOT.TF1("H3a","[0]*x^(-1*[1])*exp(-x/(4.0*1000000))+[3]*x^(-1*[4])*exp(-x/(30.0*1000000))",emin,emax)
            func.SetParameter(0,7860)
            func.SetParameter(1,1.66)
            func.SetParameter(2,20.0)
            func.SetParameter(3,1.4)
            FuncGood = True
        elif primary == 'Helium':
            func = ROOT.TF1("H3a","[0]*x^(-1*[1])*exp(-x/(8.0*1000000))+[3]*x^(-1*[4])*exp(-x/(60.0*1000000))",emin,emax)
            func.SetParameter(0,3550)
            func.SetParameter(1,1.58)
            func.SetParameter(2,20.0)
            func.SetParameter(3,1.4)
            FuncGood = True
        elif primary == 'Carbon' or primary == 'Nitrogen' or primary == 'Oxygen':
            func = ROOT.TF1("H3a","[0]*x^(-1*[1])*exp(-x/(28.0*1000000))+[3]*x^(-1*[4])*exp(-x/(210.0*1000000))",emin,emax)
            func.SetParameter(0,2200)
            func.SetParameter(1,1.63)
            func.SetParameter(2,13.4)
            func.SetParameter(3,1.4)
            FuncGood = True
        elif primary == 'Magnesium' or primary == 'Silicon':
            func = ROOT.TF1("H3a","[0]*x^(-1*[1])*exp(-x/(52.0*1000000))+[3]*x^(-1*[4])*exp(-x/(390.0*1000000))",emin,emax)
            func.SetParameter(0,1430)
            func.SetParameter(1,1.67)
            func.SetParameter(2,13.4)
            func.SetParameter(3,1.4)
            FuncGood = True
        elif primary == 'Iron':
            func = ROOT.TF1("H3a","[0]*x^(-1*[1])*exp(-x/(104.0*1000000))+[3]*x^(-1*[4])*exp(-x/(780.0*1000000))",emin,emax)
            func.SetParameter(0,2120)
            func.SetParameter(1,1.63)
            func.SetParameter(2,13.4)
            func.SetParameter(3,1.4)
            FuncGood = True
        else:
            print "%s currently not supported with this parameterisation. Please choose another"%primary
            FuncGood = False

        func.SetParName(0,"a1")
        func.SetParName(1,"y1")
        func.SetParName(2,"a2")
        func.SetParName(3,"y2")
        
    return func, FuncGood

def CreateFluxMultiGraph(data,omit_key=None):

    mg = ROOT.TMultiGraph()
    leg = ROOT.TLegend(0.5, 0.6, .89, .89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    for data_key in data.keys():
        if omit_key is None or omit_key != data_key:
            graph = ROOT.TGraphAsymmErrors()
            for j in range(0,len(data[data_key]['en'])):
                graph.SetPoint(j,
                               data[data_key]['en'][j],
                               data[data_key]['flux'][j])
                graph.SetPointError(j,0,0,
                                    data[data_key]['flux_low'][j],
                                    data[data_key]['flux_high'][j])
            graph.SetLineColor(data[data_key]['line_colour'])
            graph.SetMarkerColor(data[data_key]['marker_colour'])
            graph.SetLineWidth(1)
            graph.SetMarkerSize(1)
            graph.SetMarkerStyle(data[data_key]['marker_style'])
            graph.SetFillColor(0)
            mg.Add(graph)
            leg.AddEntry(graph,"%s %s"%(data[data_key]['expname'],
                                        data[data_key]['years']))

    return mg, leg

def CreateFluxSingleGraph(data,include_key=None):

    graph = ROOT.TGraphAsymmErrors()
    leg = ROOT.TLegend(0.5, 0.6, .89, .89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    print include_key
    
    for data_key in data.keys():
        if data_key == include_key:
            for j in range(0,len(data[data_key]['en'])):
                graph.SetPoint(j,
                               data[data_key]['en'][j],
                               data[data_key]['flux'][j])
                graph.SetPointError(j,0,0,
                                    data[data_key]['flux_low'][j],
                                    data[data_key]['flux_high'][j])
            graph.SetLineColor(data[data_key]['line_colour'])
            graph.SetMarkerColor(data[data_key]['marker_colour'])
            graph.SetLineWidth(1)
            graph.SetMarkerSize(1)
            graph.SetMarkerStyle(data[data_key]['marker_style'])
            graph.SetFillColor(0)
            leg.AddEntry(graph,"%s %s"%(data[data_key]['expname'],
                                        data[data_key]['years']))

    return graph, leg

def CreateDeviationGraphs(data,fit_func,omit_key=None,U2=None,U1=None):

    graphs = []
    leg = ROOT.TLegend(0.16, 0.15, 0.6, 0.3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(2)

    if U2 is not None:
        leg.AddEntry(U1,"1 sigma contour")
        leg.AddEntry(U2,"2 sigma contour")

    for data_key in data.keys():
        if omit_key is None or omit_key != data_key:
            graph = ROOT.TGraphAsymmErrors()
            for j in range(0,len(data[data_key]['en'])):
                fvalue = fit_func.Eval(data[data_key]['en'][j])
                yvalue = (data[data_key]['flux'][j] - fvalue)/fvalue*100.0
                yerrorup = np.abs((data[data_key]['flux'][j]+data[data_key]['flux_high'][j]-fvalue)/fvalue*100.0 - yvalue)
                yerrorlow = np.abs((data[data_key]['flux'][j]-data[data_key]['flux_low'][j]-fvalue)/fvalue*100.0 - yvalue)
                graph.SetPoint(j,data[data_key]['en'][j],yvalue)
                graph.SetPointError(j,0,0,yerrorlow,yerrorup)
            graph.SetLineColor(data[data_key]['line_colour'])
            graph.SetMarkerColor(data[data_key]['marker_colour'])
            graph.SetLineWidth(1)
            graph.SetMarkerSize(1)
            graph.SetMarkerStyle(data[data_key]['marker_style'])
            graph.SetFillStyle(0)
            graph.SetFillColor(0)
            graphs.append(graph)
            leg.AddEntry(graph,"%s %s"%(data[data_key]['expname'],
                                        data[data_key]['years']))

    return graphs, leg

def CreateDeviationGraph(data,fit_func,plot_key):

    leg = ROOT.TLegend(0.16, 0.15, 0.6, 0.3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    for data_key in data.keys():
        if plot_key == data_key:
            graph = ROOT.TGraphAsymmErrors()
            for j in range(0,len(data[data_key]['en'])):
                fvalue = fit_func.Eval(data[data_key]['en'][j])
                yvalue = (data[data_key]['flux'][j] - fvalue)/fvalue*100.0
                yerrorup = np.abs((data[data_key]['flux'][j]+data[data_key]['flux_high'][j]-fvalue)/fvalue*100.0 - yvalue)
                yerrorlow = np.abs((data[data_key]['flux'][j]-data[data_key]['flux_low'][j]-fvalue)/fvalue*100.0 - yvalue)
                graph.SetPoint(j,data[data_key]['en'][j],yvalue)
                graph.SetPointError(j,0,0,yerrorlow,yerrorup)
            graph.SetLineColor(data[data_key]['line_colour'])
            graph.SetMarkerColor(data[data_key]['marker_colour'])
            graph.SetLineWidth(1)
            graph.SetMarkerSize(1)
            graph.SetMarkerStyle(data[data_key]['marker_style'])
            graph.SetFillStyle(0)
            graph.SetFillColor(0)
            leg.AddEntry(graph,"%s %s"%(data[data_key]['expname'],
                                        data[data_key]['years']))

    return graph, leg

def Method1FitInspection(graph,leg,primary,data_to_plot,outdir,fname,emin,emax,fit_func):

    canvas = ROOT.TCanvas("c2")
    graph.SetTitle('Experimental Cosmic Ray %s Flux'%(primary))
    graph.Draw("AP")
    xmin = '1E%i'%(int(np.log10(emin)))
    xmax = '1E%i'%(int(np.log10(emax))+1)
    graph.GetXaxis().SetLimits(float(xmin),float(xmax))
    ymin = '1E%i'%(int(np.log10(fit_func.Eval(emax)))-1)
    ymax = '1E%i'%(int(np.log10(fit_func.Eval(emin)))+1)
    graph.GetYaxis().SetRangeUser(float(ymin),float(ymax))

    # H3a is parameterised as total cosmic ray energy
    if fname == 'H3a':
        graph.GetXaxis().SetTitle('Primary Energy [GeV]')
    # GSHL and power law are in terms of energy per nucleon
    else:
        graph.GetXaxis().SetTitle('Primary Energy [GeV/n]')
    graph.GetYaxis().SetTitle('Flux [ (GeV/n m^{2} s sr)^{-1} ]')
    canvas.SetLogy()
    canvas.SetLogx()
    leg.Draw()
    pt = ROOT.TPaveText(0.15,0.15,0.50,0.30,"NBNDC")
    pt.AddText("Parameters fixed to Global")
    SaveName = 'Global'
    pt.SetFillStyle(0)
    pt.SetLineColor(0)
    pt.Draw()
    if outdir == os.getcwd():
        canvas.SaveAs("%s/%sGlobal%sFitParameters.pdf"%(outdir,data_to_plot,primary))
    else:
        canvas.SaveAs("%s/%s/%s/ToleranceMethod1/TolerancePlots/%sGlobal%sFitParameters.pdf"%(outdir,fname,primary,data_to_plot,primary))
    canvas.Close()

def Method1DeviationInspection(total_data_dict,fit_func,data_to_plot,primary,emin,emax,outdir,fname,Global=False):
    
    dg, dleg = CreateDeviationGraph(total_data_dict,fit_func,data_to_plot)

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
    c1.SetTopMargin(0.10)
    c1.cd()

    # Draw the data
    dg.Draw("APE")
    dg.SetTitle("Experimental Cosmic Ray %s Flux"%primary)
    dg.GetYaxis().SetRangeUser(-100.0,100.0)
    dg.GetYaxis().SetTitle("Deviation from central fit value (%)")
    dg.GetXaxis().SetLimits(emin,emax)
    # H3a is parameterised as total cosmic ray energy
    if fname == 'H3a':
        dg.GetXaxis().SetTitle('Primary Energy [GeV]')
    # GSHL and power law are in terms of energy per nucleon
    else:
        dg.GetXaxis().SetTitle('Primary Energy [GeV/n]')
    dg.GetXaxis().CenterTitle()
    dg.GetYaxis().CenterTitle()

    # Draw the legend
    dleg.Draw()

    # Draw a zero line to help guide the eye
    ZeroLine = ROOT.TLine(emin,0,emax,0)
    ZeroLine.Draw()

    pt = ROOT.TPaveText(0.15,0.73,0.50,0.88,"NBNDC")
    pt.AddText("Parameters fixed to Global")
    pt.SetFillStyle(0)
    pt.SetLineColor(0)
    pt.Draw()

    # Redraw the axes so the ticks aren't covered
    ROOT.gPad.RedrawAxis()

    if outdir == os.getcwd():
        c1.SaveAs("%s/%sGlobal%sFitParametersDeviation.pdf"%(outdir,data_to_plot,SaveName,primary))
    else:
        c1.SaveAs("%s/%s/%s/ToleranceMethod1/TolerancePlots/%sGlobal%sFitParametersDeviation.pdf"%(outdir,fname,primary,data_to_plot,primary))
    c1.Close()

def Method2FitInspection(mg,leg,primary,data_to_omit,outdir,fname,emin,emax,fit_func,Global=False):

    canvas = ROOT.TCanvas("c2")
    mg.SetTitle('Experimental Cosmic Ray %s Flux'%(primary))
    mg.Draw("AP")
    xmin = '1E%i'%(int(np.log10(emin)))
    xmax = '1E%i'%(int(np.log10(emax))+1)
    mg.GetXaxis().SetLimits(float(xmin),float(xmax))
    ymin = '1E%i'%(int(np.log10(fit_func.Eval(emax)))-1)
    ymax = '1E%i'%(int(np.log10(fit_func.Eval(emin)))+1)
    mg.GetYaxis().SetRangeUser(float(ymin),float(ymax))

    # H3a is parameterised as total cosmic ray energy
    if fname == 'H3a':
        mg.GetXaxis().SetTitle('Primary Energy [GeV]')
    # GSHL and power law are in terms of energy per nucleon
    else:
        mg.GetXaxis().SetTitle('Primary Energy [GeV/n]')
    mg.GetYaxis().SetTitle('Flux [ (GeV/n m^{2} s sr)^{-1} ]')
    canvas.SetLogy()
    canvas.SetLogx()
    leg.Draw()
    pt = ROOT.TPaveText(0.15,0.15,0.50,0.30,"NBNDC")
    pt.AddText("All data minus %s"%data_to_omit)
    if Global:
        pt.AddText("Parameters fixed to Global")
        SaveName = 'Global'
    else:
        pt.AddText("Parameters freely fit")
        SaveName = 'Local'
    pt.SetFillStyle(0)
    pt.SetLineColor(0)
    pt.Draw()
    if outdir == os.getcwd():
        canvas.SaveAs("%s/AllDataMinus%s%s%sFitParameters.pdf"%(outdir,data_to_omit,SaveName,primary))
    else:
        canvas.SaveAs("%s/%s/%s/ToleranceMethod2/TolerancePlots/AllDataMinus%s%s%sFitParameters.pdf"%(outdir,fname,primary,data_to_omit,SaveName,primary))
    canvas.Close()

def Method2DeviationInspection(total_data_dict,fit_func,data_to_omit,primary,emin,emax,outdir,fname,Global=False):
    
    dgs, dleg = CreateDeviationGraphs(total_data_dict,fit_func,omit_key=data_to_omit)

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
    c1.SetTopMargin(0.10)
    c1.cd()

    # Draw all of the data
    for i,dg in enumerate(dgs):
        if i == 0:
            dg.Draw("APE")
            dg.SetTitle("Experimental Cosmic Ray %s Flux"%primary)
            dg.GetYaxis().SetRangeUser(-100.0,100.0)
            dg.GetYaxis().SetTitle("Deviation from central fit value (%)")
            dg.GetXaxis().SetLimits(emin,emax)
            # H3a is parameterised as total cosmic ray energy
            if fname == 'H3a':
                dg.GetXaxis().SetTitle('Primary Energy [GeV]')
            # GSHL and power law are in terms of energy per nucleon
            else:
                dg.GetXaxis().SetTitle('Primary Energy [GeV/n]')
            dg.GetXaxis().CenterTitle()
            dg.GetYaxis().CenterTitle()
        else:
            dg.Draw("PEsame")

    # Draw the legend
    dleg.Draw()

    # Draw a zero line to help guide the eye
    ZeroLine = ROOT.TLine(emin,0,emax,0)
    ZeroLine.Draw()

    pt = ROOT.TPaveText(0.15,0.73,0.50,0.88,"NBNDC")
    pt.AddText("All data minus %s"%data_to_omit)
    if Global:
        pt.AddText("Parameters fixed to Global")
        SaveName = 'Global'
    else:
        pt.AddText("Parameters freely fit")
        SaveName = 'Local'
    pt.SetFillStyle(0)
    pt.SetLineColor(0)
    pt.Draw()

    # Redraw the axes so the ticks aren't covered
    ROOT.gPad.RedrawAxis()

    if outdir == os.getcwd():
        c1.SaveAs("%s/AllDataMinus%s%sFitParametersDeviation.pdf"%(outdir,data_to_omit,SaveName,primary))
    else:
        c1.SaveAs("%s/%s/%s/ToleranceMethod2/TolerancePlots/AllDataMinus%s%s%sFitParametersDeviation.pdf"%(outdir,fname,primary,data_to_omit,SaveName,primary))
    c1.Close()
        
if __name__ == '__main__':
    
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
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
                        help="Flag if performing paper a la Barr et al.")
    parser.add_argument('--emin',type=float,default='0.5',
                        help="""Minimum energy value for fitter""")
    parser.add_argument('--emax',type=float,default='250000',
                        help="""Maximum energy value for fitter""")
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

    args = parser.parse_args()

    # Check requested function is valid. If not, exit.
    allowed_functions = ['GSHL','AMSGSHL','Simple','H3a',]
    if args.fname not in allowed_functions:
        print "Requested fitting function is invalid. Please select from:"
        for allowed_function in allowed_functions:
            print "    %s"%allowed_function
        sys.exit()

    #########################################################
    ###                                                   ###
    ### Make output directory structure if required       ###
    ###                                                   ###
    #########################################################

    if args.outdir is None:
        # Save everything in current working directory
        outdir = os.getcwd()
    else:
        # Set outdir to user's choice
        outdir = args.outdir
        # Create it if needs be
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        # Go on to create the subdirectories
        # Start with the parameterisation
        if not os.path.isdir(outdir+'/'+args.fname):
            os.makedirs(outdir+'/'+args.fname)
        # Then add the primary
        if not os.path.isdir(outdir+'/'+args.fname+'/'+args.primary):
            os.makedirs(outdir+'/'+args.fname+'/'+args.primary)
        # Then we make subdirectories for each of the tolerance methods
        primdirs = ['ToleranceMethod1','ToleranceMethod2']
        # And for what is output from those
        subdirs = ['ToleranceData','TolerancePlots']
        for primdir in primdirs:
            # If it doesn't exist make it
            if not os.path.isdir(outdir+'/'+args.fname+'/'+args.primary+'/'+primdir):
                os.makedirs(outdir+'/'+args.fname+'/'+args.primary+'/'+primdir)
            for subdir in subdirs:
                # Again make it if it doesn't exist
                if not os.path.isdir(outdir+'/'+args.fname+'/'+args.primary+'/'+primdir+'/'+subdir):
                    os.makedirs(outdir+'/'+args.fname+'/'+args.primary+'/'+primdir+'/'+subdir)
                # Delete the files which exist there, if any
                for outimg in os.listdir(outdir+'/'+args.fname+'/'+args.primary+'/'+primdir+'/'+subdir):
                    os.unlink(os.path.join(outdir+'/'+args.fname+'/'+args.primary+'/'+primdir+'/'+subdir,outimg))

    #########################################################
    ###                                                   ###
    ### Load all Data in to a Dictionary                  ###
    ###                                                   ###
    #########################################################

    total_data_dict = {}
    fh = json.load(open(args.data))

    for data_key in fh.keys():
        # 'lowcut' refers to cut on solar modulation. If in correct range,
        # there is no low energy cut. Otherwise, it is 10 GeV.
        # If matching Barr paper (given by args.Barr) then there is a special
        # cut on the BESS-TeV data of 120 GeV.
        total_data_dict[data_key], DataGood, DataExtracted = LoadData(fh[data_key],args.emin,args.emax,args.smin,args.smax,10.0,Barr=args.Barr)
        if not DataGood:
            sys.exit()
        if not DataExtracted:
            del total_data_dict[data_key]

    #########################################################
    ###                                                   ###
    ### Perform Initial Global Fit                        ###
    ###                                                   ###
    #########################################################

    # First initialise the desired fitting function
    fit_func, FuncGood = InitialiseFittingFunction(args.fname,args.emin,args.emax,args.primary)
    if not FuncGood:
        sys.exit()

    # Load all data from dictionary in to a multigraph object
    mg,leg = CreateFluxMultiGraph(total_data_dict)

    # Fit with user specified function
    mg.Fit(args.fname,"R")

    # Save these best fit parameters in a list
    globalFitPar = []
    for i in range(0,fit_func.GetNpar()):
        globalFitPar.append(fit_func.GetParameter(i))

    #########################################################
    ###                                                   ###
    ### Perform Tolerance Calculation Method 1            ###
    ###                                                   ###
    #########################################################

    if args.outdir is None:
        outfile1 = open('%s/Method1ChiVsN.dat'%outdir,'w')
    else:
        outfile1 = open('%s/%s/%s/ToleranceMethod1/ToleranceData/Method1ChiVsN.dat'%(outdir,args.fname,args.primary),'w')
        
    outfile1.write('# ExperimentNumber ExperimentFile NumData Chi DeltaChiN Err\n')

    for v,data_key in enumerate(total_data_dict.keys()):
        graph, leg = CreateFluxSingleGraph(total_data_dict,data_key)
        # Fit this dataset with parameters fixed to global fit
        for i in range(0,fit_func.GetNpar()):
            fit_func.FixParameter(i,globalFitPar[i])
        FitResult = graph.Fit(args.fname,"SR")
        num_datapoints = FitResult.Ndf() - FitResult.NFreeParameters()
        chi = fit_func.GetChisquare()

        # Make a plot of this to inspect
        Method1FitInspection(graph,leg,args.primary,data_key,outdir,args.fname,args.emin,args.emax,fit_func)

        # Make a plot of the deviations to also inspect
        Method1DeviationInspection(total_data_dict,fit_func,data_key,args.primary,args.emin,args.emax,outdir,args.fname)
        
        outfile1.write('%i %s %i %.4f %.4f %.4f\n'%(v,
                                                    data_key,
                                                    num_datapoints,
                                                    chi,
                                                    np.abs(num_datapoints-chi),
                                                    np.sqrt(2*num_datapoints)))

    outfile1.close()

    #########################################################
    ###                                                   ###
    ### Perform Tolerance Calculation Method 2            ###
    ###                                                   ###
    #########################################################

    if len(total_data_dict.keys()) > 1:

        if args.outdir is None:
            outfile2 = open('%s/Method2DeltaChi.dat'%outdir,'w')
        else:
            outfile2 = open('%s/%s/%s/ToleranceMethod2/ToleranceData/Method2DeltaChi.dat'%(outdir,args.fname,args.primary),'w')
        outfile2.write('# ExperimentNumber ExperimentName Chi0 ChiMin DeltaChi\n')

        for v,data_to_omit in enumerate(total_data_dict.keys()):

            # Create MultiGraph with one set of data excluded
            mg, leg = CreateFluxMultiGraph(total_data_dict,data_to_omit)

            # Fit this dataset with parameters fixed to global fit
            for i in range(0,fit_func.GetNpar()):
                fit_func.FixParameter(i,globalFitPar[i])
            mg.Fit(args.fname,"R")
            chi0 = fit_func.GetChisquare()

            # Make a plot of this to inspect
            Method2FitInspection(mg,leg,args.primary,data_to_omit,outdir,args.fname,args.emin,args.emax,fit_func,Global=True)

            # Make a plot of the deviations to also inspect
            Method2DeviationInspection(total_data_dict,fit_func,data_to_omit,args.primary,args.emin,args.emax,outdir,args.fname,Global=True)

            # Fit this dataset freely
            for i in range(0,fit_func.GetNpar()):
                fit_func.ReleaseParameter(i)
            mg.Fit(args.fname,"R")
            chiMin = fit_func.GetChisquare()

            # Make a plot of this to inspect
            Method2FitInspection(mg,leg,args.primary,data_to_omit,outdir,args.fname,args.emin,args.emax,fit_func)

            # Make a plot of the deviations to also inspect
            Method2DeviationInspection(total_data_dict,fit_func,data_to_omit,args.primary,args.emin,args.emax,outdir,args.fname)

            # Calculate the chi2 difference
            deltachi = chi0-chiMin

            # Write all of this out to file
            outfile2.write('%i %s %.4f %.4f %.4f\n'%(v,
                                                     data_to_omit,
                                                     chi0,
                                                     chiMin,
                                                     deltachi))

        outfile2.close()

    else:

        print "Only one dataset has been passed to the global fitter, so the second tolerance calculation method cannot be performed."
