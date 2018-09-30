import ROOT as rt
from ROOT import *
import os
import sys
import re


fitMethod = ['1gauss','2gauss']
fitm = '1gauss'
if not os.path.exists('1gauss'): 
    os.makedirs('1gauss')
histos = {'M_JPsiPicut4.1-4.2_total':'4.1-4.2','M_JPsiPicut4.2-4.25_total':'4.2-4.25','M_JPsiPicut4.25-4.3_total':'4.25-4.3','M_JPsiPicut4.3-4.4_total':'4.3-4.4','M_JPsiPicut4.4-4.7_total':'4.4-4.7','M_JPsiPicut4.7-5.0_total':'4.7-5.0'}
hval= []
herrHi =[]
herrLo =[]
xval= [4.225,4.275,4.85,4.55,4.35,4.15]
xerrHi= [0,0,0,0,0,0]
xerrLo=[0,0,0,0,0,0]

w = RooWorkspace("w")
##This is the initial version
def main():
    if os.path.exists(fitm+"/signalval.txt"): 
        os.remove(fitm+"/signalval.txt")
    sigfile = open(fitm+"/signalval.txt",'a')
    f =  TFile(path)
    for h in histos:
        w.factory("Voigtian::bwgauss1(mass[3.9,3.6,4.04],mean3900[3.886,3.86,3.91],width3900[0.0282,0.0281,0.0283],sigma3900[0.011,0.009,0.013])")
        w.factory("Voigtian::bwgauss2(mass[3.9,3.6,4.05],mean3730[3.73,3.7,3.76],width3730[0.0282,0.026,0.03],sigma3730[0.011,0.001,0.013])")
        if '4.1-4.2' in h:
            bin =4.15
            w.var('mean3900').setRange(3.86,3.95)
            w.var('mass').setRange(3.6,4.04)
            w.var('sigma3900').setRange(0.005,0.05)
        if '4.2-4.25' in h:
            bin = 4.225
            w.var('mass').setRange(3.62,4.06)
        if '4.25-4.3' in h:
            bin = 4.275
            w.var('mass').setRange(3.63,4.1)
        if '4.3-4.4' in h:
            bin = 4.35
            w.var('mass').setRange(3.72,4.18)
            w.var('sigma3900').setRange(0.005,0.05)
        if '4.4-4.7' in h:
            bin = 4.55
            w.var('mass').setRange(3.74,4.4)
        if '4.7-5.0' in h:
            bin =4.85
            w.var('mass').setRange(3.8,4.4)
        c1 =  RooRealVar("a1","a1",-20000,20000);
        c2 =  RooRealVar("a2","a2",-20000,20000);
        c3 =  RooRealVar("a3","a3",-20000,20000);
        c4 =  RooRealVar("a4","a4",-20000,20000);
        #a5 =  RooRealVar("a5","a5",-20000,20000);
        #a6 =  RooRealVar("a6","a6",-20000,20000);
        #w.factory("Chebychev::ch(mass[3.9,3.5,4.3],RooArgList(a1,a2,a3,a4,a5,a6))")
        ch = RooChebychev("ch","ch",w.var('mass'),RooArgList(c1,c2,c3,c4))
        getattr(w,'import')(ch)
        if fitm is '2gauss' :
            w.factory("SUM::modelsum(nsig3900[20,0,2000]*bwgauss1,nsig3730[20,0,2000]*bwgauss2,nbkg[1000,0,50000]*ch)")
        if fitm is '1gauss' :
            w.factory("SUM::modelsum(nsig3900[20,0,5000]*bwgauss1,nbkg[1000,0,50000]*ch)")
        h1 =  TH1F()
        f.GetObject('histos/'+h,h1)
        print h1
        data =  RooDataHist("data","mydata", RooArgList(w.var('mass')),h1)
        print data
        result = w.pdf("modelsum").fitTo(data,rt.RooFit.Save())
        print result
        result.Print()
        paralist = result.floatParsFinal()
        signal3900 = paralist.at(paralist.index('nsig3900'))
        print signal3900
        hval.append(signal3900.getValV())
        herrHi.append(signal3900.getErrorHi())
        herrLo.append(signal3900.getErrorLo())
        print signal3900.getErrorLo()
        print fitm
        xframe = w.var('mass').frame()
        xframe2= w.var('mass').frame()
        xframe.SetMinimum(0)     # set minimum y-axis value
        xframe.SetMaximum(80)   # set maximum y-axis value
        xframe.SetTitle("M_{J/#psi #pi} with M_{J/#psi #pi^{+}#pi^{-}} cut"+histos[h])
        xframe2.SetTitle("Parameters Table")
        xframe2.SetAxisColor(0,"X")
        xframe2.SetAxisColor(0,"Y")
        xframe2.SetLabelColor(0,"X")
        xframe2.SetLabelColor(0,"Y")
        xframe2.SetXTitle("")
        data.plotOn(xframe)
        w.pdf("modelsum").plotOn(xframe)
        w.pdf("modelsum").paramOn(xframe2,rt.RooFit.Layout(0.1, 0.99,0.9))
        w.pdf("modelsum").plotOn(xframe,rt.RooFit.Components("ch"),rt.RooFit.LineStyle(kDashed),rt.RooFit.LineColor(kBlue),rt.RooFit.Name("bkg."))
        w.pdf("modelsum").plotOn(xframe,rt.RooFit.Components("bwgauss1"),rt.RooFit.LineStyle(kDashed),rt.RooFit.LineColor(kRed),rt.RooFit.Name("sig1."))
        if '2' in fitm :
            w.pdf("modelsum").plotOn(xframe,rt.RooFit.Components("bwgauss2"),rt.RooFit.LineStyle(kDashed),rt.RooFit.LineColor(kRed-1),rt.RooFit.Name("sig2."))
        c1=TCanvas("c1","c1",800,400)
        c1.Divide(2)
        c1.cd(1)
        xframe.Draw()
        c1.cd(2)
        xframe2.Draw()
        c1.SaveAs(fitm+"/fit"+histos[h]+".pdf")
    #for i in [5,0,1,4,3,2]:
        sigfile.write(str(bin)+'\t'+str(signal3900.getValV())+'\t'+str(0)+'\t'+str(0)+'\t'+str(-signal3900.getErrorLo())+'\t'+str(signal3900.getErrorHi())+'\t')
        #sigfile.write(str(xval[i])+'\t'+str(hval[i])+'\t'+str(0)+'\t'+str(0)+'\t'+str(-herrLo[i])+'\t'+str(herrHi[i])+'\t')
        sigfile.write('\n')
    sigfile.close()
    sigfiler = open(fitm+"/signalval.txt",'r')
    line = sigfiler.readline()
    while line:
        print line
        line = sigfiler.readline()
    sigfiler.close()
    c2=TCanvas("c2","A Simple Graph with asymmetric error bars",200,10,700,500)
    #graph = rt.TGraphAsymmErrors(6,xval,hval,xerrLo,xerrHi,herrLo,herrHi)
    #graph = rt.TGraphAsymmErrors(fitm+"/signalval.txt", "%lg %lg %lg %lg %lg %lg")
    #print graph.GetErrorYlow(3)
    #graph.SetTitle("Z_{c}(3900) signal events")
    #graph.SetMarkerColor(1)
    #graph.SetMarkerStyle(20)
    #graph.GetXaxis().SetTitle("M_{J/#psi#pi#pi} (GeV)");
    #graph.GetYaxis().SetTitle("Events /50 MeV");
    #graph.GetYaxis().SetRangeUser(0.0,300);
    #graph.Draw("AP")
    #c2.SaveAs(fitm+"/signal3900.pdf")
if __name__ == '__main__':
    path = sys.argv[2]
    for arg in sys.argv:  
        print arg
    main()

