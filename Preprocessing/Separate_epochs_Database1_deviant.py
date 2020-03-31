# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import CartoolFiles as cf
import numpy as np
import glob
import os
import CartoolAnalyse as ca
import xlrd
from scipy import signal
import shutil
import tables


def Averaging(Threshold,InterpedData,EpochsTF,Event,n,Res,Name):
    # 2 step first using threshold to determine TVA
    # second use TVA for ERP
    Value=InterpedData.copy()
    TVA=[]
    for i,e in enumerate(Event):
        BaselinePeriod=Value[:,e-102:e]
        BaselineValue=np.mean(BaselinePeriod, axis=1) 
        InterestPeriod=Value[:,e-EpochsTF[0]:e+EpochsTF[1]]
        Bl_InterestPeriod=(InterestPeriod.transpose() - BaselineValue).transpose()
        Diff=np.array([np.abs(Bl_InterestPeriod.max(axis=1)),np.abs(Bl_InterestPeriod.min(axis=1))])
        Diff=Diff.max(axis=0)
##        avg=Diff[-1]
##        Diff=Diff[0:-1]
        if (Diff<Threshold).all():
            TVA.append(1)
        else:
            if (Diff<Threshold).sum() ==1:
                print('ok')
                TVA.append(1)
            else:
                TVA.append(0)
    TVA=np.array(TVA)
    InterestData=[]
    Accepted=0
    tmp=cf.Eph('')
    tmp.Electrodes=64
    tmp.TF=614
    tmp.Fs=1024
    for i,e in enumerate(Event):
        if TVA[i]==1:
            BaselinePeriod=Value[:,e-102:e]
            BaselineValue=np.mean(BaselinePeriod, axis=1) 
            InterestPeriod=Value[:,e-EpochsTF[0]:e+EpochsTF[1]]
            Bl_InterestPeriod=(InterestPeriod.transpose() - BaselineValue).transpose()
            Val=Bl_InterestPeriod
            tmp.Data=Val.T
            tmp.write(Res+r'/%s.%.4d.eph'%(Name,n))
            
            InterestData.append(Val)
            Accepted+=1
            n+=1
    InterestData=np.array(InterestData)
    ERP=InterestData.sum(axis=0)
    return ERP,Accepted,n

def BandPass(Frequencies2BandPass,Data,FS):
    Value=Data.copy()
    Order=2
    nyq=FS/2.
    w=Frequencies2BandPass/(nyq)
    FilteredData=np.zeros(Data.shape)
    b,a=signal.butter(Order,w,btype='bandpass')
    for i,d in enumerate(Value):
        FilteredData[i,:]=signal.filtfilt(b,a,d)
    return FilteredData
def RmvDC(Data):
    Data=Data.copy()
    Base=Data.mean(axis=1)
    Base=Base.reshape((len(Base),1))
    Base=Base.repeat(Data.shape[1],axis=1)
    Data=Data-Base
    return Data

def NotchFilter(Data,Freq,FS):
    Value=Data.copy()
    order=2
    nyq=FS/2.
    Freq = Freq / nyq
    b, a = signal.butter(order, [Freq - 1. / nyq, Freq + 1. / nyq],btype='bandstop')
    FilterData=np.zeros(Data.shape)
    for i,d in enumerate(Value):
        Val=signal.lfilter(b,a,d)
        FilterData[i,:]=signal.filtfilt(b,a,d)
    return FilterData
def ReadXls(XlsFile):
    wb=xlrd.open_workbook(XlsFile)
    sh=wb.sheet_by_index(0)
    Sbj=sh.col_values(0) 
    Int=sh.col_values(1)
    return Sbj,Int

XlsFile=r'/Users/laura/Documents/EPFL/Projets_Master/PdM/Data/Interpole-NewRaw.xls'
BdfFolder=r'/Volumes/Volume/EPFL/Projets/PdM/RawData'
#BdfFolder=r'/Users/laura/Documents/EPFL/Projets_Master/PdM/Data/RawData'
TvaFolder=r'/Users/laura/Documents/EPFL/Projets_Master/PdM/Data/csv/tva'
ResBeforeInt=r'/Users/laura/Documents/EPFL/Projets_Master/PdM/Data/EPH/BeforeInt'
ResAfterInt=r'/Volumes/Volume/EPFL/Projets/PdM/EPH/dev'
#ResAfterInt=r'/Users/laura/Documents/EPFL/Projets_Master/PdM/Data/EPH/dev'

# For interpoloation
xyzFile=r'/Users/laura/Documents/EPFL/Projets_Master/PdM/Data/Cap/64-Biosemi.xyz'

Sbj,Int=ReadXls(XlsFile)
Int=[i.upper() for i in Int]

for i,s in enumerate(Sbj):
    print(s)
    bdf=glob.glob('%s/%s*.bdf'%(BdfFolder,s))
    n=1
    try:
        os.mkdir(ResBeforeInt+'/%s'%s)
    except:
        pass
    Res=ResBeforeInt+'/%s'%s
    acc=0
    Erp=[]
    for b in bdf:
        name=b.split('\\')[-1]
##        Tva=r'%s\%s'%(TvaFolder,name.replace('.bdf','.tva'))
##        Tva=np.int64(np.float64(np.array(cf.Tva(Tva).Data)))
        BdfData=cf.Bdf(b)
        BdfData.ExtractMrk(WriteFile=False)
        Trig=np.unique(BdfData.MrkTrig)
        #count=[]
        #for t in Trig:
        #    count.append((BdfData.MrkTrig==t).sum())
        #count=np.array(count)
        #TrigStd=Trig[np.argmax(count)]
        Events=[]
        Events=np.concatenate([BdfData.MrkTime[BdfData.MrkTrig==64], BdfData.MrkTime[BdfData.MrkTrig==96]])
        Events=np.concatenate([Events,BdfData.MrkTime[BdfData.MrkTrig==128]])
        #Events=BdfData.MrkTime[BdfData.MrkTrig==TrigStd]
        H5=tables.open_file(BdfData.BdfFile.replace('.bdf','.h5'))
        FS=float(H5.get_node('/FS')[0])
        Data=H5.get_node('/RawData')[0:64,:]
        H5.close()
        FilteredData=RmvDC(Data)
        FilteredData1=BandPass(np.array([2,40]),FilteredData,FS)
        FilteredData2=NotchFilter(FilteredData1,50,FS)
        Interpole=ca.Interpolation(xyzFile,FilteredData2.T,Int[i].split(' '),eph_file=False,writing=False)
        Interpole.CalclM2('Int')
        
        try:
            os.mkdir(ResAfterInt+'/%s'%s)
        except:
            pass
        Res=ResAfterInt+r'/%s'%s
        # shutil.move(src,Res)
        Sum,accepted,n=Averaging(80,Interpole.InterpoletedData.T,[102,512],Events,n,Res,s)
        Erp.append(Sum)
        acc+=accepted
        BdfData.RemoveH5()
        
    eph=glob.glob(ResBeforeInt+r'\%s\*.eph'%s)
    print('Number of accepted trials:')
    print(acc)
    
