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
        BaselinePeriod=Value[:,e-204:e]
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
            TVA.append(0)
#            if (Diff<Threshold).sum() ==1:
#                print('ok')
#                TVA.append(1)
#            else:
#                TVA.append(0)
    TVA=np.array(TVA)
    InterestData=[]
    Accepted=0
    tmp=cf.Eph('')
    tmp.Electrodes=64
    tmp.TF=716
    tmp.Fs=1024
    for i,e in enumerate(Event):
        if TVA[i]==1:
            BaselinePeriod=Value[:,e-204:e]
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
    Order=4
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
    order=4
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

XlsFile=r'/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Elec2Interpolate_Database1.xls'
BdfFolder=r'/Users/mip/Documents/PdM/Data/BDFs/Database1_auditory'
TvaFolder=r'/Users/mip/Documents/PdM/Data/csv/tva'
#ResBeforeInt=r'/Users/laura/Documents/EPFL/Projets_Master/PdM/Data/EPH/BeforeInt'
ResAfterInt=r'/Users/mip/Documents/PdM/Data/ERPs/Dataset1'
# For interpoloation
xyzFile=r'/Users/mip/Documents/PdM/Data/Cap/64-Biosemi.xyz'

Sbj,Int=ReadXls(XlsFile)
Int=[i.upper() for i in Int]

for i,s in enumerate(Sbj):
    print(s)
    bdf=glob.glob('%s/%s*.bdf'%(BdfFolder,s))
    n=1
    acc=0
    Erp=[]
    for b in bdf:
        name=b.split('\\')[-1]
##        Tva=r'%s\%s'%(TvaFolder,name.replace('.bdf','.tva'))
##        Tva=np.int64(np.float64(np.array(cf.Tva(Tva).Data)))
        BdfData=cf.Bdf(b)
        try:
            BdfData.ExtractMrk(WriteFile=False)
        except AttributeError:
            print('Bdf instance has no attribute Statut')
            break
        Trig=np.unique(BdfData.MrkTrig)
        count=[]
        for t in Trig:
            count.append((BdfData.MrkTrig==t).sum())
        count=np.array(count)
        try:
            TrigStd=Trig[np.argmin(count)]
        except ValueError:
            print('Bdf problem:')
            print('b in Bdf:')
            print(b)
            print('count:')
            print(count)
            print('Trig:')
            print(Trig)
            print('TrigStd:')
            print(TrigStd)
        Events=BdfData.MrkTime[BdfData.MrkTrig==TrigStd]
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
        Sum,accepted,n=Averaging(80,Interpole.InterpoletedData.T,[204,614],Events,n,Res,s)
        Erp.append(Sum)
        acc+=accepted
        BdfData.RemoveH5()
        
    eph=glob.glob(ResBeforeInt+r'\%s\*.eph'%s)
    print('Number of accepted trials:')
    print(acc)
    # for e in eph:
    #     Interpole=ca.Interpolation(xyzFile,e,Int[i].split(' '))
    #     Interpole.CalclM2('Int')
    #     src=e.replace('.eph','.Int.eph')
    #     try:
    #         os.mkdir(ResAfterInt+'\%s'%s)
    #     except:
    #         pass
    #     dst=ResAfterInt+r'/%s/%s'%(s,e.split('\\')[-1])
    #     shutil.move(src,dst)
    # for b in bdf:
    #     src=b
    #     dst=BdfFolder+r'/Done/%s'%b.split('\\')[-1]
    #     shutil.move(src,dst)
    
##    Erp=np.array(Erp)
##    Erp=(Erp.sum(axis=0))/float(acc)
##    import matplotlib.pyplot as plt
##    plt.plot(Erp[0,:])
##    plt.show()


