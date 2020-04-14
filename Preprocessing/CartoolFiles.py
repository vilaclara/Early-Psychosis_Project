import numpy as np
import tables
import struct
import os
import re
import shutil
class Lffile:
    def __init__(self,FileName):
        with open(FileName,'rb') as fid:
            self.numelec=np.fromfile(fid,dtype=np.int32,count=1)
            self.numsolutionpoints=np.fromfile(fid,dtype=np.int32,count=1)/3
            self.Data=np.fromfile(fid,dtype=np.float64)
            self.Data=self.Data.reshape((self.numelec,self.numsolutionpoints,3))
class Isfile:
    def __init__(self,FileName):
        with open(FileName,'rb') as fid:
            magic=fid.read(4)
            self.numelec=np.fromfile(fid,dtype=np.int32,count=1)
            self.numsolutionpoints=np.fromfile(fid,dtype=np.int32,count=1)
            self.numreg=np.fromfile(fid,dtype=np.int32,count=1)
            self.isinversescalar =np.fromfile(fid,dtype=np.bool,count=1)
            self.elecnames=[fid.read(32).replace('\x00','') for e in range(self.numelec)]
            self.spinames=[fid.read(16).replace('\x00','') for e in range(self.numsolutionpoints)]
            self.regvalues =np.fromfile(fid,dtype=np.float64,count=self.numreg)
            self.regnames  = [fid.read(32).replace('\x00','') for e in range(self.numreg)]
            self.data=np.fromfile(fid,dtype=np.float32)
            self.data=self.data.reshape((self.numreg,self.numsolutionpoints,3,self.numelec))
class Vmrk:
    def __init__(self,FileName):
        with open(FileName,'r') as VmrkData:
            Stim=[]
            Header=[]
            for r in VmrkData.readlines():
                if r.find('Stimulus')!=-1:
                    Trial=r.split(',')
                    TrigNumber=int(Trial[1].split(" ")[-1])
                    TimeCode=int(Trial[2])
                    Stim.append([TrigNumber,TimeCode])
                    Typeline=r
                else:
                    Header.append(r)
            self.Stim=Stim
            self.Header=Header
    def write(self,FileName):
        with open(FileName,'w') as ResData:
            for h in self.Header:
                ResData.write(h)
            for i,s in enumerate(self.Stim):
                ResData.write("".join(['Mk',str(i+2),'=Stimulus,S ',str(np.absolute(s[0])),',',str(s[1]),',1,0\n']))
class Tva:
    def __init__(self,FileName):
        if FileName=='':
            self.Data=[]
        else:
            Tva= open(FileName, "r")
            Header=Tva.readline()
            Data=[]
            for ligne in Tva:
                ligne=ligne.replace('\n','')
                Data.append(ligne.split('\t'))
            self.Data=Data
    def write(self,FileName):
        fichier=open(FileName,"w")
        fichier.write('TV01')
        fichier.write('\n')
        for l in self.Data:
            #ecriture ligne par ligne
            fichier.write("\t".join(l))
            #saut de ligne
            fichier.write('\n')
        fichier.close()
class Mrk:
    def __init__(self,FileName):
        if FileName=='':
            self.Time=np.array([])
            self.Trig=[]
        else:
            MrkData=open(FileName,'r')
            Head=MrkData.readline()
            Data=MrkData.readlines()
            Time=[]
            Trig=[]
            for d in Data:
                info=d.split('\t')
                Time.append(int(info[1]))
                Trig.append(info[-1].replace('\n',''))
            Time=np.array(Time)
            MrkData.close()
            self.Time=Time
            self.Trig=Trig
    def write(self,FileName):
        MrkFile=open(FileName,'w')
        MrkFile.write('TL02')
        MrkFile.write('\n')
        for i,v in enumerate(self.Trig):
            MrkFile.write('\t')
            MrkFile.write(str(self.Time[i]))
            MrkFile.write('\t')
            MrkFile.write(str(self.Time[i]))
            MrkFile.write('\t')
            MrkFile.write("".join(['"',str(v),'"']))
            MrkFile.write('\n')
        MrkFile.close()
class Vrb:
    def __init__(self,FileName):
        VrbFile= open(FileName, "r")
        TriggerValue=[]
        TriggerNumber=[]
        for ligne in VrbFile:
            if ligne.find("Trigger:")!=-1:
                Value=ligne.split(" ")[-1]
                if Value!="*\n":
                    TriggerValue.append(int(Value))
            if ligne.find("Total number of triggers:")!=-1:
                Number=ligne.split(" ")[-1]
                if Number!="*\n":
                    TriggerNumber.append(int(Number))
        self.Value=np.array(TriggerValue)
        self.Number=np.array(TriggerNumber)
        VrbFile.close()

class Ep: #lire les Ep, puis faire les claluls dessus (GFP, ST, ..:)
    """lecture: path eph/eph name""" 
    def __init__(self,FileName):
        """ on initialise l'objet ep
        soit on lis des EP 1 parametres 
        1) FileName = Eph full path if '' create an obj only
        """
        if FileName=='':
            self.Data=np.array([])
        else:
            self.Data=np.loadtxt(FileName)
    def write(self,FileName): # ecrire l'objet 
        """ ecritrue de l'eph : 
        nom de l'ep = name_eph  
        path de l'ep= path_result"""
        #creation du fichier
        fichier=open(FileName,"w")
        # boucle sur les time chaque ligne est un temps
        for time in self.Data:
            #ecriture ligne par ligne
            time.tofile(fichier,sep='\t',format="%s")
            #saut de ligne
            fichier.write('\n')
        fichier.close()


class Eph: #lire les Eph, puis faire les claluls dessus (GFP, ST, ..:)
    """lecture: path eph/eph name""" 
    def __init__(self,FileName):
        """ on initialise l'objet eph
        soit on lis des EPH 1 parametres 
        1) FileName = Eph full path if '' create an obj only
        """
        if FileName=='':
            self.TF=1
            self.Electrodes=1
            self.Fs=1
            self.Data=np.array([])
            self.GFP=np.array([])
        else:
            header = open(FileName).readline()
            header=header.split('\t')
            self.Electrodes =int(header[0])
            self.TF=int(header[1])
            self.Fs=int(header[2])
            self.Data=np.loadtxt(FileName,skiprows=1)
            self.Data=self.Data.reshape((self.TF,self.Electrodes))
            """ calcul le gfp de l'objet et revois sa valeur"""
            try:
                GFP=self.Data.std(1)
                self.GFP=GFP
            except:
                self.GFP=False
    def write(self,FileName): # ecrire l'objet 
        """ ecritrue de l'eph : 
        nom de l'eph = name_eph  
        path de l'eph = path_result"""
        #creation du fichier
        fichier=open(FileName,"w")
        # on prend le header
        header=[str(self.Electrodes),'\t',str(self.TF),'\t',str(self.Fs),'\n']
        #ecrtiture du header
        fichier.write("".join(header))
        # boucle sur les time chaque ligne est un temps
        for time in self.Data:
            #ecriture ligne par ligne
            time.tofile(fichier,sep='\t',format="%s")
            #saut de ligne
            fichier.write('\n')
        fichier.close()
class Xyz:
    def __init__(self,FileName):
        if FileName=='':
            self.Coord=np.array([])
            self.ElectrodeName=''
        else:
            XyzData=open(FileName,'r')
            self.NbElectrodes=int(XyzData.readline())
            Coord=[]
            ElectrodeName=[]
            for line in XyzData.readlines():
                line=line.split('\t')
                if len(line)==4:
                    if line[-1][-1]=='\n':
                        ElectrodeName.append(line[-1][0:-1].upper())
                    else:
                        ElectrodeName.append(line[-1].upper())
                    tmp=[]
                    for l in line[0:-1]:
                        tmp.append(float(l))
                    Coord.append(tmp)
            self.Coord=np.array(Coord)
            self.ElectrodeName=ElectrodeName
    def write(self,FileName):
        fichier=open(FileName,"w")
        # on prend le header
        header=str(len(self.ElectrodeName))
        #ecrtiture du header
        fichier.write("".join(header))
        fichier.write('\n')
        # boucle sur les time chaque ligne est un temps
        for i,c in enumerate(self.Coord):
            #ecriture ligne par ligne
            c.tofile(fichier,sep='\t',format="%s")
            fichier.write('\t')
            fichier.write(self.ElectrodeName[i])
            #saut de ligne
            fichier.write('\n')
        fichier.close()
        
class Bdf:
    def __init__(self,BdfFile,InfoOnly=False):
        H5Name=BdfFile.replace('.bdf','.h5')
        if os.path.isfile(H5Name):
            OutPutFile=tables.open_file(BdfFile.replace('.bdf','.h5'),'r')
            ChanelName=OutPutFile.get_node('/ChanelName').read()
            self.BdfFile=BdfFile
        else:
            
            f = open(BdfFile, "rb")
            self.BdfFile=BdfFile
          
            # first 236 bytes information not usefull
            f.read(168)
            Date=f.read(8)
            Hour=f.read(8)
            f.read(52)
            #Data Record
            DataRec=int(f.read(8))
            # Duration of Data record
            Duration=int(f.read(8))
            # number of recorder chanel
            NbChanel=f.read(4)
            NbChanel=int(NbChanel)
           
            # name of each recorded chanel
            ChanelName=[]
            for i in enumerate(range(NbChanel)):
                    ChanelName.append(f.read(16))
                    
            #try:
            #    ChanelName=[x.replace(' ','') for x in ChanelName]
            #except TypeError:
            #    print('Weird error in cartool file Bdf function about ChanelName')
            ChanelName=[x.replace(' ','') for x in ChanelName]
            # Tansducteur type not use full
            f.read(NbChanel*80)
            # Physical dimention not usefull
            f.read(NbChanel*8)
            #physical min and max not usefull
            f.read(NbChanel*16)
            # Digital min and max not usefull
            f.read(NbChanel*16)
            # prefiltering not usefull
            f.read(NbChanel*80)
            # sampling rate
            FS=[]
            for i in enumerate(range(NbChanel)):
                    FS.append(int(f.read(8)))
            
            # Reserve not usefull
            f.read(NbChanel*32)
            # RawData
            if InfoOnly:
                self.Info=[Date,Hour]
                #print('Date')
                print(Date)
                #print('Hour')
                #print(Hour)
            else:
                # h5file with same name and location as bdf file
                OutPutFile=tables.open_file(BdfFile.replace('.bdf','.h5'),'w')  
                OutPutFile.create_array('/','NbChanel',NbChanel)
                OutPutFile.create_array('/','ChanelName',ChanelName)
                OutPutFile.create_array('/','FS',FS)
                RawData=OutPutFile.create_earray('/','RawData',tables.Float64Atom(),(NbChanel,0))
                AvgRef=OutPutFile.create_earray('/','AvgRef',tables.Float64Atom(),(NbChanel,0))
                for d in range(DataRec):
                        Data=[]
                        for c in range(NbChanel):
                                Chanel=[]
                                for p in range(FS[c]):
                                    Bin=(f.read(3))
                                    Chanel.append(struct.unpack('<l','\xc3'+Bin)[0]/float(256*32))
                                Data.append(np.array(Chanel))
                        Data=np.array(Data)
                        Avg=Data.mean(axis=0)
                        Avg=Avg.reshape((1,Avg.shape[0]))
                        Avg=Avg.repeat(Data.shape[0],axis=0)
                        RawData.append(Data)
                        AvgRef.append(Data-Avg)
                trig=ChanelName.index('Status')
                if trig!=[]:
                    self.Statut=OutPutFile.get_node('/RawData')[trig,:]
                OutPutFile.close()
    def ExtractMrk(self,WriteFile=True):
        Time=(np.diff(np.sign(np.diff(self.Statut)))<0).nonzero()[0]+1
        Trig=(self.Statut[Time+1]-self.Statut[Time-1])/0.03125
##        if (Trig[1::2]<0).all()==False:
##            Error=np.nonzero(Trig[1::2]<0)[0]
##            Trig[(2*Error)+1]=-Trig[(2*Error)+2]
        Trig=Trig[1::2]
        Time=Time[1::2]
        # reference value for 1 =0.03125
        Trig=np.int32(Trig/0.03125)
        # write mrk file
        if WriteFile:
            MrkFile=Mrk('')
            MrkFile.Time=Time
            MrkFile.Trig=Trig
            MrkFile.write("".join([self.BdfFile,'.mrk']))
        self.MrkTime=Time
        self.MrkTrig=Trig
    def RemoveH5(self):
        os.remove(self.BdfFile.replace('.bdf','.h5'))
        
class BrainVoyager:
     def __init__(self,DatFile):
        with open(DatFile,'r') as RawData:
            Data=[]
            ElectrodeName=[]
            for line in RawData.readlines():
                 line=line.split(' ')
                 line=np.array(line)
                 line=line[line!='']
                 ElectrodeName.append(line[0])
                 tmp=line[1:]
                 Data.append(np.float32(tmp))
            Data=np.array(Data)
        self.Data=Data
        self.ElectrodeName=ElectrodeName
             
class Spi:
    def __init__(self,SpiFile):
        if SpiFile=='':
            self.Label=['']
            self.Coordonee=np.array([])
        else:
            Data=np.loadtxt(SpiFile,dtype='string')
            if Data.shape[1]>3:
                self.Label=Data[:,3]
            else:
                self.Label=['']
            self.Coordonee=np.float64(Data[:,0:3])
         
    def write(self,FileName):
        fid=open(FileName)
        for i,c in enumerate(self.Coordonee):
            c.tofile(fid,sep='\t')
            fid.write('\t'+self.Label+'\n')
        fid.close()
        
    def MatrixDist(self):
        NbPoint=self.Coordonee.shape[0]
        MatrixDist=np.zeros((NbPoint,NbPoint))
        for v in range(NbPoint):
            dist=self.Coordonee-self.Coordonee[v,:]
            dist=dist*dist
            dist=dist.sum(1)
            dist=np.sqrt(dist)
            MatrixDist[v,:]=dist
        self.Distance=MatrixDist
class BrainVision:
##    Read .eeg file, with this file (in the smme folder
##    it must have a .vmrk and .vhdr with the same name as the .eeg file
##    the Vhdr file is unsefull for the file information
##    the Vmrk File is use full for event(trigger)   
    def __init__(self,FileDotEEG):
        self.FileName=FileDotEEG
        self.Info=self.__ReadVhdrFile__()
        self.__ReadVmrkFile__()
    def Read(self):
        Format = { 'INT_16':np.int16, 'IEEE_FLOAT_32':np.float32}
        DataType=Format[self.Info['Format']]
        # Data read memap Format
        Data = np.memmap(self.FileName, DataType, 'r')
        self.NbTimeFrame = int(Data.size/self.Info['NumberOfChannels'])
        Data=Data.reshape(self.NbTimeFrame, self.Info['NumberOfChannels'])
        self.MemapData=Data
    def AvgCalculation(self,Track):
        if type(Track)==int:
            Avg=self.MemapData[:,Track]
        else:
            InterestData=self.MemapData[:,Track]
            Avg=InterestData.mean(axis=1)
        Avg=Avg.reshape((len(Avg),1))
        Avg=Avg.repeat(self.MemapData.shape[1],axis=1)
        self.AvgMemapData=self.MemapData-Avg
    def RemoveOneChanel(self,ChannelValue,Name):
        NbChanel=self.Info['NumberOfChannels']
        el=range(NbChanel-1)
        el.remove(ChannelValue)
        self.MemapData=self.MemapData[:,el]
        self.__ChannelNameRemovedName__=self.Info['ChanelName'][ChannelValue]
        self.WriteDotEEG(Name,WritingVhdr=True)
    def WriteDotEEG(self,Name,WritingVmrk=False,WritingVhdr=False):
        Format = { 'INT_16':np.int16, 'IEEE_FLOAT_32':np.float32}
        DataType=Format[self.Info['Format']]
        print(self.MemapData.shape)
        Data = np.memmap(self.FileName.replace('.eeg','.'+Name+'.eeg'), DataType, 'w+',shape=self.MemapData.size)
        Data[:]=self.MemapData.reshape(np.array(self.MemapData.size).prod())
        del Data
        if WritingVhdr:
            self.__WriteWhdr__(self.FileName.replace('.eeg','.'+Name+'.vhdr'),self.__ChannelNameRemovedName__)
        else:
            src=self.FileName.replace('.eeg','.vhdr')
            dst=self.FileName.replace('.eeg','.'+Name+'.vhdr')
            shutil.copyfile(src,dst)
        if WritingVmrk:
            self.WriteVmrk(self.FileName.replace('.eeg','.'+Name+'.vmrk'))
        else:
            src=self.FileName.replace('.eeg','.vmrk')
            dst=self.FileName.replace('.eeg','.'+Name+'.vmrk')
            shutil.copyfile(src,dst)

    def __WriteWhdr__(self,Name,ChannelNameRemoved):
        Impedence=False
        ChanelName=self.Info['ChanelName']
        ChanelName.remove(ChannelNameRemoved)
        with open(Name,'w') as f:
            c=0
            for l in self.__RawVhdr__:
                if l.find('\xc2\xb5V')!=-1:
                    if c<len(ChanelName):
                        if ChanelName[c].find('Aux')!=-1:
                            Text='Ch'+str(c+1)+'='+ChanelName[c]+',,1.0,\xc2\xb5V\n'
                        else:
                            Text='Ch'+str(c+1)+'='+ChanelName[c]+',REF,1.0,\xc2\xb5V\n'
                        f.write(Text)
                        c+=1
                elif l.find('Reference channel:')!=-1:
                    f.write('Reference channel: '+ChannelNameRemoved+'\n')
                elif l.find('Impedance [KOhm]')!=-1:
                    Impedence=True
                elif Impedence:
                    f.write(l.replace(ChannelNameRemoved+':','REF_'+ChannelNameRemoved+':'))
                else:
                    f.write(l)
    def __ReadVhdrFile__(self):
        FileName=self.FileName
        section = None
        all_info = { }
        Lines=[]
        for line in open(FileName.replace('.eeg','.vhdr'), 'rU'):
            Lines.append(line)
            line = line.strip('\n').strip('\r')
            if line.startswith('['):
                section = re.findall('\[([\S ]+)\]', line)[0]
                all_info[section] = { }
                continue
            if line.startswith(';'):
                continue
            if '=' in line and len(line.split('=')) ==2:
                k,v = line.split('=')
                all_info[section][k] = v
        self.__RawVhdr__=Lines
        RelevantInfo={}
        RelevantInfo['NumberOfChannels']=int(all_info['Common Infos']['NumberOfChannels'])
        print(all_info['Common Infos']['SamplingInterval'])
        RelevantInfo['Fs']=1000000./int(all_info['Common Infos']['SamplingInterval'])
        RelevantInfo['Format']=all_info['Binary Infos']['BinaryFormat']
        ChanelName=[]
        for c in range(RelevantInfo['NumberOfChannels']):
            ChanelName.append(all_info['Channel Infos']['Ch%d' % (c+1,)].split(',')[0])
        RelevantInfo['ChanelName']=ChanelName
        return RelevantInfo

    def __ReadVmrkFile__(self):
        FileName=self.FileName
        with open(FileName.replace('.eeg','.vmrk'),'r') as VmrkData:
            Event=[]
            Header=[]
            for r in VmrkData.readlines():
                if r.find('Stimulus')!=-1:
                    Trial=r.split(',')
                    TrigNumber=Trial[1]
                    TimeCode=int(Trial[2])
                    Event.append([TrigNumber,TimeCode])
                    Typeline=r
                elif r.find('Trigger')!=-1:
                    Trial=r.split(',')
                    TrigNumber=Trial[1]
                    TimeCode=int(Trial[2])
                    Event.append([TrigNumber,TimeCode])
                    Typeline=r
                else:
                    Header.append(r)
        self.__HeaderVmrk__= Header
        self.Event=Event
    def WriteVmrk(self,FileName):
        with open(FileName,'w') as f:
            for h in self.__HeaderVmrk__:
                f.write(h)
            for i,s in enumerate(self.Event):
                if type(s[0])!='float':
                    f.write("".join(['Mk',str(i+2),'=Stimulus,',str(s[0]),',',str(s[1]),',1,0\n']))
                else:
                    f.write("".join(['Mk',str(i+2),'=Stimulus,S ',str(np.absolute(s[0])),',',str(s[1]),',1,0\n']))


        
