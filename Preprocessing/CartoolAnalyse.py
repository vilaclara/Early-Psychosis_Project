import numpy as np
import CartoolFiles
class Interpolation:
    def __init__(self,FileXyz,FileEph,BadElectrodesName,eph_file=True,writing=True):
        self.writing=writing
        self.XyzData=CartoolFiles.Xyz(FileXyz)
        if eph_file:
            self.EphData=CartoolFiles.Eph(FileEph)
        else:
            self.EphData=CartoolFiles.Eph('')
            self.EphData.Data=FileEph
            
        # Extraself.EphDataction des numeros des electrodes
        self.InterpolateElectrodes=[]
        self.CorrectElectrode=range(self.XyzData.NbElectrodes)
        self.BadElectrodes=BadElectrodesName
        for i in BadElectrodesName:
            if len(i)>1:
                try:
                    self.InterpolateElectrodes.append(self.XyzData.ElectrodeName.index(i+'\r'))
                    self.CorrectElectrode.remove(self.XyzData.ElectrodeName.index(i+'\r'))
                except ValueError:
                    self.InterpolateElectrodes.append(self.XyzData.ElectrodeName.index(i))
                    self.CorrectElectrode.remove(self.XyzData.ElectrodeName.index(i))
        self.FileName=FileEph
    def CalclM2(self,Name):
        # Calcul des coefficiants pour interpolation m=2 comme carto
        Potential=np.concatenate((self.EphData.Data[:,self.CorrectElectrode],np.zeros((self.EphData.Data.shape[0],4))),axis=1)
        MatriceE=np.array([])
        for c in self.XyzData.Coord[self.CorrectElectrode,:]:
            tmp=np.array(1)
            tmp=np.append(tmp,c)
            MatriceE=np.append(MatriceE,tmp,axis=0)
        MatriceE=MatriceE.reshape((MatriceE.shape[0]/4,4))
        MatriceK=np.zeros((len(self.CorrectElectrode),len(self.CorrectElectrode)))
        for i,r1 in enumerate(self.XyzData.Coord[self.CorrectElectrode,:]):
            for j,r2 in enumerate(self.XyzData.Coord[self.CorrectElectrode,:]):
                if i==j:
                    MatriceK[i,j]=0
                else:
                    Diff=(np.square(r1-r2)).sum()
                    MatriceK[i,j]=Diff*np.log(Diff)
        Mat=np.concatenate((MatriceK,MatriceE),axis=1)            
        tmp=np.concatenate((MatriceE.T,np.zeros((4,4))),axis=1)
        Mat=np.concatenate((Mat,tmp),axis=0)
        MatInv=np.linalg.inv(Mat)
        Coef=[]
        for v in Potential:
            Coef.append(np.dot(MatInv,v))
        Coef=np.array(Coef)
        # calcul interpolation des bad Chanel
        NewData=self.EphData.Data
        for b in self.InterpolateElectrodes:
            r=self.XyzData.Coord[b].reshape((1,3))
            Q=np.array(1)
            Q=np.append(Q,r)
            CorrectCoord=self.XyzData.Coord[self.CorrectElectrode]
            r=r.repeat(CorrectCoord.shape[0],axis=0)
            Diff=(np.square(r-CorrectCoord)).sum(axis=1)
            K=Diff*np.log(Diff)
            K[np.isnan(K)]=0
            IntData=np.dot(Coef,np.concatenate((K,Q), axis=0))
            NewData[:,b]=IntData
        # Ecrire New Chanel
        self.EphData.Data=NewData
        if self.writing:
            self.EphData.write(self.FileName.replace('.eph',"".join(['.',Name,'.eph'])))
        else:
            self.InterpoletedData=NewData #.reshape((self.EphData.TF,self.EphData.Electrodes))
