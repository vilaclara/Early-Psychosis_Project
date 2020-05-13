clc
clear all
close all

db=2;

OutDir=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Dataset%d/',db);
Bad_elec_file=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Elec2Interpolate_Database%d.xls',db);
[~,~,raw] = xlsread(Bad_elec_file);
Subj_names=raw(:,1);
Bad_elec=raw(:,2);

for subj=1:size(Subj_names,1)
    Subj=Subj_names{subj};
    
    trig_type_name{1}='std';
    trig_type_name{2}='dev1';
    trig_type_name{3}='dev2';
    trig_type_name{4}='dev3';
    trig_type_name{5}='IC1';
    trig_type_name{6}='NIC1';
    trig_type_name{5}='IC2';
    trig_type_name{6}='NIC2';

    for trig=1:size(trig_type_name,2)
        OutDir1=sprintf('%s/%s/epochs',OutDir,trig_type_name{trig});
        OutDir2=sprintf('%s/%s/epochs/%s',OutDir,trig_type_name{trig},Subj);
        
        try
            load(sprintf('%s/nb_epochs.mat',OutDir2));
            count_epoch(subj,trig)=count_e-1;
            count_rej(subj,trig)=count_r-1;
            count_ratio(subj,trig)=(count_e-1)/(count_e+count_r-2);
        catch
            count_epoch(subj,trig)=0;
            count_rej(subj,trig)=0;
            count_ratio(subj,trig)=0;
        end
        
    end
    
end

Subj2check1=Subj_names(sum(count_ratio<0.3,2)>=1,:)

disp('Auditory \n')

Subj2check2=Subj_names(sum(count_epoch(:,2:4),2)<100,:)

disp('Visual \n')

Subj2check2=Subj_names(sum(count_epoch(:,5:6),2)<100,:)
