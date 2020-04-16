clc
clear all
close all

task='auditory';

OutDir='/Users/mip/Documents/PdM/Data/ERPs/Dataset1/';
Bad_elec_file=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Elec2Interpolate_Database1.xls');
[~,~,raw] = xlsread(Bad_elec_file);
Subj_names=raw(:,1);
Bad_elec=raw(:,2);

for subj=1:size(Subj_names,1)
    Subj=Subj_names{subj};
    if strcmp(task,'auditory')
        trig_type_name{1}='std';
        trig_type_name{2}='dev1';
        trig_type_name{3}='dev2';
        trig_type_name{4}='dev3';
    end

    for trig=1:size(trig_type_name,2)
        OutDir1=sprintf('%s/%s/epochs',OutDir,trig_type_name{trig});
        OutDir2=sprintf('%s/%s/epochs/%s',OutDir,trig_type_name{trig},Subj);
        
        load(sprintf('%s/nb_epochs.mat',OutDir2));
        count_epoch(subj,trig)=count_e;
        count_rej(subj,trig)=count_r;
        count_ratio(subj,trig)=count_e/(count_e+count_r);
    end
    
end

Subj2check1=Subj_names(sum(count_ratio<0.5,2)>=1,:)
Subj2check2=Subj_names(sum(count_epoch(:,2:4),2)<100,:)
