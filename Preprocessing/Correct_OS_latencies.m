clc
clear all
close all

% Manual computed latencies:
OS_lat{1}=0; %xp - nestle
OS_lat{2}=-15; %win7 - nestle
OS_lat{3}=-25; %linux - nestle
OS_lat{4}=26; %win7 - Cery ancien labo
OS_lat{5}=22; %win7 - Cery nouveau labo


% db=1;
% load('OS_groups_db1.mat','Subj_names')
% OS_groups{1}={}; %xp - nestle
% OS_groups{2}={}; %win7 - nestle
% OS_groups{3}={}; %linux - nestle
%Tous les sujets NAC ont été enregistré à Nestlé, mais check les dates d'enregistrement eeg pour savoir si c'était sous Windows XP ou Win7 (2 groupes différents car certainement pas les mêmes latences...)
%Labo EEG Hôpital Nestlé (data JFK, Chrysa ; Dir. M.Murray):
%Windows XP : datasets enregistrés avant ~01.01.2013
%Win7 : datasets enregistrés dès ~01.01.2013 jusqu'à ~01.03.2017
%Linux : datasets enregistrés dès ~01.03.2017

db=2;
load('OS_groups_db2.mat','Subj_names')
OS_groups{1}={}; %xp - nestle
OS_groups{2}={}; %win7 - nestle
OS_groups{3}={}; %linux - nestle
OS_groups{4}={}; %win7 - Cery ancien labo
OS_groups{5}={}; %win7 - Cery nouveau labo
% Pour les sujets NCCR-Biomarqueurs enregistrés à Nestlé (JFK ou Chrysa):
% Les datasets nommés "S1000_LS_...bdf" ont été enregistrés à Nestlé sous Windows XP (les no à partir de S1000, S1001, S1002, etc.)
% Les datasets nommés "S2000_Ls_...bdf" ont été enregistrés à Nestlé sous Windows 7 ou sous Linux: tu dois regarder la date enregistrement pour déterminer cela
% Win7 : datasets enregistrés avant le ~01.03.2017
% Linux : datasets enregistrés dès ~01.03.2017
% Pour les sujets NCCR-Biomarqueurs enregistrés à Cery (par moi):
% Les datasets nommés "S3000_Cery_...bdf" ont été enregistrés à Cery (Win7) dans l'ancien labo EEG (ancien CNP)
% Les datasets nommés "S4000_Cery_...bdf" ont été enregistrés à Cery (Win 7 toujours) dans le nouveau labo EEG (nouveau CNP)

for i=1:size(OS_groups,2)
    indexes=find(cell2mat(Subj_names(:,2))==i);
    OS_groups{i}=Subj_names(indexes,1);
end

%trig='std';
trig='dev';
% trig_type_name{2}='dev';
% trig_type_name{3}='IC';
% trig_type_name{4}='NIC';

baseline_tp=204;
OutDir=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Dataset%d/',db);
load('/Users/mip/Documents/PdM/Data/Cap/chanlocs.mat')
for chan=1:size(EEG.chanlocs,2)
    chan_names{chan}=EEG.chanlocs(chan).labels;
end   
elec=30; %30 for visual  % 32 CPz
%Cz(48) seems to have biggest peak

for os=1:size(OS_groups,2)
    
    OutDir0=sprintf('%s%s',OutDir,trig);
    
    Pt_ERPs=[];
    Ct_ERPs=[];
    
    Npt(os)=0;
    Nct(os)=0;
    
    for subj=1:size(OS_groups{os},1)
    
        Subj=OS_groups{os}{subj};
        
        try
            load(sprintf('%s/%s_ERP.mat',OutDir0,Subj));
            error_loading=0;
        catch
            error_loading=1;
        end

        if ~error_loading
            if Subj(1)=='F'
                Npt(os)=Npt(os)+1;
                Pt_ERPs(:,Npt(os))=ERP(elec,102-OS_lat{os}:102-OS_lat{os}+614);
            else
                Nct(os)=Nct(os)+1;
                Ct_ERPs(:,Nct(os))=ERP(elec,102-OS_lat{os}:102-OS_lat{os}+614);
            end
        end
        a=1;
        clear count_epoch count_rej
    end
    
    Avg_Pt=mean(Pt_ERPs,2);
    Avg_Ct=mean(Ct_ERPs,2);
    Avg_All=mean([Pt_ERPs Ct_ERPs],2);
    
    save(sprintf('%s/OS_%d_elec_%s_ERP_corr.mat',OutDir0,os,chan_names{elec}),'Avg_Pt','Avg_Ct','Avg_All');
    clear Pt_pocs Ct_pocs n_poc_pt n_poc_ct Npoc
end

%% plot

figure
for os=1:size(OS_groups,2)
    OutDir0=sprintf('%s%s',OutDir,trig);
    load(sprintf('%s/OS_%d_elec_%s_ERP_corr.mat',OutDir0,os,chan_names{elec}),'Avg_Pt','Avg_Ct','Avg_All');
    plot(Avg_Ct,'LineWidth',2)
    hold on
    a=1; 
end
xticks([0 102 204.8 307.2])
xticklabels({'-100','0','100','200'})
xline(baseline_tp-102);
xline(baseline_tp-102+51.2);
xline(baseline_tp-102+102.4);
l=legend('xp','win7','linux','cery win7 1','cery win7 2','baseline','50','100');
l.FontSize=16;

a=1;



