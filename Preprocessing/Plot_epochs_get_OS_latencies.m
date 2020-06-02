clc
clear all
close all

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
trig='std';
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
    
    OutDir1=sprintf('%s%s/epochs',OutDir,trig);
    
    n_poc_pt=1;
    n_poc_ct=1;
    
    Pt_pocs=[];
    Ct_pocs=[];
    
    for subj=1:size(OS_groups{os},1)
    
        Subj=OS_groups{os}{subj};
        if Subj(1)=='F'
            pt=1;
        else
            pt=0;
        end
        Npoc=1;
        
        OutDir2=sprintf('%s/%s/epochs/%s',OutDir,trig,Subj);

        try
            load(sprintf('%s/nb_epochs.mat',OutDir2));
            count_epoch=count_e-1;
            count_rej=count_r-1;
        catch
            count_epoch=0;
            count_rej=0;
        end

        if count_epoch
            files=dir(sprintf('%s/%s*.mat',OutDir2,Subj));
            
            if pt
                for poc=Npoc:size(files,1)
                    load(sprintf('%s/%s',OutDir2,files(poc-Npoc+1).name));
                    Pt_pocs(:,n_poc_pt)=double(epoch(elec,:)-mean(epoch(elec,100:baseline_tp),2));%-mean(epoch(elec,100:baseline_tp),2)
                    Npoc=poc+1;
                    n_poc_pt=n_poc_pt+1;
                end
            else
                for poc=Npoc:size(files,1)
                    load(sprintf('%s/%s',OutDir2,files(poc-Npoc+1).name));
                    Ct_pocs(:,n_poc_ct)=double(epoch(elec,:)-mean(epoch(elec,100:baseline_tp),2));%-mean(epoch(elec,100:baseline_tp),2)
                    Npoc=poc+1;
                    n_poc_ct=n_poc_ct+1;
                end
            end
        end
        a=1;
        clear count_epoch count_rej
    end
    
    Avg_Pt=mean(Pt_pocs,2);
    Avg_Ct=mean(Ct_pocs,2);
    Avg_All=mean([Pt_pocs Ct_pocs],2);
    
    save(sprintf('%s/OS_%d_elec_%s_epoch_avg.mat',OutDir0,os,chan_names{elec}),'Avg_Pt','Avg_Ct','Avg_All');
    
    fig=figure;
    plot(Avg_Ct,'LineWidth',2)
    hold on
    plot(Avg_Pt,'LineWidth',2)
    hold on
    plot(Avg_All,'LineWidth',2)
    l=legend('Ct','Pt','All');
    l.FontSize=16;
    text1=sprintf('OS %d elec %s (N=%d)',os,chan_names{elec},size(OS_groups{os},1));
    t=title(text1);
    t.FontSize=18;
    xticks([0 204 306.4 408.8])
    xticklabels({'-200','0','100','200'})
    saveas(fig,sprintf('%s/OS_%d_elec_%s_plot.png',OutDir0,os,chan_names{elec}));
    a=1;
    clear Pt_pocs Ct_pocs n_poc_pt n_poc_ct Npoc
end

%% plot

figure
for os=1:size(OS_groups,2)
    OutDir0=sprintf('%s%s',OutDir,trig);
    load(sprintf('%s/OS_%d_elec_%s_ERP.mat',OutDir0,os,chan_names{elec}),'Avg_Pt','Avg_Ct','Avg_All');
    plot(Avg_Ct,'LineWidth',2)
    hold on
    a=1; 
end
xticks([0 204 306.4 408.8])
xticklabels({'-200','0','100','200'})
xline(baseline_tp);
xline(baseline_tp+51.2);
xline(baseline_tp+102.4);
l=legend('xp','win7','linux','cery win7 1','cery win7 2','baseline','50','100');
l.FontSize=16;

a=1;

