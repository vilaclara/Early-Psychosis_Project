clc
clear all
close all

db=2;

%trig='std';
trig='IC';
% trig_type_name{2}='dev';
% trig_type_name{3}='IC';
% trig_type_name{4}='NIC';

baseline_tp=102;
OutDir=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Dataset%d/',db);
load('/Users/mip/Documents/PdM/Data/Cap/chanlocs.mat')
for chan=1:size(EEG.chanlocs,2)
    chan_names{chan}=EEG.chanlocs(chan).labels;
end   
if strcmp(trig,'std') || strcmp(trig,'dev')
    elec=48;% 32 CPz Cz(48) seems to have biggest peak
else
    elec=30; %30 for visual  
end

Bad_elec_file=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Elec2Interpolate_Database%d.xls',db);
[~,~,raw] = xlsread(Bad_elec_file);
Subj_names=raw(:,1);

OutDir01=sprintf('%s%s/ERPs',OutDir,trig);
OutDir02=sprintf('%s%s/GFPs',OutDir,trig);
    
Pt_ERPs=[];
Ct_ERPs=[];

Npt=0;
Nct=0;

for subj=1:size(Subj_names,1)

    Subj=Subj_names{subj};

    try
        load(sprintf('%s/%s_ERP.mat',OutDir01,Subj));
        error_loading=0;
    catch
        error_loading=1;
    end

    if ~error_loading
        if Subj(1)=='F' || strcmp(Subj(1:2),'Ln')
            Npt=Npt+1;
            Pt_ERPs(:,Npt)=ERP(elec,:);
            Pt_names{:,Npt}=Subj;
        else
            Nct=Nct+1;
            Ct_ERPs(:,Nct)=ERP(elec,:);
            Ct_names{:,Nct}=Subj;
        end
    end
    a=1;
    clear count_epoch count_rej
end

Avg_Pt=mean(Pt_ERPs,2);
Avg_Ct=mean(Ct_ERPs,2);
All_ERPs=[Pt_ERPs Ct_ERPs];
All_names=[Pt_names Ct_names];
Avg_All=mean(All_ERPs,2);

%save(sprintf('%s/OS_%d_elec_%s_ERP.mat',OutDir01,os,chan_names{elec}),'Avg_Pt','Avg_Ct','Avg_All');

fig=figure;
plot(Avg_Ct,'LineWidth',2)
hold on
plot(Avg_Pt,'LineWidth',2)
hold on
plot(Avg_All,'LineWidth',2)
l=legend(sprintf('Ct (%d)',Nct),sprintf('Pt (%d)',Npt),sprintf('All (%d)',Nct+Npt));
l.FontSize=16;
text1=sprintf('Mean ERP %s elec %s db%d',trig,chan_names{elec},db);
t=title(text1);
t.FontSize=18;
xticks([0 102 204.4 306.6 409])
xticklabels({'-100','0','100','200','300'})
saveas(fig,sprintf('%s/Plot_mean_elec_%s_%s_plot.png',OutDir,trig,chan_names{elec}));
a=1;

%% plot

fig=figure;
for su=1:size(All_ERPs,2)
    plot(All_ERPs(:,su),'LineWidth',1)
    fprintf('\n %s\n',All_names{su})
    hold on
end
text1=sprintf('ERPs %s elec %s db%d',trig,chan_names{elec},db);
t=title(text1);
t.FontSize=18;
xticks([0 102 204.4 306.6 409])
xticklabels({'-100','0','100','200','300'})
saveas(fig,sprintf('%s/Plot_each_elec_%s_%s_plot.png',OutDir,trig,chan_names{elec}));
a=1;



