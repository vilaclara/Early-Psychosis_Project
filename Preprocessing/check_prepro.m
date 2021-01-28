clc
clear all
close all

CompSameSubj=1;


trigsa={'std','dev','IC','NIC'};
trigsb={'std','dev','visu_illu','visu_nn_illu'};
baseline_tp=102;

load('/Users/mip/Documents/PdM/Data/Cap/chanlocs.mat')
for chan=1:size(EEG.chanlocs,2)
    chan_names{chan}=EEG.chanlocs(chan).labels;
end   

for trig_i=3:4

    triga=trigsa{trig_i};
    trigb=trigsb{trig_i};
    
    if strcmp(triga,'std') || strcmp(triga,'dev')
        elec=48;% 32 CPz Cz(48) seems to have biggest peak
    else
        elec=26; %  POz(30) PO3(26) PO4(63) for visual 
    end

    % a = now
    % b = before
    % 1 = NAC
    % 2 = Synapsy
    
    Dira1=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Dataset1/');
    Dira2=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Dataset2/');
    Dirb1=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Nac/');
    Dirb2=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Synapsy/');

    MainDira1=sprintf('%s%s/ERPs',Dira1,triga);
    MainDira2=sprintf('%s%s/ERPs',Dira2,triga);
    MainDirb1=sprintf('%s%s_Avg_and_GFP',Dirb1,trigb);
    MainDirb2=sprintf('%s%s_Avg_and_GFP',Dirb2,trigb);
    
    files_pt_a1=dir(sprintf('%s/Lnac*_ERP.mat',MainDira1));
    files_ct_a1=dir(sprintf('%s/Ctrl*_ERP.mat',MainDira1));
    files_pt_b1=dir(sprintf('%s/Lnac*_int_CAR_mean.mat',MainDirb1));
    files_ct_b1=dir(sprintf('%s/Ctrl*_int_CAR_mean.mat',MainDirb1));
    files_pt_a2=dir(sprintf('%s/F*_ERP.mat',MainDira2));
    files_ct_a2=dir(sprintf('%s/L*_ERP.mat',MainDira2));
    files_pt_b2=dir(sprintf('%s/F*_int_CAR_mean.mat',MainDirb2));
    files_ct_b2=dir(sprintf('%s/L*_int_CAR_mean.mat',MainDirb2));
    
    %Patients
    Pt_nac_before=[];
    Pt_nac_before_names=[];
    for subj=1:size(files_pt_b1,1)
        TheName=extractBefore(files_pt_b1(subj).name,"_");
        Pt_nac_before_names=[Pt_nac_before_names [TheName,' ']];
        load(sprintf('%s/%s',files_pt_b1(subj).folder,files_pt_b1(subj).name));
        Pt_nac_before=[Pt_nac_before;mean_of_epochs(:,elec)'];
    end
    Pt_nac_now=[];
    Pt_nac_now_names=[];
    for subj=1:size(files_pt_a1,1)
        TheName=extractBefore(files_pt_a1(subj).name,"_");
        Pt_nac_now_names=[Pt_nac_now_names [TheName,' ']];
        load(sprintf('%s/%s',files_pt_a1(subj).folder,files_pt_a1(subj).name));
        Pt_nac_now=[Pt_nac_now;ERP(elec,2:end)];
    end
    Pt_sy_before=[];
    Pt_sy_before_names=[];
    for subj=1:size(files_pt_b2,1)
        TheName=extractBefore(files_pt_b2(subj).name,"_");
        Pt_sy_before_names=[Pt_sy_before_names [TheName,' ']];
        load(sprintf('%s/%s',files_pt_b2(subj).folder,files_pt_b2(subj).name));
        Pt_sy_before=[Pt_sy_before;mean_of_epochs(:,elec)'];
    end
    Pt_sy_now=[];
    Pt_sy_now_names=[];
    for subj=1:size(files_pt_a2,1)
        TheName=extractBefore(files_pt_a2(subj).name,"_");
        Pt_sy_now_names=[Pt_sy_now_names [TheName,' ']];
        load(sprintf('%s/%s',files_pt_a2(subj).folder,files_pt_a2(subj).name));
        Pt_sy_now=[Pt_sy_now;ERP(elec,2:end)];
    end
    
    %Controls
    Ct_nac_before=[];
    Ct_nac_before_names=[];
    for subj=1:size(files_ct_b1,1)
        TheName=extractBefore(files_ct_b1(subj).name,"_");
        Ct_nac_before_names=[Ct_nac_before_names [TheName,' ']];
        load(sprintf('%s/%s',files_ct_b1(subj).folder,files_ct_b1(subj).name));
        Ct_nac_before=[Ct_nac_before;mean_of_epochs(:,elec)'];
    end
    Ct_nac_now=[];
    Ct_nac_now_names=[];
    for subj=1:size(files_ct_a1,1)
        TheName=extractBefore(files_ct_a1(subj).name,"_");
        Ct_nac_now_names=[Ct_nac_now_names [TheName,' ']];
        load(sprintf('%s/%s',files_ct_a1(subj).folder,files_ct_a1(subj).name));
        Ct_nac_now=[Ct_nac_now;ERP(elec,2:end)];
    end
    Ct_sy_before=[];
    Ct_sy_before_names=[];
    for subj=1:size(files_ct_b2,1)
        TheName=extractBefore(files_ct_b2(subj).name,"_");
        Ct_sy_before_names=[Ct_sy_before_names [TheName,' ']];
        load(sprintf('%s/%s',files_ct_b2(subj).folder,files_ct_b2(subj).name));
        Ct_sy_before=[Ct_sy_before;mean_of_epochs(:,elec)'];
    end
    Ct_sy_now=[];
    Ct_sy_now_names=[];
    for subj=1:size(files_ct_a2,1)
        TheName=extractBefore(files_ct_a2(subj).name,"_");
        Ct_sy_now_names=[Ct_sy_now_names [TheName,' ']];
        load(sprintf('%s/%s',files_ct_a2(subj).folder,files_ct_a2(subj).name));
        Ct_sy_now=[Ct_sy_now;ERP(elec,2:end)];
    end
    
    %___________ find same subjects before and now
    
    Ct_sy_before_names=strsplit(Ct_sy_before_names);
    Ct_sy_now_names=strsplit(Ct_sy_now_names);
    for ind=1:size(Ct_sy_now_names,2)-1
        ind_ct_sy(1,ind)=ind;
        a = find(strcmp(Ct_sy_now_names{ind},Ct_sy_before_names)==1);
        if ~isempty(a)
            ind_ct_sy(2,ind)=a;
        else
            ind_ct_sy(2,ind)=0;
        end
    end
    
    Ct_nac_before_names=strsplit(Ct_nac_before_names);
    Ct_nac_now_names=strsplit(Ct_nac_now_names);
    for ind=1:size(Ct_nac_now_names,2)-1
        ind_ct_nac(1,ind)=ind;
        a = find(strcmp(Ct_nac_now_names{ind},Ct_nac_before_names)==1);
        if ~isempty(a)
            ind_ct_nac(2,ind)=a;
        else
            ind_ct_nac(2,ind)=0;
        end
    end
    
    Pt_sy_before_names=strsplit(Pt_sy_before_names);
    Pt_sy_now_names=strsplit(Pt_sy_now_names);
    for ind=1:size(Pt_sy_now_names,2)-1
        ind_pt_sy(1,ind)=ind;
        a = find(strcmp(Pt_sy_now_names{ind},Pt_sy_before_names)==1);
        if ~isempty(a)
            ind_pt_sy(2,ind)=a;
        else
            ind_pt_sy(2,ind)=0;
        end
    end
    
    Pt_nac_before_names=strsplit(Pt_nac_before_names);
    Pt_nac_now_names=strsplit(Pt_nac_now_names);
    for ind=1:size(Pt_nac_now_names,2)-1
        ind_pt_nac(1,ind)=ind;
        a = find(strcmp(Pt_nac_now_names{ind},Pt_nac_before_names)==1);
        if ~isempty(a)
            ind_pt_nac(2,ind)=a;
        else
            ind_pt_nac(2,ind)=0;
        end
    end
    
    ind_ct_sy=ind_ct_sy(:,find(ind_ct_sy(2,:)~=0));
    ind_pt_sy=ind_pt_sy(:,find(ind_pt_sy(2,:)~=0));
    ind_ct_nac=ind_ct_nac(:,find(ind_ct_nac(2,:)~=0));
    ind_pt_nac=ind_pt_nac(:,find(ind_pt_nac(2,:)~=0));
    
    % NAC
    
    mean_ct_nac_before=mean(Ct_nac_before,1);
    mean_pt_nac_before=mean(Pt_nac_before,1);
    mean_ct_nac_now=mean(Ct_nac_now,1);
    mean_pt_nac_now=mean(Pt_nac_now,1);

    % standard error = std/sqrt(N-1)

    stdE_ct_nac_before=std(Ct_nac_before,1)/sqrt(size(Ct_nac_before,1)-1);
    stdE_pt_nac_before=std(Pt_nac_before,1)/sqrt(size(Pt_nac_before,1)-1);
    stdE_ct_nac_now=std(Ct_nac_now,1)/sqrt(size(Ct_nac_now,1)-1);
    stdE_pt_nac_now=std(Pt_nac_now,1)/sqrt(size(Pt_nac_now,1)-1);
    
    fig=figure('Position',[10 10 1800 550]);
    filename=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Compare_old_new_prepro_nac_%s_%s.png',triga,chan_names{elec});

    subplot(2,1,1)
    a=mean_ct_nac_before + stdE_ct_nac_before;
    b=mean_ct_nac_before - stdE_ct_nac_before;
    dim=614;
    [asc, ord]=areapiena(a,b,dim);
    hp_ctrl_nac=patch(asc,ord, 'black');
    set(hp_ctrl_nac, 'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none');
    hold on
    a= mean_ct_nac_now + stdE_ct_nac_now;
    b=mean_ct_nac_now - stdE_ct_nac_now;
    [asc, ord]=areapiena(a,b,dim);
    hp_ctrl_sy=patch(asc,ord, 'green');
    set(hp_ctrl_sy, 'facecolor', [0.8 1 0.8], 'edgecolor', 'none');
    hold on
    p1=plot(mean_ct_nac_before,'black');
    hold on 
    p2=plot(mean_ct_nac_now,'green');
    min_all=min([min(mean_ct_nac_before - stdE_ct_nac_before) min(mean_pt_nac_before - stdE_pt_nac_before) min(mean_ct_nac_now - stdE_ct_nac_now) min(mean_pt_nac_now - stdE_pt_nac_now)]);
    max_all=max([max(mean_ct_nac_before + stdE_ct_nac_before) max(mean_pt_nac_before + stdE_pt_nac_before) max(mean_ct_nac_now + stdE_ct_nac_now) max(mean_pt_nac_now + stdE_pt_nac_now)]);
    gray = [0.1 0.1 0.1];
    dark = [0.4 0.4 0.4];
    % gray significant
    [h_avg_ctrl,p_avg_ctrl] = ttest2(Ct_nac_before(:,:),Ct_nac_now(:,:),'Alpha',0.05);
    flag=0;
    x_coo_nac=[0 0 0 0];
    for i=1:size(mean_ct_nac_before,2)
        if h_avg_ctrl(:,i)==1 && flag==0
            flag=1;
            x_coo_nac=[i 0 0 i];
        elseif h_avg_ctrl(:,i)==0 && flag==1
            x_coo_nac(2:3)=[i i];
            if (x_coo_nac(2)-x_coo_nac(1))>=15
                patch(x_coo_nac,[min_all min_all max_all max_all],gray,'EdgeColor','none','FaceAlpha',.1)
            end
            flag=0;
            x_coo_sy=[0 0 0 0];
        end
    end
    legend([p1 p2],{sprintf('Nac ct before (%d)',size(Ct_nac_before,1)),sprintf('Nac ct now (%d)',size(Ct_nac_now,1))},'FontSize',18)
    axis([0 614 min_all max_all])
    xticks([0 101 204 306 409 511])
    xticklabels({'-100','0','100','200','300','400','500'})
    xlabel('time (ms)','FontSize',18)
    ylabel('intensity (uV)','FontSize',18)
    set(gca,'FontSize',18);
    title(sprintf('Nac Controls %s Electrode %d (%s)',triga,elec,chan_names{elec}),'FontSize',18)

    subplot(2,1,2)
    a= mean_pt_nac_before + stdE_pt_nac_before;
    b=mean_pt_nac_before - stdE_pt_nac_before;
    [asc, ord]=areapiena(a,b,dim);
    hp_pt_nac=patch(asc,ord, 'black');
    set(hp_pt_nac, 'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none');
    hold on
    a= mean_pt_nac_now + stdE_pt_nac_now;
    b=mean_pt_nac_now - stdE_pt_nac_now;
    [asc, ord]=areapiena(a,b,dim);
    hp_pt_sy=patch(asc,ord, 'green');
    set(hp_pt_sy, 'facecolor', [0.8 1 0.8], 'edgecolor', 'none');
    hold on
    p1=plot(mean_pt_nac_before,'black');
    hold on
    p2=plot(mean_pt_nac_now,'green');
    hold on
    % gray significant
    [h_avg_sy,p_avg_sy] = ttest2(Pt_nac_before(:,:),Pt_nac_now(:,:),'Alpha',0.05);
    flag=0;
    x_coo_sy=[0 0 0 0];
    for i=1:size(mean_ct_nac_now,2)
        if h_avg_sy(:,i)==1 && flag==0
            flag=1;
            x_coo_sy=[i 0 0 i];
        elseif h_avg_sy(:,i)==0 && flag==1
            x_coo_sy(2:3)=[i i];
            if (x_coo_sy(2)-x_coo_sy(1))>=15
                patch(x_coo_sy,[min_all min_all max_all max_all],gray,'EdgeColor','none','FaceAlpha',.1)
            end
            flag=0;
            x_coo_sy=[0 0 0 0];
        end
    end
    legend([p1 p2],{sprintf('Nac pt before (%d)',size(Pt_nac_before,1)),sprintf('Nac pt now (%d)',size(Pt_nac_now,1))},'FontSize',18)
    title(sprintf('Nac Patients %s Electrode %d (%s)',triga,elec,chan_names{elec}),'FontSize',18)
    axis([0 614 min_all max_all])
    xticks([0 101 204 306 409 511])
    xticklabels({'-100','0','100','200','300','400','500'})
    xlabel('time (ms)','FontSize',18)
    ylabel('intensity (uV)','FontSize',18)
    set(gca,'FontSize',18);
    saveas(fig,filename)
    
    if CompSameSubj
        for ind=1:size(ind_ct_nac,2)
            su=Ct_nac_now_names{ind_ct_nac(1,ind)};
            fig=figure('Position',[10 10 1800 550]);
            filename=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Compare_old_new_prepro_nac_%s_%s_%s.png',triga,chan_names{elec},su);
            p1=plot(Ct_nac_now(ind_ct_nac(1,ind),:),'red');
            hold on
            p2=plot(Ct_nac_before(ind_ct_nac(2,ind),:),'blue');
            legend([p1 p2],{sprintf('%s (Nac ct) now',su),sprintf('%s (Nac ct) before',su)},'FontSize',18)
            title(sprintf('Nac Control %s now vs before %s Electrode %d (%s)',su,triga,elec,chan_names{elec}),'FontSize',18)
            %axis([0 614 min_all max_all])
            xticks([0 101 204 306 409 511])
            xticklabels({'-100','0','100','200','300','400','500'})
            xlabel('time (ms)','FontSize',18)
            ylabel('intensity (uV)','FontSize',18)
            set(gca,'FontSize',18);
            saveas(fig,filename)
            close(fig)
        end
        for ind=1:size(ind_pt_nac,2)
            su=Pt_nac_now_names{ind_pt_nac(1,ind)};
            fig=figure('Position',[10 10 1800 550]);
            filename=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Compare_old_new_prepro_nac_%s_%s_%s.png',triga,chan_names{elec},su);
            p1=plot(Pt_nac_now(ind_pt_nac(1,ind),:),'red');
            hold on
            p2=plot(Pt_nac_before(ind_pt_nac(2,ind),:),'blue');
            legend([p1 p2],{sprintf('%s (Nac pt) now',su),sprintf('%s (Nac pt) before',su)},'FontSize',18)
            title(sprintf('Nac Patient %s now vs before %s Electrode %d (%s)',su,triga,elec,chan_names{elec}),'FontSize',18)
            %axis([0 614 min_all max_all])
            xticks([0 101 204 306 409 511])
            xticklabels({'-100','0','100','200','300','400','500'})
            xlabel('time (ms)','FontSize',18)
            ylabel('intensity (uV)','FontSize',18)
            set(gca,'FontSize',18);
            saveas(fig,filename)
            close(fig)
        end
    end
    
    % Synapsy
    
    mean_ct_sy_before=mean(Ct_sy_before,1);
    mean_pt_sy_before=mean(Pt_sy_before,1);
    mean_ct_sy_now=mean(Ct_sy_now,1);
    mean_pt_sy_now=mean(Pt_sy_now,1);

    % standard error = std/sqrt(N-1)

    stdE_ct_sy_before=std(Ct_sy_before,1)/sqrt(size(Ct_sy_before,1)-1);
    stdE_pt_sy_before=std(Pt_sy_before,1)/sqrt(size(Pt_sy_before,1)-1);
    stdE_ct_sy_now=std(Ct_sy_now,1)/sqrt(size(Ct_sy_now,1)-1);
    stdE_pt_sy_now=std(Pt_sy_now,1)/sqrt(size(Pt_sy_now,1)-1);

    fig=figure('Position',[10 10 1800 550]);
    filename=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Compare_old_new_prepro_sy_%s_%s.png',triga,chan_names{elec});

    subplot(2,1,1)
    a=mean_ct_sy_before + stdE_ct_sy_before;
    b=mean_ct_sy_before - stdE_ct_sy_before;
    dim=614;
    [asc, ord]=areapiena(a,b,dim);
    hp_ctrl_sy=patch(asc,ord, 'black');
    set(hp_ctrl_sy, 'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none');
    hold on
    a= mean_ct_sy_now + stdE_ct_sy_now;
    b=mean_ct_sy_now - stdE_ct_sy_now;
    [asc, ord]=areapiena(a,b,dim);
    hp_ctrl_sy=patch(asc,ord, 'green');
    set(hp_ctrl_sy, 'facecolor', [0.8 1 0.8], 'edgecolor', 'none');
    hold on
    p1=plot(mean_ct_sy_before,'black');
    hold on 
    p2=plot(mean_ct_sy_now,'green');
    min_all=min([min(mean_ct_sy_before - stdE_ct_sy_before) min(mean_pt_sy_before - stdE_pt_sy_before) min(mean_ct_sy_now - stdE_ct_sy_now) min(mean_pt_sy_now - stdE_pt_sy_now)]);
    max_all=max([max(mean_ct_sy_before + stdE_ct_sy_before) max(mean_pt_sy_before + stdE_pt_sy_before) max(mean_ct_sy_now + stdE_ct_sy_now) max(mean_pt_sy_now + stdE_pt_sy_now)]);
    gray = [0.1 0.1 0.1];
    dark = [0.4 0.4 0.4];
    % gray significant
    [h_avg_ctrl,p_avg_ctrl] = ttest2(Ct_sy_before(:,:),Ct_sy_now(:,:),'Alpha',0.05);
    flag=0;
    x_coo_nac=[0 0 0 0];
    for i=1:size(mean_ct_sy_before,2)
        if h_avg_ctrl(:,i)==1 && flag==0
            flag=1;
            x_coo_nac=[i 0 0 i];
        elseif h_avg_ctrl(:,i)==0 && flag==1
            x_coo_nac(2:3)=[i i];
            if (x_coo_nac(2)-x_coo_nac(1))>=15
                patch(x_coo_nac,[min_all min_all max_all max_all],gray,'EdgeColor','none','FaceAlpha',.1)
            end
            flag=0;
            x_coo_sy=[0 0 0 0];
        end
    end
    legend([p1 p2],{sprintf('Synapsy ct before (%d)',size(Ct_sy_before,1)),sprintf('Synapsy ct now (%d)',size(Ct_sy_now,1))},'FontSize',18)
    axis([0 614 min_all max_all])
    xticks([0 101 204 306 409 511])
    xticklabels({'-100','0','100','200','300','400','500'})
    xlabel('time (ms)','FontSize',18)
    ylabel('intensity (uV)','FontSize',18)
    set(gca,'FontSize',18);
    title(sprintf('Synapsy Controls %s Electrode %d (%s)',triga,elec,chan_names{elec}),'FontSize',18)

    subplot(2,1,2)
    a= mean_pt_sy_before + stdE_pt_sy_before;
    b=mean_pt_sy_before - stdE_pt_sy_before;
    [asc, ord]=areapiena(a,b,dim);
    hp_pt_nac=patch(asc,ord, 'black');
    set(hp_pt_nac, 'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none');
    hold on
    a= mean_pt_sy_now + stdE_pt_sy_now;
    b=mean_pt_sy_now - stdE_pt_sy_now;
    [asc, ord]=areapiena(a,b,dim);
    hp_pt_sy=patch(asc,ord, 'green');
    set(hp_pt_sy, 'facecolor', [0.8 1 0.8], 'edgecolor', 'none');
    hold on
    p1=plot(mean_pt_sy_before,'black');
    hold on
    p2=plot(mean_pt_sy_now,'green');
    hold on
    % gray significant
    [h_avg_sy,p_avg_sy] = ttest2(Pt_sy_before(:,:),Pt_sy_now(:,:),'Alpha',0.05);
    flag=0;
    x_coo_sy=[0 0 0 0];
    for i=1:size(mean_ct_sy_now,2)
        if h_avg_sy(:,i)==1 && flag==0
            flag=1;
            x_coo_sy=[i 0 0 i];
        elseif h_avg_sy(:,i)==0 && flag==1
            x_coo_sy(2:3)=[i i];
            if (x_coo_sy(2)-x_coo_sy(1))>=15
                patch(x_coo_sy,[min_all min_all max_all max_all],gray,'EdgeColor','none','FaceAlpha',.1)
            end
            flag=0;
            x_coo_sy=[0 0 0 0];
        end
    end
    legend([p1 p2],{sprintf('Synapsy pt before (%d)',size(Pt_sy_before,1)),sprintf('Synapsy pt now (%d)',size(Pt_sy_now,1))},'FontSize',18)
    title(sprintf('Synapsy Patients %s Electrode %d (%s)',triga,elec,chan_names{elec}),'FontSize',18)
    axis([0 614 min_all max_all])
    xticks([0 101 204 306 409 511])
    xticklabels({'-100','0','100','200','300','400','500'})
    xlabel('time (ms)','FontSize',18)
    ylabel('intensity (uV)','FontSize',18)
    set(gca,'FontSize',18);
    saveas(fig,filename)

    if CompSameSubj
        for ind=1:size(ind_ct_sy,2)
            su=Ct_sy_now_names{ind_ct_sy(1,ind)};
            fig=figure('Position',[10 10 1800 550]);
            filename=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Compare_old_new_prepro_sy_%s_%s_%s.png',triga,chan_names{elec},su);
            p1=plot(Ct_sy_now(ind_ct_sy(1,ind),:),'red');
            hold on
            p2=plot(Ct_sy_before(ind_ct_sy(2,ind),:),'blue');
            legend([p1 p2],{sprintf('%s (Sy ct) now',su),sprintf('%s (Sy ct) before',su)},'FontSize',18)
            title(sprintf('Synapsy Control %s now vs before %s Electrode %d (%s)',su,triga,elec,chan_names{elec}),'FontSize',18)
            %axis([0 614 min_all max_all])
            xticks([0 101 204 306 409 511])
            xticklabels({'-100','0','100','200','300','400','500'})
            xlabel('time (ms)','FontSize',18)
            ylabel('intensity (uV)','FontSize',18)
            set(gca,'FontSize',18);
            saveas(fig,filename)
            close(fig)
        end
        for ind=1:size(ind_pt_sy,2)
            su=Pt_sy_now_names{ind_pt_sy(1,ind)};
            fig=figure('Position',[10 10 1800 550]);
            filename=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Compare_old_new_prepro_sy_%s_%s_%s.png',triga,chan_names{elec},su);
            p1=plot(Pt_sy_now(ind_pt_sy(1,ind),:),'red');
            hold on
            p2=plot(Pt_sy_before(ind_pt_sy(2,ind),:),'blue');
            legend([p1 p2],{sprintf('%s (Sy pt) now',su),sprintf('%s (Sy pt) before',su)},'FontSize',18)
            title(sprintf('Synapsy Patient %s now vs before %s Electrode %d (%s)',su,triga,elec,chan_names{elec}),'FontSize',18)
            %axis([0 614 min_all max_all])
            xticks([0 101 204 306 409 511])
            xticklabels({'-100','0','100','200','300','400','500'})
            xlabel('time (ms)','FontSize',18)
            ylabel('intensity (uV)','FontSize',18)
            set(gca,'FontSize',18);
            saveas(fig,filename)
            close(fig)
        end
    end
    clearvars -except CompSameSubj trigsa trigsb baseline_tp EEG chan_names trig_i
end
