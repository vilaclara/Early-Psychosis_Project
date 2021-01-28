clc
clear all
close all

% Path to HOSVD toolbox
addpath('/Users/mip/Documents/PdM/Code/Matlab/tensor_toolbox-master');
%T = hosvd(X,tol,varargin)


Types=[{'std'}, {'dev'}, {'IC'}, {'NIC'}];%, {'IC_shape1'}, {'IC_shape2'}, {'NIC_shape1'}, {'NIC_shape2'}

Recompute=1;

for cas=1:size(Types,2)
    
    % Database 
    % 1 for Nac ct and pt
    % 2 for Sy ct and pt
    % 3 for Nac ct and Sy pt 
    % 4 for Sy ct and Nac pt
    % 5 for all ct and all pt

    DB=2;

    switch DB
        case 1
            ct_nac=1;
            ct_sy=0;
            pt_nac=1;
            pt_sy=0;
        case 2
            ct_nac=0;
            ct_sy=1;
            pt_nac=0;
            pt_sy=1;
    end

    data_names=[{'Standard'},{'Deviant'},{'Illusory Contour 1'},{'Illusory Contour 2'},{'Non Illusory Contour 1'},{'Non Illusory Contour 2'},{'Illusory Contour'},{'Non Illusory Contour'}];
    Types=[{'std'}, {'dev'}, {'IC'}, {'NIC'}];%, {'IC_shape1'}, {'IC_shape2'}, {'NIC_shape1'}, {'NIC_shape2'}

    CurrentPath = addpath('/Users/mip/Documents/PdM/Code/Matlab/');
    
    data_type=Types{cas};
    switch data_type
        case 'std'
            data_name=data_names{1};
        case 'dev'
            data_name=data_names{2};
        case 'IC_shape1'
            data_name=data_names{3};
        case 'IC_shape2'
            data_name=data_names{4};
        case 'NIC_shape1'
            data_name=data_names{5};
        case 'NIC_shape2'
            data_name=data_names{6};
        case 'IC'
            data_name=data_names{7};
        case 'NIC'
            data_name=data_names{8};
    end
    
    OutputDir = sprintf('/Users/mip/Documents/PdM/Results/HOSVD_results');
    
    switch DB
        case 1
            theDataIs='Nac';
        case 2
            theDataIs='Sy';
        case 5
            theDataIs='Both';
    end
    
    n_comp = 2;
    
    %% Take off outter channels
    
    All_Channel_num=1:64;
    All_Channel_names={'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1',...
    'C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9',...
    'PO7','PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4',...
    'Afz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2',...
    'C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8',...
    'PO4','O2'};

    %Outter_channels=[1 2 7 8 15 16 23:25 27:29 33:35 42:43 52:53 60:62 64];
    %Channel_num=All_Channel_num(setdiff(1:end,Outter_channels));
    %Channel_names=All_Channel_names(Channel_num);
    
    filename=sprintf('/Users/mip/Documents/PdM/Code/Matlab/Inner_channels.mat');
    load(filename); %,'Channel_num','Channel_names'
    
    n_channel = size(Channel_num,2);
    
    EpochPath_nac = '/Users/mip/Documents/PdM/Data/ERPs/Dataset1';
    EpochPath_sy = '/Users/mip/Documents/PdM/Data/ERPs/Dataset2';
    
    MainDir_sy = sprintf('%s/%s/epochs',EpochPath_sy,data_type);
    AvgDir_sy = sprintf('%s/%s/ERPs',EpochPath_sy,data_type);
    MainDir_nac = sprintf('%s/%s/epochs',EpochPath_nac,data_type);
    AvgDir_nac = sprintf('%s/%s/ERPs',EpochPath_nac,data_type);

    %% Load Data

    Subj_names={};
    Groups=[];
    X=[];

    Avg_Controls=[];

    nb_ctl_avg=1;
    Ctl_code_dist=[];

    if ct_nac

        span_ct = 1:60; 
        for Subj_num = span_ct
            Subj_name = sprintf('Ctrl%03d',Subj_num);
            nom_avg=sprintf('%s_ERP.mat',Subj_name);
            if exist(sprintf('%s/%s',AvgDir_nac,nom_avg), 'file')
                load(sprintf('%s/%s',AvgDir_nac,nom_avg));
                Ctl_code_dist=[Ctl_code_dist Subj_num];
                Subj_names=[Subj_names {Subj_name}];
                Groups=[Groups 1];
                i=1;
                for elec = Channel_num %some have 65 electrodes???????
                    Avg_Controls{nb_ctl_avg}(i,:) = ERP(elec,:)';
                    i=i+1;
                end
                nb_ctl_avg = nb_ctl_avg+1;
            end
        end

    end

    if ct_sy

        span_ct=350:570;
        for Subj_num = span_ct
            Subj_name = sprintf('L%03d',Subj_num);
            nom_avg=sprintf('%s_ERP.mat',Subj_name);
            if exist(sprintf('%s/%s',AvgDir_sy,nom_avg), 'file')
                load(sprintf('%s/%s',AvgDir_sy,nom_avg));
                Ctl_code_dist=[Ctl_code_dist Subj_num];
                Subj_names=[Subj_names {Subj_name}];
                Groups=[Groups 1];
                i=1;
                for elec = Channel_num %some have 65 electrodes???????
                    Avg_Controls{nb_ctl_avg}(i,:) = ERP(elec,:)';
                    i=i+1;
                end
                nb_ctl_avg = nb_ctl_avg+1;
            end

        end
    end


    Avg_Patients=[];
    GFP_Patients=[];
    nb_pt_avg=1;
    nb_pt_gfp=1;
    Pt_code_dist=[];

    if pt_nac

        span_pt = 1:60;
        for Subj_num = span_pt
            Subj_name = sprintf('Lnac%03d',Subj_num);
            nom_avg=sprintf('%s_ERP.mat',Subj_name);
            if exist(sprintf('%s/%s',AvgDir_nac,nom_avg), 'file')
                Pt_code_dist=[Pt_code_dist Subj_num];
                load(sprintf('%s/%s',AvgDir_nac,nom_avg));
                Subj_names=[Subj_names, {Subj_name}];
                Groups=[Groups 2];
                i=1;
                for elec = Channel_num %some have 65 electrodes???????
                    Avg_Patients{nb_pt_avg}(i,:) = ERP(elec,:)';
                    i=i+1;
                end
                nb_pt_avg = nb_pt_avg+1;
            end
        end

    end

    if pt_sy

        span_pt=100:230;
        for Subj_num = span_pt
            Subj_name = sprintf('F%03d',Subj_num);
            nom_avg=sprintf('%s_ERP.mat',Subj_name);
            %nom_gfp=sprintf('%s_GFP.mat',Subj_name);

            if exist(sprintf('%s/%s',AvgDir_sy,nom_avg), 'file')
                Pt_code_dist=[Pt_code_dist Subj_num];
                load(sprintf('%s/%s',AvgDir_sy,nom_avg));
                Subj_names=[Subj_names, {Subj_name}];
                Groups=[Groups 2];
                i=1;
                for elec = Channel_num %some have 65 electrodes???????
                    Avg_Patients{nb_pt_avg}(i,:) = ERP(elec,:)';
                    i=i+1;
                end
                nb_pt_avg = nb_pt_avg+1;
            end
        end

    end



    %% Concatenate
    
    clear X_mat X T U Y U1 U2 U3 F_U1 F_U2 F_U3

    X=tenrand(size(Avg_Controls,2)+size(Avg_Patients,2),size(Avg_Controls{1},1),size(Avg_Controls{1},2));
    X_mat=zeros(size(Avg_Controls,2)+size(Avg_Patients,2),size(Avg_Controls{1},1),size(Avg_Controls{1},2));
    % [subjects x channels x time]

    for subj=1:size(Avg_Controls,2) % nb controls
        X(subj,:,:)=Avg_Controls{subj};
        X_mat(subj,:,:)=Avg_Controls{subj};
    end

    for subj=1:size(Avg_Patients,2) % nb patient
        X(subj+size(Avg_Controls,2),:,:)=Avg_Patients{subj};
        X_mat(subj+size(Avg_Controls,2),:,:,1)=Avg_Patients{subj};
    end
    
    % zscore
%    X_mat=zscore(X_mat,0,1);
%    X=tensor(X_mat);

%     grand_avg=squeeze(mean(X_mat,2));
%     
%     for s=1:size(X,1)
%         figure
%         plot(grand_avg(s,:))
%         title(sprintf('%s',Subj_names{s}))
%         a=1;
%         close all
%         
%     end

    % Try putting everything positive !
    %OutputDir = sprintf('/Users/mip/Desktop/Results for Elvira feb 2020/Trying_to_solve_the_HOSVD_comp_direction/Shift_so_everything_pos');
    %the_min=min(min(min(X_mat)));
    %X_mat=X_mat-the_min;
    %X=X-the_min;

    clear Avg_Patients Avg_Patients

    n_subj_all = size(X,1);
    n_chan_all = size(X,2);
    n_time_all = size(X,3);
    
% ___________________    % test indiv svd direction   ___________

%     if (strcmp(data_type,'std')||strcmp(data_type,'dev')) 
%         elec=20;
%     elseif (strcmp(data_type,'visu_illu')||strcmp(data_type,'visu_nn_illu'))
%         elec=18;%elec=15;
%     end
%     for n=1:n_subj_all
%         [U,S,V]=svd(squeeze(X_mat(n,:,:)));
%         figure
%         plot(V(:,1))
%         
%         a=1;
%         close all
%         
%     end

    %% Singular value decomposition of XTX

    tol=10^(-7);
    [T, U, Y] = hosvd(X,tol);

    %[U,S,V]=svd(XTX,'econ');
    % U is topographical map = PCA
    % V temporal course of topography
    % S the eigen value (the importance of each topography)

%         if mean(U{1}(:,1))<0
%             for i=1:3
%                 U{i}(:,1)=-1*U{i}(:,1);
%                 U{i}(:,2)=-1*U{i}(:,2);
%             end
%         end

    U1=U{1};
    U2=U{2};
    % U2 is topographical map 
    U3=U{3};
    % U3 temporal course of topography

    PC1_svd=U{3}(:,1);

%         len = size(U{3});
%         save(sprintf('%s/sizeAll_%s.mat',OutputDir,data_type), 'len');
%         save(sprintf('%s/S_for_screeplot_%s.mat',OutputDir,data_type), 'S');
%         S_diag = diag(S);
%         ScreePlot(S_diag(1:10),n_channel);
%         title(sprintf('Screeplot %s',data_type))
%         save([OutputDir '/ScreePlot_SVD_' data_type '.fig'])
%         
%         OutputDir = sprintf('/Volumes/MyPassportforMac/PdM/Results_sy_ct_sy_pt/PCA_SVD_results');
%         data_type='std'
%         load(sprintf('%s/S_for_screeplot_%s.mat',OutputDir,data_type), 'S');
%         SS = diag(S);
%         EV=cumsum(SS.^2/sum(SS.^2))*100;
%         EV(1)
%         EV(2)-EV(1)

    for nc1 = 1:n_comp
        for nc2 = 1:n_comp
            corr_all(nc1,nc2) = corr(U3(:,nc1),U3(:,nc2));
        end
    end
    save(sprintf('%s/CorrAll_%d_%s.mat',OutputDir,DB,data_type), 'corr_all');

    clear corr_all

    F_U3 = U3(:,1:n_comp);
    F_U2 = U2(:,1:n_comp);
    F_U1 = U1(:,1:n_comp);
    % Take the absolute value of the subjects
    %F_U1 = abs(U1(:,1:n_comp));
    
    % If everything put positive, the 1st component contains about the mean
%     F_U3 = U3(:,2:n_comp+1);
%     F_U2 = U2(:,2:n_comp+1);
%     F_U1 = U1(:,2:n_comp+1);

    my_hosvd_plot_saliencies(F_U1,F_U2,F_U3,data_name,theDataIs,OutputDir,Subj_names)
    
     [S,~]=hosvd_compute_variance(X,F_U1,F_U2,F_U3,n_comp);
     
     % _________      Try to look at the reconstruction of the data from the PCs       __________
     
%     [S1,recon1]=hosvd_compute_variance(X,U1,U2,U3,1);
%     [S2,recon2]=hosvd_compute_variance(X,U1,U2,U3,2);
%     [S3,recon3]=hosvd_compute_variance(X,U1,U2,U3,3);
%     
%     recon1=double(recon1);
%     recon2=double(recon2);
%     recon3=double(recon3);
%     
%     % try reconstruction to look at direction
%     if (strcmp(data_type,'std')||strcmp(data_type,'dev')) 
%         elec=20;
%     elseif (strcmp(data_type,'visu_illu')||strcmp(data_type,'visu_nn_illu'))
%         elec=18;%elec=15;
%     end
%     for n=1:n_subj_all
%         % correlation PC1
%         fig=figure
%         plot(squeeze(X_mat(n,elec,:)),'LineWidth',2)
%         hold on
%         plot(squeeze(recon1(n,elec,:)),'LineWidth',2)
%         hold on
%         plot(squeeze(recon2(n,elec,:)),'LineWidth',2)
%         hold on
%         plot(squeeze(recon3(n,elec,:)),'LineWidth',2)
%         legend('X','Reconstruction 1PC','Reconstruction 2PC','Reconstruction 3PC')
%         xticks([0 102 153 204 255 306 357 408 511 612])
%         xticklabels({'-100','0','','100','','200','','300','350','400','500'})
%         xlabel('Time (ms)','FontSize',18);
%         title(sprintf('%s %s Usu=%.1f Utopo=%.1f',Channel_names{1,elec}, Subj_names{1,n},F_U1(n,1),F_U2(elec,1)))
%         filename=(sprintf('%s/HOSVD_reconstruction/%s_%s_%d.png',OutputDir,data_type,Subj_names{1,n},DB))
%         saveas(fig,filename);
%         close all
%     end


    save(sprintf('%s/ProjectedActivity_%d_%s.mat',OutputDir,DB,data_type), 'F_U1', 'F_U2', 'F_U3','S','-v7.3');

%         len_sum = [0 cumsum(len_all)];
%         k_all = 1;

%         for sub = 1:nb_ctl_avg-1+nb_pt_avg-1
%             base_use = X(1+n_channel*(sub-1):n_channel*sub,:);
%             U_proc = ((F_V/diag(diag(S(:,1:n_comp))))'*base_use')';
%             V_proc = ((F_U3/diag(diag(S(:,1:n_comp))))'*base_use')';
%             k_all = k_all+1;
%             save([OutputDir '/' Subj_names{sub} '_U_and_V_proc_' data_type '.mat'], 'U_proc', 'V_proc')
%             clear base_use U_proc V_proc
%         end

    %% Effect of PC2 on PC1 all

%         suc=[{'all'},{'ct'},{'pt'}];
%         % all
%         std_coeff_PC2{1}=std(F_U1(:,2));
%         mean_coeff_PC1{1}=mean(F_U1(:,1));
%         mean_coeff_PC2{1}=mean(F_U1(:,2));
% 
%         %ct
%         std_coeff_PC2{2}=std(F_U1(Groups==1,2));
%         mean_coeff_PC1{2}=mean(F_U1(Groups==1,1));
%         mean_coeff_PC2{2}=mean(F_U1(Groups==1,2));
%         
%         %pt
%         std_coeff_PC2{3}=std(F_U1(Groups==2,2));
%         mean_coeff_PC1{3}=mean(F_U1(Groups==2,1));
%         mean_coeff_PC2{3}=mean(F_U1(Groups==2,2));
%         n_curves=12;
%         for i=1:3
%             span_alpha{i}=[mean_coeff_PC2{i}-std_coeff_PC2{i} mean_coeff_PC2{i}+std_coeff_PC2{i}];
%             alpha{i}=linspace(span_alpha{i}(1),span_alpha{i}(2),n_curves-1);
%             k=1;
%             for n=1:n_curves-2
%                 if alpha{i}(n)*alpha{i}(n+1)<0
%                     new_alpha(k)=alpha{i}(n);
%                     new_alpha(k+1)=0;
%                     k=k+1;
%                 else
%                     new_alpha(k)=alpha{i}(n);
%                 end
%                 k=k+1;
%             end
%             new_alpha(k)=alpha{i}(n+1);
%             alpha{i}=new_alpha;
%         end
%         
%         fig=figure;
%         
%         for k=1:3
%             subplot(3,1,k)
%             for i=1:n_curves
%                 % positive
%                 Comb_PC1_and2=mean_coeff_PC1{k}*F_U3(:,1)+alpha{k}(i)*F_U3(:,2);
%                 if alpha{k}(i)==0
%                     color='black';
%                 elseif alpha{k}(i)>0
%                     color=[0.4+i/(2*n_curves) 0.4 0.4];
%                 elseif alpha{k}(i)<0
%                     color=[0.4 0.4 0.4+i/(2*n_curves)];
%                 end
%                 plot(Comb_PC1_and2,'color',color)
%                 hold on
%             end
%             title(sprintf('Effect of PC2 on PC1 %s %s %s',theDataIs,data_name,suc{k}),'FontSize',18)
%             ylim([-0.03 0.03])
%         end
%         filename=sprintf('/Users/laura/Desktop/HOSVD_results/Effect of PC2 on PC1 %s %s.png',theDataIs,data_name);
%         saveas(fig,filename)
%         

    %% Norm

    PC1_ct_su_comp{cas}=U1(Groups==1,1);
    PC2_ct_su_comp{cas}=U1(Groups==1,2);

    PC1_pt_su_comp{cas}=U1(Groups==2,1);
    PC2_pt_su_comp{cas}=U1(Groups==2,2);

    %[h1_PC1(ta),p_PC1]=ttest2(PC1_pt_norm{ta},PC1_ct_norm{ta});
    [p_PC1(cas),h1_PC1(cas)]=ranksum(PC1_pt_su_comp{cas},PC1_ct_su_comp{cas});
    fprintf('Norm %s %s: %d \n',data_name,theDataIs,h1_PC1(cas))

    %[h1_PC2(ta),p_PC2]=ttest2(PC2_pt_norm{ta},PC2_ct_norm{ta});
    [p_PC2(cas),h1_PC2(cas)]=ranksum(PC2_pt_su_comp{cas},PC2_ct_su_comp{cas});
    fprintf('Norm %s %s: %d \n',data_name,theDataIs,h1_PC1(cas))


     save(sprintf('%s/Subject_names_%d_%s.mat',OutputDir,DB,data_type), 'Subj_names','Groups');

    clearvars -except p_PC1 h1_PC1 h1_PC2 p_PC2 PC1_ct_su_comp PC2_ct_su_comp  PC1_pt_su_comp PC2_pt_su_comp n_channel n_comp cas Recompute DB theDataIs data_names Types data_type data_name OutputDir

end

fig=figure;
model_series=[mean(PC1_pt_su_comp{1}) mean(PC1_ct_su_comp{1}) ; mean(PC1_pt_su_comp{2}) mean(PC1_ct_su_comp{2}) ; mean(PC1_pt_su_comp{3}) mean(PC1_ct_su_comp{3}) ; mean(PC1_pt_su_comp{4}) mean(PC1_ct_su_comp{4})];
model_error=[std(PC1_pt_su_comp{1})/sqrt(size(PC1_pt_su_comp{1},1)) std(PC1_ct_su_comp{1})/sqrt(size(PC1_ct_su_comp{1},1)) ; std(PC1_pt_su_comp{2})/sqrt(size(PC1_pt_su_comp{2},1)) std(PC1_ct_su_comp{2})/sqrt(size(PC1_ct_su_comp{2},1)) ; std(PC1_pt_su_comp{3})/sqrt(size(PC1_pt_su_comp{3},1)) std(PC1_ct_su_comp{3})/sqrt(size(PC1_ct_su_comp{3},1)) ; std(PC1_pt_su_comp{4})/sqrt(size(PC1_pt_su_comp{4},1)) std(PC1_ct_su_comp{4})/sqrt(size(PC1_ct_su_comp{4},1))];
bar(model_series)
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'LineWidth', 2, 'linestyle', 'none');
end
for ta=1:4
    if h1_PC1(ta)
        text(ta,mean(PC1_ct_su_comp{ta}),'*','FontSize',19)
    end
end
legend('pt','ctrl')
xticklabels({'Standard','Deviant','IC','NIC'})
set(gca,'FontSize',18)
%xlabel('Norm','FontSize',18);
%ylabel('Mean/Std','FontSize',18);
title(sprintf('Subj comp 1 no hungarian %s\n\n',theDataIs),'FontSize',18);
%filename=sprintf('/Users/laura/Desktop/Barplot_Norm/Norm_PC1_%s.png',DB_name);
filename=sprintf('%s/Subj_bar_PC1_%s.png',OutputDir,theDataIs);
saveas(fig,filename)

fig=figure;
model_series=[mean(PC2_pt_su_comp{1}) mean(PC2_ct_su_comp{1}) ; mean(PC2_pt_su_comp{2}) mean(PC2_ct_su_comp{2}) ; mean(PC2_pt_su_comp{3}) mean(PC2_ct_su_comp{3}) ; mean(PC2_pt_su_comp{4}) mean(PC2_ct_su_comp{4})];
model_error=[std(PC2_pt_su_comp{1})/sqrt(size(PC2_pt_su_comp{1},1)) std(PC2_ct_su_comp{1})/sqrt(size(PC2_ct_su_comp{1},1)) ; std(PC2_pt_su_comp{2})/sqrt(size(PC2_pt_su_comp{2},1)) std(PC2_ct_su_comp{2})/sqrt(size(PC2_ct_su_comp{2},1)) ; std(PC2_pt_su_comp{3})/sqrt(size(PC2_pt_su_comp{3},1)) std(PC2_ct_su_comp{3})/sqrt(size(PC2_ct_su_comp{3},1)) ; std(PC2_pt_su_comp{4})/sqrt(size(PC2_pt_su_comp{4},1)) std(PC2_ct_su_comp{4})/sqrt(size(PC2_ct_su_comp{4},1))];
bar(model_series)
hold on
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'LineWidth', 2, 'linestyle', 'none');
end
for ta=1:4
    if h1_PC2(ta)
        text(ta,mean(PC2_ct_su_comp{ta}),'*','FontSize',19)
    end
end
legend('pt','ctrl')
xticklabels({'Standard','Deviant','IC','NIC'})
set(gca,'FontSize',18)
title(sprintf('Subj comp 2 no hungarian %s\n\n',theDataIs),'FontSize',18);
%filename=sprintf('/Users/laura/Desktop/Barplot_Norm/Norm_PC2_%s.png',DB_name);
%filename=sprintf('/Users/laura/Desktop/HOSVD_results/Subj_bar_PC2_%s.png',theDataIs);
filename=sprintf('%s/Subj_bar_PC2_%s.png',OutputDir,theDataIs);
saveas(fig,filename)

% Plot

fig=figure;
subplot(2,2,1)
scatter(PC1_pt_su_comp{1},PC2_pt_su_comp{1},'or')
hold on
scatter(PC1_ct_su_comp{1},PC2_ct_su_comp{1},'ob')
legend('pt','ct')
title('Standard')
set(gca,'FontSize',14)
subplot(2,2,2)
scatter(PC1_pt_su_comp{2},PC2_pt_su_comp{2},'or')
hold on
scatter(PC1_ct_su_comp{2},PC2_ct_su_comp{2},'ob')
legend('pt','ct')
title('Deviant')
set(gca,'FontSize',14)
subplot(2,2,3)
scatter(PC1_pt_su_comp{3},PC2_pt_su_comp{3},'or')
hold on
scatter(PC1_ct_su_comp{3},PC2_ct_su_comp{3},'ob')
legend('pt','ct')
title('IC')
set(gca,'FontSize',14)
subplot(2,2,4)
scatter(PC1_pt_su_comp{4},PC2_pt_su_comp{4},'or')
hold on
scatter(PC1_ct_su_comp{4},PC2_ct_su_comp{4},'ob')
legend('pt','ct')
title('NIC')
set(gca,'FontSize',14)
%xlabel('Norm','FontSize',18);
%ylabel('Mean/Std','FontSize',18);
%title(sprintf('Plot of distribution %s\n\n',theDataIs),'FontSize',18);
%filename=sprintf('/Users/laura/Desktop/Barplot_Norm/Norm_PC1_%s.png',DB_name);
filename=sprintf('%s/Subj_dist_PC1_PC2_%s.png',OutputDir,theDataIs);
saveas(fig,filename)