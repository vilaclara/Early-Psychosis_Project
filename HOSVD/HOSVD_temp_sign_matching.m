clc
clear all
close all

Types=[{'std'}, {'dev'}, {'visu_illu'}, {'visu_nn_illu'}];%, {'IC_shape1'}, {'IC_shape2'}, {'NIC_shape1'}, {'NIC_shape2'}
Types_names=[{'std'}, {'dev'}, {'IC'}, {'NIC'}];%, {'IC_shape1'}, {'IC_shape2'}, {'NIC_shape1'}, {'NIC_shape2'}

% OutputDir_nac = sprintf('/Volumes/MyPassportforMac/PdM/Results_nac_ct_nac_pt/PCA_HOSVD_results');
% OutputDir_sy = sprintf('/Volumes/MyPassportforMac/PdM/Results_sy_ct_sy_pt/PCA_HOSVD_results');

OutputDir = sprintf('/Users/mip/Documents/PdM/Results/HOSVD_results');
    
fig=figure('units','normalized','outerposition',[0 0 1 1]);
nb_PC=2;

for cas=1:size(Types,2)

    data_type=Types{cas};

    load(sprintf('%s/ProjectedActivity_0_%s.mat',OutputDir,data_type)); % , 'F_U1', 'F_U2', 'F_U3'
    U_subj_nac{cas}=F_U1;
    U_space_nac{cas}=F_U2;
    U_time_nac{cas}=F_U3;
    S_nac{cas}=S;
    load(sprintf('%s/Subject_names_0_%s.mat',OutputDir,data_type));%, 'Subj_names','Groups'
    Subj_names_nac{cas}=Subj_names;
    
    load(sprintf('%s/ProjectedActivity_1_%s.mat',OutputDir,data_type)); % , 'F_U1', 'F_U2', 'F_U3'
    U_subj_sy{cas}=F_U1;
    U_space_sy{cas}=F_U2;
    U_time_sy{cas}=F_U3;
    S_sy{cas}=S;
    load(sprintf('%s/Subject_names_1_%s.mat',OutputDir,data_type));%, 'Subj_names','Groups'
    Subj_names_sy{cas}=Subj_names;
    
    %PC1-2-3
    for i=1:nb_PC
        subplot(4,nb_PC,nb_PC*cas-nb_PC+i)
        plot(U_time_nac{cas}(:,i),'black','LineWidth',1.2);
        hold on
        plot(U_time_sy{cas}(:,i),'blue','LineWidth',1.2);
        legend(sprintf('nac ev=%.2f',S_nac{cas}(i,i)),sprintf('sy ev=%.2f',S_sy{cas}(i,i)))
        xticks([0 100 150 200 250 300 350 400 500 600])
        xticklabels({'-100','0','','100','','200','','300','400','500'})
        title(sprintf('%s PC%d',Types_names{cas},i))
    end
    
    a=1;
end

clear F_U1 F_U2 F_U3 OutputDir_nac OutputDir_sy

OutputDir = '/Users/mip/Documents/PdM/Code/Matlab/HOSVD_sign_matched_results';
%OutputDir_sy = OutputDir;
%OutputDir_nac=OutputDir;

fileName=sprintf('%s/Before_matching.png',OutputDir);
saveas(fig,fileName)

%% Which to inverse (PCs 1 2 and 3 NOT coupled)

% Inv_mat_nac=[+1 +1 +1 ; ... % std
%              +1 +1 +1 ; ... % dev
%              +1 +1 +1 ; ... % IC
%              +1 +1 +1 ];    % NIC
%          
% Inv_mat_syn=[+1 +1 +1 ; ... % std
%              +1 +1 +1 ; ... % dev
%              +1 +1 +1 ; ... % IC
%              +1 +1 +1 ];    % NIC

%(PCs 1 2 and 3 coupled)

Inv_mat_nac=[-1 -1 -1 ; ... % std
             -1 -1 -1 ; ... % dev
             +1 +1 +1 ; ... % IC
             +1 +1 +1 ];    % NIC
             
         
Inv_mat_syn=[-1 -1 -1 ; ... % std
             -1 -1 -1 ; ... % dev
             -1 -1 -1 ; ... % IC
             -1 -1 -1 ];    % NIC
         
%% inverse

fig=figure('units','normalized','outerposition',[0 0 1 1]);

for cas=1:size(Types,2)

    data_type=Types{cas};
    
    %PC1-2-3
    for i=1:nb_PC
        
        %nac_U1{cas}(:,i)=U_subj_nac{cas}(:,i)*Inv_mat_nac(cas,i); 
        nac_U1{cas}(:,i)=U_subj_nac{cas}(:,i); 
        % change if not taking the absolute value!!!
        nac_U2{cas}(:,i)=U_space_nac{cas}(:,i); 
        nac_U3{cas}(:,i)=U_time_nac{cas}(:,i)*Inv_mat_nac(cas,i);
        
        %sy_U1{cas}(:,i)=U_subj_sy{cas}(:,i)*Inv_mat_syn(cas,i);
        sy_U1{cas}(:,i)=U_subj_sy{cas}(:,i)
        % change if not taking the absolute value!!!
        sy_U2{cas}(:,i)=U_space_sy{cas}(:,i); 
        sy_U3{cas}(:,i)=U_time_sy{cas}(:,i)*Inv_mat_syn(cas,i);
        
        [rho_nac_sy(cas,i) p_nac_sy(cas,i)] = corr(nac_U3{cas}(:,i),sy_U3{cas}(:,i),'type','Spearman');
    
        subplot(4,nb_PC,nb_PC*cas-nb_PC+i)
        plot(nac_U3{cas}(:,i),'black','LineWidth',1.2);
        hold on
        plot(sy_U3{cas}(:,i),'blue','LineWidth',1.2);
        legend(sprintf('nac ev=%.2f',S_nac{cas}(i,i)),sprintf('sy ev=%.2f',S_sy{cas}(i,i)))
        xticks([0 102 153 204 255 306 408 511 612])
        xticklabels({'-100','0','','100','','200','','300','350','400','500'})
        title(sprintf('%s PC%d (sign matched) rho=%f',Types_names{cas},i,rho_nac_sy(cas,i)))
    
    end
    
    F_U1=nac_U1{cas};
    F_U2=nac_U2{cas};
    F_U3=nac_U3{cas};
    Subj_names=Subj_names_nac{cas};
    S=S_nac;
    save(sprintf('%s/ProjectedActivity_sign_matched_%s_nac.mat',OutputDir,data_type), 'F_U1', 'F_U2', 'F_U3','Subj_names','S','-v7.3');

    F_U1=sy_U1{cas};
    F_U2=sy_U2{cas};
    F_U3=sy_U3{cas};
    Subj_names=Subj_names_sy{cas};
    S=S_sy;
    save(sprintf('%s/ProjectedActivity_sign_matched_%s_sy.mat',OutputDir,data_type), 'F_U1', 'F_U2', 'F_U3','Subj_names','S','-v7.3');

    %clear nac_U1 nac_U2 nac_U3 sy_U1 sy_U2 sy_U3
    
end

% for i=1:3
%     [rho_std_dev_nac(:,i) p_std_dev_nac(:,i)] = corr(nac_U3{1}(:,i),nac_U3{2}(:,i),'type','Spearman');
%     [rho_IC_NIC_nac(:,i) p_std_dev_nac(:,i)] = corr(nac_U3{3}(:,i),nac_U3{4}(:,i),'type','Spearman');
%     [rho_std_dev_sy(:,i) p_std_dev_sy(:,i)] = corr(sy_U3{1}(:,i),sy_U3{2}(:,i),'type','Spearman');
%     [rho_IC_NIC_sy(:,i) p_std_dev_sy(:,i)] = corr(sy_U3{3}(:,i),sy_U3{4}(:,i),'type','Spearman');
% end
% 
% rho_std_dev_nac
% rho_std_dev_sy
% rho_IC_NIC_nac
% rho_IC_NIC_sy

fileName=sprintf('%s/After_matching.png',OutputDir);
saveas(fig,fileName)


