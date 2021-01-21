clc
clear all
close all

db=1;

baseline_tp=102;
OutDir=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Dataset%d/',db);
load('/Users/mip/Documents/PdM/Data/Cap/chanlocs.mat')
for chan=1:size(EEG.chanlocs,2)
    chan_names{chan}=EEG.chanlocs(chan).labels;
end   
Bad_elec_file=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Elec2Interpolate_Database%d.xls',db);
[~,~,raw] = xlsread(Bad_elec_file);
Subj_names=raw(:,1);
Bad_elec_detect(2:size(Subj_names,1)+1,:)=Subj_names;

All_Channel_num=1:64;
Outter_channels=[1 2 7 8 15 16 23:25 27:29 33:35 42:43 52:53 60:62 64];
Channel_num=All_Channel_num(setdiff(1:end,Outter_channels));

Bad_elec_detect(:,2:size(Channel_num,2)+1)=num2cell(zeros([size(Bad_elec_detect,1) size(Channel_num,2)]));
Bad_elec_detect(1,2:end)=chan_names(Channel_num);

tri=1;

for trig={'std','dev','IC','NIC'}
    fprintf('%s \n',trig{1})
    
    OutDir01=sprintf('%s%s/ERPs',OutDir,trig{1});
    OutDir02=sprintf('%s%s/GFPs',OutDir,trig{1});
    
    Pt_ERPs=[];
    Ct_ERPs=[];

    Npt=0;
    Nct=0;
    
    for subj=1:size(Bad_elec_detect,1)

        Subj=Bad_elec_detect{subj};

        try
            load(sprintf('%s/%s_ERP.mat',OutDir01,Subj));
            error_loading=0;
        catch
            error_loading=1;
        end

        if ~error_loading
            if Subj(1)=='F' || strcmp(Subj(1:2),'Ln')
                Npt=Npt+1;
                ind=0;
                for elec=Channel_num
                    ind=ind+1;
                    Pt_ERPs{ind}(:,Npt)=ERP(elec,:);
                end
                Pt_names{:,Npt}=Subj;
            else
                Nct=Nct+1;
                ind=0;
                for elec=Channel_num
                    ind=ind+1;
                    Ct_ERPs{ind}(:,Nct)=ERP(elec,:);
                end
                Ct_names{:,Nct}=Subj;
            end
        end
        a=1;
        clear count_epoch count_rej
    end
    
    for elec=1:size(Channel_num,2)
        Avg_Pt=mean(Pt_ERPs{elec},2);
        Avg_Ct=mean(Ct_ERPs{elec},2);
        All_ERPs=[Pt_ERPs{elec} Ct_ERPs{elec}];
        All_names=[Pt_names Ct_names];
        Avg_All=mean(All_ERPs,2);
        
        for su=1:size(All_ERPs,2)
    
            Corrval=corr(All_ERPs(:,su),Avg_All);
            [index,b]=find(strcmp(All_names{su},Bad_elec_detect)==1);
            
            if Corrval>0.4
                Bad_elec_detect{index,elec+1}=Bad_elec_detect{index,elec+1}+1;
                if Bad_elec_detect{index,elec+1}>4
                    disp('what?')
                end
            end
        end
    
    end
    
    tri=tri+1;
    
    clear Ct_ERPs Pt_ERPs Avg_All Avg_Ct Avg_Pt All_names All_ERPs Pt_names Ct_names
    
end

figure
imagesc(cell2mat(Bad_elec_detect(2:end,2:end)'))
xticks(1:size(Bad_elec_detect,1)-1)
xticklabels(Bad_elec_detect(2:end,1))
xtickangle(90)
yticks(1:size(Bad_elec_detect,2)-1)
yticklabels(Bad_elec_detect(1,2:end))
colorbar
%set(gca,'FontSize',8)
a=1;

for su=1:size(Subj_names,1)
    indexes=find(cell2mat(Bad_elec_detect(su+1,2:end)')<=2==1)+1;
    str='';
    for it=1:size(indexes,1)
        str=sprintf('%s %s',str,Bad_elec_detect{1,indexes(it,1)});
    end
    Subj_names{su,2}=str;
end

