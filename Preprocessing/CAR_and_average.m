clc
clear all
close all

db=2; % 1 or 2

OutDir=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Dataset%d/',db);

Bad_elec_file=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Elec2Interpolate_Database%d.xls',db);
[~,~,raw] = xlsread(Bad_elec_file);
Subj_names=raw(:,1);

All_Channel_num=1:64;
Outter_channels=[1 2 7 8 15 16 23:25 27:29 33:35 42:43 52:53 60:62 64];
Channel_num=All_Channel_num(setdiff(1:end,Outter_channels));

trig_type_name{1}='std';
trig_type_name{2}='dev1';
trig_type_name{3}='dev2';
trig_type_name{4}='dev3';
trig_type_name{5}='dev';
trig_type_name{6}='IC';
trig_type_name{7}='NIC';
trig_type_name{8}='IC1';
trig_type_name{9}='NIC1';
trig_type_name{10}='IC2';
trig_type_name{11}='NIC2';

baseline_tp=102;
  
for trig=1:size(trig_type_name,2)
    
    OutDir01=sprintf('%s/%s/ERPs',OutDir,trig_type_name{trig});
    if ~exist(OutDir01) 
        cd(sprintf('%s/%s/',OutDir,trig_type_name{trig}))
        mkdir('ERPs')
        cd('/Users/mip/Documents/Early-psychosis_Project/Preprocessing')
    end
    OutDir02=sprintf('%s/%s/GFPs',OutDir,trig_type_name{trig});
    if ~exist(OutDir02) 
        cd(sprintf('%s/%s/',OutDir,trig_type_name{trig}))
        mkdir('GFPs')
        cd('/Users/mip/Documents/Early-psychosis_Project/Preprocessing')
    end
    OutDir1=sprintf('%s/%s/epochs',OutDir,trig_type_name{trig});
    
    switch trig_type_name{trig}
        case 'std'
            all_trig={'std'};
        case 'dev'
            all_trig={'dev1','dev2','dev3'};
        case 'dev1'
            all_trig={'dev1'};
        case 'dev2'
            all_trig={'dev2'};
        case 'dev3'
            all_trig={'dev3'};
        case 'IC'
            all_trig={'IC1','IC2'};
        case 'NIC'
            all_trig={'NIC1','NIC2'};
        case 'IC1'
            all_trig={'IC1'};
        case 'NIC1'
            all_trig={'NIC1'};
        case 'IC2'
            all_trig={'IC2'};
        case 'NIC2'
            all_trig={'NIC2'};
    end
    
    for subj=1:size(Subj_names,1)
    
        Subj=Subj_names{subj};
        Npoc=1;
        for typ=1:size(all_trig,2)
            OutDir2=sprintf('%s/%s/epochs/%s',OutDir,all_trig{typ},Subj);

            try
                load(sprintf('%s/nb_epochs.mat',OutDir2));
                count_epoch(typ)=count_e-1;
                count_rej(typ)=count_r-1;
            catch
                count_epoch(typ)=0;
                count_rej(typ)=0;
            end

            if count_epoch(typ)
                files=dir(sprintf('%s/%s*.mat',OutDir2,Subj));
                for poc=Npoc:size(files,1)
                    load(sprintf('%s/%s',OutDir2,files(poc-Npoc+1).name));
                    % CAR (only inner channels)
                    for tp=1:size(epoch,2)
                        car_epoch(:,tp)=double(epoch(:,tp)-mean(epoch(Channel_num,tp)));
                    end
                    % Baseline correction
                    All_pocs(:,:,poc)=car_epoch-mean(car_epoch(:,1:baseline_tp),2);
                end
                Npoc=poc+1;
            end
            a=1;
        end
        count_ratio=sum(count_epoch)/(sum(count_epoch)+sum(count_rej));
        
        %if exist('All_pocs')
        if count_ratio>0.1
            ERP=mean(All_pocs,3);
            save(sprintf('%s/%s_ERP.mat',OutDir01,Subj),'ERP');
            
            GFP=std(ERP);
            save(sprintf('%s/%s_GFP.mat',OutDir02,Subj),'GFP');
        else
            fprintf('%s lack trials for %s \n',Subj,trig_type_name{trig})
        end
        clear All_pocs count_epoch count_rej
    end
       a=1;     
end

