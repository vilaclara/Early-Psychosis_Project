clc
clear all
close all

OutDir='/Users/mip/Documents/PdM/Data/ERPs/Dataset1';
Bad_elec_file=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Elec2Interpolate_Database1.xls');
[~,~,raw] = xlsread(Bad_elec_file);
Subj_names=raw(:,1);
Bad_elec=raw(:,2);

trig_type_name{1}='std';
trig_type_name{2}='dev';
trig_type_name{3}='IC1';
trig_type_name{4}='NIC1';

baseline_tp=204;
  
for trig=1:size(trig_type_name,2)
    
    OutDir0=sprintf('%s/%s',OutDir,trig_type_name{trig});
    
    OutDir1=sprintf('%s/%s/epochs',OutDir,trig_type_name{trig});
    
    switch trig_type_name{trig}
        case 'std'
            all_trig={'std'};
        case 'dev'
            all_trig={'dev1','dev2','dev3'};
        case 'IC1'
            all_trig={'IC1'};
        case 'NIC1'
            all_trig={'NIC1'};
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

            if count_epoch
                files=dir(sprintf('%s/%s*.mat',OutDir2,Subj));
                for poc=Npoc:size(files,1)
                    load(sprintf('%s/%s',OutDir2,files(poc-Npoc+1).name));
                    All_pocs(:,:,poc)=double(epoch-mean(epoch(:,1:baseline_tp),2));
                end
                Npoc=poc+1;
            end
            a=1;
        end
        count_ratio=sum(count_epoch)/(sum(count_epoch)+sum(count_rej));
        
        if count_ratio>0.2
            ERP=mean(All_pocs,3);
            figure
            plot(ERP')
            save(sprintf('%s/%s_ERP.mat',OutDir0,Subj),'ERP');
        end
        clear All_pocs count_epoch count_rej
    end
       a=1;     
end

Subj2check1=Subj_names(sum(count_ratio<0.3,2)>=1,:)

disp('Auditory \n')

Subj2check2=Subj_names(sum(count_epoch(:,2:4),2)<100,:)

disp('Visual \n')

Subj2check2=Subj_names(sum(count_epoch(:,5:6),2)<100,:)
