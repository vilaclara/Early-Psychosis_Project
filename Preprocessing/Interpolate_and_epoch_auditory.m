%Interpolate and epoch

clc
clear all
close all

db=2;

c_dir=cd;

task='auditory';
std = [1 11 65281];
dev1 = [2 12 65282]; % Dev1 = Location
dev2 = [3 13 65283]; % Dev2 = Pitch
dev3 = [4 14 65284]; % Dev3 = Duration

baseline_tp=102;
epoch_tp=512;

OS_lat{1}=0; %xp - nestle
OS_lat{2}=189-244; %win7 - nestle
OS_lat{3}=189-174; %linux - nestle
OS_lat{4}=189-268; %win7 - Cery ancien labo
OS_lat{5}=189-262; %win7 - Cery nouveau labo

OutDir=sprintf('/Users/mip/Documents/PdM/Data/ERPs/Dataset%d/',db);
%% Future input arguments and global variables

raw_path=sprintf('/Users/mip/Documents/PdM/Data/BDFs/Database%d_auditory',db);

addpath(genpath('/Users/laura/Documents/EPFL/Projets_Master/PdM/Code/Matlab/eeglab14_1_1b'));
%eeglab

load('/Users/mip/Documents/PdM/Data/Cap/chanlocs.mat')
for chan=1:size(EEG.chanlocs,2)
    chan_names{chan}=EEG.chanlocs(chan).labels;
end    

% Load bad electrode info

Bad_elec_file=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Elec2Interpolate_Database%d.xls',db);
[~,~,raw] = xlsread(Bad_elec_file);
Subj_names=raw(:,1);
Bad_elec=raw(:,2);
Latency=raw(:,3);

sr=1024;
sr_new=1024;
low_cf=0.1;
high_cf=40;
filt_order=8;

Nelec=64;

%% load Bdf data

start_i=1;

for subj=1:size(Subj_names,1)

    Subj=Subj_names{subj};
    OSlat=OS_lat{Latency{subj}};
    
    files=dir(sprintf('%s/*%s*.bdf',raw_path,Subj_names{subj}));
    files = {files.name}';  
    count_epo(1:4)=1;  % nb of accepted epochs
    count_rej(1:4)=1;  % nb of rejected epochs
    
    for bdf=start_i:numel(files)
        
            fprintf('Loading data from subject %s (bdf=%d)\n',Subj,bdf);

            filename=sprintf('%s/%s',raw_path,files{bdf});
            %[dta,The_data.nbchan,labels,txt,The_data.srate,gain,prefiltering,ChanDim] = eeg_read_bdf(filename,'all','n');
            data=pop_biosig(filename);
            base_original=data;
            
            %% Preprocessing

            base_original.pnts=size(base_original.data,2);
            base_original.xmax=(base_original.pnts-1)/sr_new;
            %base_original.times=linspace(base_original.xmin,base_original.xmax,...
            %    size(base_original.data,2));
            %eegplot(base_original.data(1:64,:),'srate',base_original.srate);

            % Add channel position
            cd(c_dir)
            n_channel_activity = Nelec;
            for ch1 = 1:n_channel_activity
                base_original.chanlocs(ch1,1).labels=chan_names{ch1};
                base_original.chanlocs(ch1,1).X = EEG.chanlocs(1,ch1).X;
                base_original.chanlocs(ch1,1).Y = EEG.chanlocs(1,ch1).Y;
                base_original.chanlocs(ch1,1).Z = EEG.chanlocs(1,ch1).Z;
                base_original.chanlocs(ch1,1).theta = EEG.chanlocs(1,ch1).theta;
                base_original.chanlocs(ch1,1).radius = EEG.chanlocs(1,ch1).radius;
                base_original.chanlocs(ch1,1).sph_theta = EEG.chanlocs(1,ch1).sph_theta;
                base_original.chanlocs(ch1,1).sph_phi = EEG.chanlocs(1,ch1).sph_phi;
                base_original.chanlocs(ch1,1).sph_radius = EEG.chanlocs(1,ch1).sph_radius;
                base_original.chanlocs(ch1,1).type = EEG.chanlocs(1,ch1).type;
            end

            %mkdir(namefolder)
            %cd(namefolder);
            %save('base_original.mat','base_original','-v7.3');
            
            bad_chans=[Nelec+1:size(base_original.data,1)];
            base_ch_removed1=base_original;
            good_chans=setdiff(1:Nelec,bad_chans);
            n_channels=length(good_chans);
            base_ch_removed1.data=base_ch_removed1.data(good_chans,:);
            base_ch_removed1.nbchan=n_channels;
            base_ch_removed1.chanlocs=base_ch_removed1.chanlocs(good_chans);
            %base_ch_removed1.setname=['Sub' subject '_resting_ch_removed1'];
            %save('base_ch_removed1.mat','base_ch_removed1');

            % Data filtering
            disp('Filtering data...');
            base_filtered=base_ch_removed1;
            base_filtered.data = double(base_filtered.data);
            base_filtered.data=base_filtered.data-repmat(mean(base_filtered.data,2),...
                [1 size(base_filtered.data,2)]);
            h1=fdesign.highpass('N,F3dB',filt_order,low_cf,sr);
            d1=design(h1,'Butter');
            h2=fdesign.lowpass('N,F3dB',filt_order,high_cf,sr);
            d2=design(h2,'Butter');
            % Plot response
            %freqz(d1)
            %freqz(d2)
            for nn=1:n_channel_activity
                base_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
                   base_filtered.data(nn,:));
                base_filtered.data(nn,:)=filtfilt(d2.sosMatrix,d2.ScaleValues,...
                    double(base_filtered.data(nn,:)));
            end
            for nn=n_channel_activity+1:size(base_ch_removed1.data,1)
                base_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
                   base_filtered.data(nn,:));
            end
            base_filtered.setname=['Sub' Subj '_' task '_filtered'];
            %save('base_filtered.mat','base_filtered','-v7.3');
            %clear base_ch_removed1

            base_resampled=base_filtered;
            
            % Data resampling
%             disp('Resampling data...');
%             base_resampled=pop_resample(base_filtered,sr_new);        
%             base_resampled.setname=['Sub' subject '_resting_resampled'];          
%             save('base_resampled.mat','base_resampled','-v7.3');
            clear base_filtered

            
            % Channels removal
            bad_chans=[];
            chain=Bad_elec{subj};
            if isnan(chain)
                splitted{1}=[];
            else
                splitted = textscan(chain,'%s','Delimiter',' ');
            end
            for i=1:size(splitted{1},1)
                chan=splitted{1}{i};
                flag=0;
                if ~isempty(chan)
                    for j=1:size(chan_names,2)
                    	if strcmp(chan_names{j},chan)
                            bad_chans=[bad_chans j];
                            flag=1;
                        end
                    end
                    if ~flag
                        disp('Electrode not found!!!');
                    end
                end
            end    
            base_ch_removed2=base_resampled;    
            good_chans=setdiff(1:size(base_ch_removed2.data,1),bad_chans);
            n_channels=length(good_chans);
            base_ch_removed2.data=base_ch_removed2.data(good_chans,:);
            base_ch_removed2.nbchan=n_channels;
            base_ch_removed2.chanlocs=base_ch_removed2.chanlocs(good_chans);
            base_ch_removed2.setname=['Sub' Subj '_' task '_ch_removed2'];
            %eegplot(base_ch_removed2.data(1:n_channel_activity-length(bad_chans),:),'srate',sr_new);
            %save('base_ch_removed2.mat','base_ch_removed2');

            % Channel interpolation
            base_ch_interp = base_resampled;
            base_ch_interp = pop_interp(base_ch_interp,bad_chans,'spherical');
            base_ch_interp.setname=['Sub' Subj '_' task '_interp'];
            %eegplot(base_ch_interp.data(1:n_channel_activity,:),'srate',sr_new);
            %save('base_ch_interp.mat','base_ch_interp');
            
            %% EPOCH

            trigs=base_ch_interp.event;
            for i=1:size(trigs,2)
                T_type(i)=trigs(i).type;
                T_lat(i)=trigs(i).latency-OSlat;
            end
            if ~isempty(trigs)
                [trig_types,~,which_trig]=unique(T_type);
                for i=1:size(trig_types,2)
                    nb_trig(i)=sum(which_trig==i);
                end

                if sum(trig_types)<11
                    id=1;
                elseif sum(trig_types)<60
                    id=2;
                else
                    id=3;
                end
                
                trig_type_name{1}='std';
                trig_type_name{2}='dev1';
                trig_type_name{3}='dev2';
                trig_type_name{4}='dev3';
                trig_type(1)=std(id);
                trig_type(2)=dev1(id);
                trig_type(3)=dev2(id);
                trig_type(4)=dev3(id);

                for trig=1:size(trig_type,2)
                    fprintf('%s = %s \n',trig_type_name{trig},files{bdf});
                    OutDir1=sprintf('%s/%s/epochs',OutDir,trig_type_name{trig});
                    latencies=(T_lat(:,T_type==trig_type(trig)));
                    latencies=latencies(latencies>baseline_tp);
                    latencies=latencies(latencies<(size(base_ch_interp.data,2)-epoch_tp));
                    OutDir2=sprintf('%s/%s/epochs/%s',OutDir,trig_type_name{trig},Subj);
                    if ~exist(OutDir2) 
                        cd(OutDir1)
                        mkdir(Subj)
                        OutDir2=sprintf('%s/%s/epochs/%s',OutDir,trig_type_name{trig},Subj);
                        cd(OutDir2)
                    else
                        cd(OutDir2)
                        load('nb_epochs.mat')
                        count_epo(trig)=count_e;
                        count_rej(trig)=count_r;
                    end
                    epochs=zeros(Nelec,baseline_tp+epoch_tp+1);
                    epochs_CAR=epochs;
                    N=0;
                    for Nlat=1:size(latencies,2)
                        lat=latencies(Nlat);
                        epoch=base_ch_interp.data(:,lat-baseline_tp:lat+epoch_tp);
%                         maxE=max(epoch(:,102:716)');
%                         minE=min(epoch(:,102:716)');
                        maxE=max(epoch');
                        minE=min(epoch');
                        if (maxE-minE)<80
                            epochs=epochs+epoch;
                            epochs_CAR=epochs_CAR+(epoch-mean(epoch(:,1:baseline_tp),2));
                            %save(sprintf('%s/%s_epoch_%d.mat',OutDir2,Subj,count_epo(trig)),'epoch');
                            save(sprintf('%s/%s_epoch_%d_%d.mat',OutDir2,Subj,bdf,lat),'epoch');
                            count_epo(trig)=count_epo(trig)+1;
                            N=N+1;
                        else
                            count_rej(trig)=count_rej(trig)+1;
                        end
                    end
                    epochs=epochs/N;
                    epochs_CAR=epochs_CAR/N;
    %                 figure
    %                 plot(epochs_CAR','LineWidth',2)
                    count_e=count_epo(trig);
                    count_r=count_rej(trig);
                    save('nb_epochs.mat','count_e','count_r')
                    %save(sprintf('%s/%s_ERP.mat',OutDir1,Subj),'epochs');
                    %save(sprintf('%s/%s_ERP_CAR.mat',OutDir1,Subj),'epochs_CAR');
                    clear OutDir1 OutDir2 latencies epochs epochs_CAR epoch lat Nlat
                end
                
                if count_epo(1)<10
                    disp('not enough trials')
                end
            end
            
            movefile(filename,[raw_path '/Done']);
            
            clearvars -except epoch_tp baseline_tp Latency OS_lat OSlat count_epo count_rej bdf subj c_dir files Subj start_i Nelec sr filt_order high_cf low_cf sr_new task std dev1 dev2 dev3 OutDir raw_path chan_names EEG Subj_names Bad_elec
            close all
    end
end
