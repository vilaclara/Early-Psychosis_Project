%Interpolate and epoch

clc
clear all
close all

c_dir=cd;

task='std';

%% Future input arguments and global variables

raw_path='/Users/mip/Documents/PdM/Data/BDFs/Database1_auditory';

addpath(genpath('/Users/laura/Documents/EPFL/Projets_Master/PdM/Code/Matlab/eeglab14_1_1b'));
%eeglab

load('/Users/mip/Documents/PdM/Data/Cap/chanlocs.mat')
for chan=1:size(EEG.chanlocs,2)
    chan_names{chan}=EEG.chanlocs(chan).labels;
end    

% Load bad electrode info

Bad_elec_file=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Elec2Interpolate_Database1.xls');
[~,~,raw] = xlsread(Bad_elec_file);
Subj_names=raw(:,1);
Bad_elec=raw(:,2);

sr=1024;
sr_new=1024;
low_cf=1;
high_cf=40;
filt_order=8;


%% load Bdf data

start_i=1;

for subj=1:size(Subj_names,1)

    Subj=Subj_names{subj};
    files=dir(sprintf('%s/%s*.bdf',raw_path,Subj_names{subj}));
    files = {files.name}';  
    
    for bdf=start_i:numel(files)
        
            fprintf('Loading data from subject %s (bdf=%d)\n',Subj,bdf);

            filename=sprintf('%s/%s',raw_path,files{bdf});
            %[dta,The_data.nbchan,labels,txt,The_data.srate,gain,prefiltering,ChanDim] = eeg_read_bdf(filename,'all','n');
            data=pop_biosig(filename);
            base_original=data;
            
            %% Preprocessing

            
            %% To add for the "rest" files
            base_original.pnts=size(base_original.data,2);
            base_original.xmax=(base_original.pnts-1)/sr_new;
            %base_original.times=linspace(base_original.xmin,base_original.xmax,...
            %    size(base_original.data,2));
            eegplot(base_original.data(1:64,:),'srate',base_original.srate);

            % Add channel position
            cd(c_dir)
            n_channel_activity = 64;
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
            
            bad_chans=[65:size(base_original.data,1)];
            base_ch_removed1=base_original;
            good_chans=setdiff(1:64,bad_chans);
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
            freqz(d1)
            freqz(d2)
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
            splitted = textscan(chain,'%s','Delimiter',' ');
            for i=1:size(splitted{1},1)
                chan=splitted{1}{i};
                if ~isempty(chan)
                    for j=1:size(chan_names,2)
                    	if strcmp(chan_names{j},chan)
                            bad_chans=[bad_chans j];
                        end
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
            eegplot(base_ch_removed2.data(1:n_channel_activity-length(bad_chans),:),'srate',sr_new);
            %save('base_ch_removed2.mat','base_ch_removed2');

            % Channel interpolation
            base_ch_interp = base_resampled;
            base_ch_interp = pop_interp(base_ch_interp,bad_chans,'spherical');
            base_ch_interp.setname=['Sub' Subj '_' task '_interp'];
            eegplot(base_ch_interp.data(1:n_channel_activity,:),'srate',sr_new);
            %save('base_ch_interp.mat','base_ch_interp');
            
            %% EPOCH
            
            
            
            %% Common average re-referencing
            disp('Re-referencing data...');
            base_rereferenced=base_ch_interp;
            base_rereferenced.data(1:n_channel_activity,:)=base_rereferenced.data(1:n_channel_activity,:)-repmat(mean...
                (base_rereferenced.data(1:n_channel_activity,:),1),[size(base_rereferenced.data(1:n_channel_activity,:),1) 1]);
            base_rereferenced.setname=['Sub' subject '_resting_rereferenced'];
            eegplot(base_rereferenced.data(1:n_channel_activity,:),'srate',sr_new);
            %save('base_rereferenced.mat','base_rereferenced');
            channels_eyes = [1:32 34 35 36 37 40];

            % Identify part where the patient opens the eye
            TMPREJ = [];
            eegplot(base_rereferenced.data(channels_eyes,:),'srate',sr_new,'command','close');
            rm_data=[];
            for nn=1:size(TMPREJ,1)
                rm_data=[rm_data round(TMPREJ(nn,1)):round(TMPREJ(nn,2))];
            end

            
    end
end
