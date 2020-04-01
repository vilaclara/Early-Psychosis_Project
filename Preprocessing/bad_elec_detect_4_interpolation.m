
clear all
close all
c_dir=cd;


%% Future input arguments and global variables

%raw_path='/Volumes/My Passport for Mac/PdM/RawData2/Control';
%raw_path='/Volumes/My Passport for Mac/PdM/RawData2/Early-phase-1stRecording';
raw_path='/Users/mip/Documents/PdM/Data/BDFs/All_new/EEG data';

addpath(genpath('/Users/laura/Documents/EPFL/Projets_Master/PdM/Code/Matlab/eeglab14_1_1b'));

load('/Users/mip/Documents/PdM/Data/Cap/chanlocs.mat')
for i=1:size(EEG.chanlocs,2)
    chan_names{i}=EEG.chanlocs(i).labels;
end    

sr=1024;
low_cf=1;
high_cf=40;
filt_order=8;

%% Preprocessing

name_mat=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Bad_elec_new_subj.mat');
try
    load(name_mat);
    start_i=i;
catch
    warning('No previous file, we start from scatch');
    Bad_elec=cell(1,2);
    start_i=1;
end

count_su=0;
last_subj='';

% Patients
files = dir( fullfile(raw_path) );   %# list all *.eph files (eg. Ctrl001.Epoch 001.eph)
files = {files.name}';                      %'# file names

%data = cell(numel(files),1);                %# store file contents
for i=start_i:numel(files)
    
    if files{i}(1)=='S' || files{i}(1)=='s'
        Key1   = sprintf('PAT');  
        Key2   = sprintf('Pat');
        Key3   = sprintf('CTRL');
        Key4   = sprintf('Ctrl');   % We want the number right after the key
        Index1 = strfind(files{i}, Key1);
        Index2 = strfind(files{i}, Key2);
        Index3 = strfind(files{i}, Key3);
        Index4 = strfind(files{i}, Key4);
        Index=[Index1 Index2 Index3 Index4];
        Subj = sscanf(files{i}(1,Index(1) + length(Key):Index(1) + length(Key)+3), '%s', 1);
    elseif files{i}(1)=='C' || files{i}(1)=='c'
        Key   = sprintf('CONTROL');
        Index = strfind(files{i}, Key);
        Subj0 = sscanf(files{i}(1,Index(1) + length(Key):Index(1) + length(Key)+1), '%s', 1);
        Subj=sprintf('Ctrl0%s',Subj0);
    elseif files{i}(1)=='L' || files{i}(1)=='l'
        Key   = sprintf('LNAC');
        Index = strfind(files{i}, Key);
        Subj0 = sscanf(files{i}(1,Index(1) + length(Key):Index(1) + length(Key)+1), '%s', 1);
        Subj=sprintf('Lnac0%s',Subj0);
    else
        continue
    end
    
        if ~strcmp(last_subj,Subj)
            count_su=count_su+1; 
            Bad_elec{count_su,1}=Subj;
        end
        fprintf('Loading data from subject %s (i=%d)\n',Subj,i);
        
        filename=sprintf('%s/%s',raw_path,files{i});
        %[dta,The_data.nbchan,labels,txt,The_data.srate,gain,prefiltering,ChanDim] = eeg_read_bdf(filename,'all','n');
        data=pop_biosig(filename);
        The_data=pop_importdata('setname',['Sub' Subj ...
        '_resting_original'],'data',data.data,'dataformat','array','srate',...
        data.srate,'nbchan',data.nbchan,'chanlocs',data.chanlocs);

        %The_data.data=dta-repmat(mean(dta,2),[1 size(dta(:,:),2)]);
        %The_data.data=The_data.data-mean(The_data.data,2);
        %The_data.pnts=size(The_data.data,2);
        %The_data.xmax=(The_data.pnts-1)/sr;
        %The_data.xmin=0;
        %The_data.times=linspace(The_data.xmin,The_data.xmax,size(The_data.data,2));

        % Data filtering
        %disp('Filtering data...');
        base_filtered=The_data;
        base_filtered.data = double(base_filtered.data(1:64,:));
%             base_filtered.data=base_filtered.data-repmat(mean(base_filtered.data,2),...
%                 [1 size(base_filtered.data,2)]);
        h1=fdesign.highpass('N,F3dB',filt_order,low_cf,sr);
        d1=design(h1,'Butter');
        h2=fdesign.lowpass('N,F3dB',filt_order,high_cf,sr);
        d2=design(h2,'Butter');
        %freqz(d1)
        %freqz(d2)
        for nn=1:size(base_filtered.data,1)
            base_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
               base_filtered.data(nn,:));
            base_filtered.data(nn,:)=filtfilt(d2.sosMatrix,d2.ScaleValues,...
                double(base_filtered.data(nn,:)));
        end
%         for nn=size(base_filtered.data,1)+1:size(The_data.data,1)
%             base_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
%                base_filtered.data(nn,:));
%         end

        % Channels visual removal
        %disp('Channels visual removal...');
        eegplot(base_filtered.data,'srate',sr);%,'srate',sr_new
        bad_chans=input('Bad channels (vector): ');
        good_chans=setdiff(1:size(base_filtered.data,1),bad_chans);
        names=Bad_elec{count_su,2};
        for j=1:size(bad_chans,2)
            names=sprintf('%s %s(%d)',names,chan_names{bad_chans(j)},bad_chans(j));
        end
        Bad_elec{count_su,2}=names;
        
        
        close all
        last_subj=Subj;
    
        % Save file just in case
        name_mat=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Bad_elec_new_subj.mat');
        save(name_mat,'Bad_elec','i')
end

name_mat=sprintf('/Users/mip/Documents/Early-psychosis_Project/Preprocessing/Bad_elec_new_subj.mat');
save(name_mat,'Bad_elec')
fprintf('Last subject finished: (i=%d) \n',i)