
clear all
close all
c_dir=cd;


%% Future input arguments and global variables

%raw_path='/Volumes/My Passport for Mac/PdM/RawData2/Control';
%raw_path='/Volumes/My Passport for Mac/PdM/RawData2/Early-phase-1stRecording';
raw_path='/Users/mip/Documents/PdM/Data/BDFs/All_new/EEG data';

%addpath(genpath('/Users/laura/Documents/EPFL/Projets_Master/PdM/Code/Matlab/eeglab14_1_1b'));

sr=1024;
low_cf=1;
high_cf=40;
filt_order=8;

%% Preprocessing

Bad_elec={};
count_su=1;

% Patients
files = dir( fullfile(raw_path) );   %# list all *.eph files (eg. Ctrl001.Epoch 001.eph)
files = {files.name}';                      %'# file names

%data = cell(numel(files),1);                %# store file contents
for i=1:numel(files)
    
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
    
        fprintf('Loading data from subject %s (i=%d)\n',Subj,i);
        
%         Bad_elec{count_su,1}=Subj;
%         
%         bdfs = dir( fullfile(sprintf('%s/%s',raw_path,files{i}),'*.bdf'));
%         bdfs = {bdfs.name}';
%         
%         flag_EC=1;
%         flag_MMN=1;
%         flag_Kani=1;
%         
%         count_ta=2;
%         
%         for j=1:numel(bdfs)
%              
%             Key   = sprintf('%s_',Subj);  
%             Index = strfind(bdfs{j}, Key);
%             Task = sscanf(bdfs{j}(1,Index(1) + length(Key):end-4), '%s', 1);
%             fprintf('For task %s \n',Task);
%             if ((Task(1)=='e'||Task(1)=='E'||Task(1)=='R'||Task(1)=='r')&&flag_EC) || (Task(1)=='M'&&flag_MMN) || ((Task(1)=='K'||Task(1)=='k')&&flag_Kani)
%                 filename=sprintf('%s/%s/%s',raw_path,files{i},bdfs{j});
%                 %[dta,The_data.nbchan,labels,txt,The_data.srate,gain,prefiltering,ChanDim] = eeg_read_bdf(filename,'all','n');
%                 data=pop_biosig(filename);
%                 The_data=pop_importdata('setname',['Sub' Subj ...
%                 '_resting_original'],'data',data.data,'dataformat','array','srate',...
%                 data.srate,'nbchan',data.nbchan,'chanlocs',data.chanlocs);
% 
%                 %The_data.data=dta-repmat(mean(dta,2),[1 size(dta(:,:),2)]);
%                 %The_data.data=The_data.data-mean(The_data.data,2);
%                 %The_data.pnts=size(The_data.data,2);
%                 %The_data.xmax=(The_data.pnts-1)/sr;
%                 %The_data.xmin=0;
%                 %The_data.times=linspace(The_data.xmin,The_data.xmax,size(The_data.data,2));
% 
%                 % Data filtering
%                 %disp('Filtering data...');
%                 base_filtered=The_data;
%                 base_filtered.data = double(base_filtered.data(1:64,:));
%     %             base_filtered.data=base_filtered.data-repmat(mean(base_filtered.data,2),...
%     %                 [1 size(base_filtered.data,2)]);
%                 h1=fdesign.highpass('N,F3dB',filt_order,low_cf,sr);
%                 d1=design(h1,'Butter');
%                 h2=fdesign.lowpass('N,F3dB',filt_order,high_cf,sr);
%                 d2=design(h2,'Butter');
%                 %freqz(d1)
%                 %freqz(d2)
%                 for nn=1:size(base_filtered.data,1)
%                     base_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
%                        base_filtered.data(nn,:));
%                     base_filtered.data(nn,:)=filtfilt(d2.sosMatrix,d2.ScaleValues,...
%                         double(base_filtered.data(nn,:)));
%                 end
%     %             for nn=size(base_filtered.data,1)+1:size(The_data.data,1)
%     %                 base_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
%     %                    base_filtered.data(nn,:));
%     %             end
% 
%                 % Channels visual removal
%                 %disp('Channels visual removal...');
%                 eegplot(base_filtered.data,'srate',sr);%,'srate',sr_new
%                 bad_chans=input('Bad channels (vector): ');
%                 good_chans=setdiff(1:size(base_filtered.data,1),bad_chans);
%                 Bad_elec{count_su,count_ta}=bad_chans;
%                 count_ta=count_ta+1;
%                 close all
%                 switch Task(1)
%                     case 'E'
%                         flag_EC=0;
%                     case 'e'
%                         flag_EC=0;
%                     case 'R'
%                         flag_EC=0;
%                     case 'r'
%                         flag_EC=0;
%                     case 'M'
%                         flag_MMN=0;
%                     case 'K'
%                         flag_Kani=0;
%                     case 'k'
%                         flag_Kani=0;
%                 end
%             end
%         end
%         name_mat=sprintf('/Users/laura/Documents/EPFL/Projets_Master/PdM/Data/Interp_bad_elec/Bad_elec_%s.mat',Subj);
%         for m=1:size(Bad_elec,2)
%             bad_el_to_file{m}=Bad_elec{count_su,m};
%         end
%         save(name_mat,'bad_el_to_file')
%         count_su=count_su+1;
%         fprintf('Last subject finished: (i=%d) \n',i)
    
end