% montage


clear all
if ispc
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
else
  edir='/mnt/hgfs/D/audtac/eeg_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end
cd(edir)

% sub{100}='p01'; 
sub{1}='e01'; 
sub{2}='e02'; 
sub{3}='e03'; 
sub{4}='e04'; 
%%%  above is pilot, below is real
sub{5}='e05'; 
sub{6}='e06'; 
sub{7}='e07'; 
sub{8}='e08'; 
sub{9}='e09';
sub{10}='e10';
sub{11}='e11';
sub{12}='e12';
sub{13}='e13';
sub{14}='e14';
sub{15}='e15';
sub{16}='e16';
sub{17}='e17';
sub{18}='e18';
sub{19}='e19';
sub{20}='e20';
sub{21}='e21';
sub{22}='e22';
sub{23}='e23';
sub{24}='e24';
sub{25}='e25';
sub{26}='e26';
sub{27}='e27';
sub{28}='e28';
sub{29}='e29';
sub{30}='e30';
sub{31}='e31';
sub{32}='e32';

%%
% for ii=[3 5:32]
for ii=30
%   for sleep=[0 1]
  for sleep=[0]
    clearvars -except ii sleep sub edir ddir
    
    cd([edir sub{ii} ])
    
    try % if this code has been run already, then don't over-write the cc* files
      % there should only be 1 file for each if successfully removed intermediates
      files=dir('cc*s.mat');
      fileb=dir('cc*b.mat');
    catch
      if ii==1
        files=dir('cc*.mat');
        fileb=dir('spm8*b.mat');
      else
        % load into SPM format
        sfiles=dir('*s.eeg');
        bfiles=dir('*b.eeg');
        
        if ii==1
          sfiles(1)=dir('*12.eeg');
          sfiles(2:4)=dir('*s.eeg');
        end
        
        for ll=1:length(sfiles)
          S.dataset=sfiles(ll).name;
          S.continuous=1;
          D=spm_eeg_convert(S);
        end
        
        for ll=1:length(bfiles)
          S.dataset=bfiles(ll).name;
          S.continuous=1;
          D=spm_eeg_convert(S);
        end
        
        if 0 % concatenate (automate for multiple files later?)
          sfiles=dir('spm8_e*s.mat');
          for ll=1:floor(length(sfiles)/2)
            prefile=[sfiles(ll).name; sfiles(ll+1).name]
            crc_concatenate(prefile)
          end
        else
          % Do file concat in FASST in Windows Matlab
          
          % for both sitting(s) and bed(b) portions
          
          % remove unncessary intermediate files
        end
      end
    end
    
    if ii<11
      load('D:\audtac\eeg_data\sleepMontage_combo_noExG.mat');
    else
      load('D:\audtac\eeg_data\sleepMontage_combo_ExG.mat');  % still to be made
    end
    
    S=[];
    if sleep
      S.D=fileb.name;
    else
      if ~exist('CRC_wake.mat')
        load(files.name)
        CRC=D.other.CRC; % save out what was previously done, even though done using sleepmontageT(xG)
        save('CRC_wake.mat','CRC');
        clear CRC
      else
        movefile('CRC_wake.mat','CRC_wake_old.mat')
      end
      S.D=files.name;
    end
    
    S.montage=montage;
    S.keepothers=0;
    [D,montage]=spm_eeg_montage(S);
    
    if ~isdir('Mcc_fasst')
      mkdir('Mcc_fasst')
    end
    movefile(['M' S.D(1:end-3) '*'],'Mcc_fasst');
    
  end
end



