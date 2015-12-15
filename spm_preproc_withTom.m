% spm preproc: convert, concat, montage


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

ii=32;
sleep=1;

cd([edir sub{ii} ])

clearvars -except edir ddir sub ii sleep

if sleep
  fileb=dir('cc*b.mat');
else
  fileb=dir('cc*s.mat');
end

fasst

%%




load(['Mcc_fasst/M' fileb.name])
% end
CRC=D.other.CRC;
CRC.score{1}
sum(CRC.score{1}==2 |  CRC.score{1}==3)
if sleep
  save('CRC_sleep_withTom.mat','CRC');
else
  save('CRC_wake_withTom.mat','CRC');
end
load(fileb.name)
if isfield(D.other,'CRC') % can exist from artifact detection
  D.other.CRC.score=CRC.score;
else
  D.other.CRC=CRC;
end
if ispc
  save(fullfile([edir sub{ii} '\'],D.fname), 'D');
else
  D.path(strfind(D.path,'\'))='/';
  save(fullfile(['/mnt/hgfs/D' D.path(3:end)], D.fname), 'D');
end


%%% END %%%%

%%
% % copy sleep scores into original un-montaged file

if 0
  if ~done_awake_already
    % first the awake bit (which may include subjects falling asleep!)
    if isdir('oldmontage') % will apply to ii==11-15 as these were scored without use of separate ExG channels even though included in data
      load(['oldmontage/M' files.name])
    else
      load(['M' files.name])
    end
    
    CRC=D.other.CRC;
    load(files.name)
    D.other.CRC=CRC;
    if ispc
      save(fullfile(D.path,D.fname), 'D');
    else
      D.path(strfind(D.path,'\'))='/';
      save(fullfile(['/mnt/hgfs/D' D.path(3:end)], D.fname), 'D');
    end
  else
    % then the sleep
    % if isdir('oldmontage') % will apply to ii==11-15 as these were scored without use of separate ExG channels even though included in data
    %   load(['oldmontage/M' fileb.name])
    % else
    load(['M' fileb.name])
    % end
    CRC=D.other.CRC;
    save('CRC_sleep.mat','CRC');
    load(fileb.name)
    if isfield(D.other,'CRC') % can exist from artifact detection
      D.other.CRC.score=CRC.score;
    else
      D.other.CRC=CRC;
    end
    if ispc
      save(fullfile([edir sub{ii} '\'],D.fname), 'D');
    else
      D.path(strfind(D.path,'\'))='/';
      save(fullfile(['/mnt/hgfs/D' D.path(3:end)], D.fname), 'D');
    end
  end
end

%%
if 0
  % sleep montage
  % load([edir 'sleepmontage.mat']) % AASM
  load([edir 'sleepAnalyzerMontage.mat']) % from Jo/Becky
  mkdir analyzer
  copyfile([files.name(1:end-4) '*'],'analyzer')
  copyfile([fileb.name(1:end-4) '*'],'analyzer')
  cd analyzer
  
  S=[];
  S.D=files.name;
  S.montage=montage;
  S.keepothers=0;
  [D,montage]=spm_eeg_montage(S);
  
  S=[];
  S.D=[fileb.name];
  S.montage=montage;
  S.keepothers=0;
  [D,montage]=spm_eeg_montage(S);
  
  %
  % Do sleep staging in FASST in Windows Matlab
  
  % % copy sleep scores into original un-montaged file
  % first the awake bit (which may include subjects falling asleep!)
  load(['M' files.name])
  CRC=D.other.CRC;
  load(files.name)
  D.other.CRC=CRC;
  if ispc
    save(fullfile(D.path,D.fname), 'D');
  else
    D.path(strfind(D.path,'\'))='/';
    save(fullfile(['/mnt/hgfs/D' D.path(3:end)], D.fname), 'D');
  end
  
  % then the sleep
  load(['M' fileb.name])
  CRC=D.other.CRC;
  load(fileb.name)
  D.other.CRC=CRC;
  if ispc
    save(fullfile(D.path,D.fname), 'D');
  else
    D.path(strfind(D.path,'\'))='/';
    save(fullfile(['/mnt/hgfs/D' D.path(3:end)], D.fname), 'D');
  end
  
  cd ..
  
end


% go to eeg_legomagic_preproc2.m



