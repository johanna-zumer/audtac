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

sub{100}='p01'; % ma.a. 03/04/14
sub{1}='e01'; % ab.m. 21/05/14
sub{2}='e02'; % ma.a. 04/06/14
sub{3}='e03'; % ag.m. 10/06/14
sub{4}='e04'; % re.g. 17/06/14
%%%  above is pilot, below is real
sub{5}='e05'; % ma.a. 25/06/14
sub{6}='e06'; % e.u.  01/07/14
sub{7}='e07'; % a.s.  07/07/14
sub{8}='e08'; % k.t.  09/07/14
sub{9}='e09';% d.a.  14/07/14
sub{10}='e10';% k.l.  15/07/14
sub{11}='e11';% ab.m.  16/07/14  % from here on, had EOG/EMG
sub{12}='e12';% b.s.  17/07/14
sub{13}='e13';% d.t.  21/07/14
sub{14}='e14';% f.g.  22/07/14
sub{15}='e15';% r.m.  23/07/14
sub{16}='e16';% t.p.  24/07/14
sub{17}='e17';% t.t.  28/07/14
sub{18}='e18';% k.n.v.  29/07/14
sub{19}='e19';% j.b.  30/07/14
sub{20}='e20';% n.m.  31/07/14
sub{21}='e21';% l.c.  04/08/14
sub{22}='e22';% a.b.  05/08/14
sub{23}='e23';% r.c.
sub{24}='e24';% a.d.
sub{25}='e25';% j.c.
sub{26}='e26';% r.s.
sub{27}='e27';% a.p.
sub{28}='e28';% w.p.
sub{29}='e29';% i.r.
sub{30}='e30';% o.y.l.
sub{31}='e31';% r.b.
sub{32}='e32';% i.f.


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



