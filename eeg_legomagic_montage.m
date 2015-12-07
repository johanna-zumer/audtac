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



