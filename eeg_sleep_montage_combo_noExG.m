
montage.labelorg={'Fp1'    'Fp2'    'F3'    'F4'    'C3'    'C4'    'P3'    'P4'    'O1'    'O2'    'F7'    'F8'    'T7'    'T8'    'P7'    'P8'    'Fz'    'Cz'    'Pz'    'Oz'    'FC1'    'FC2'    'CP1'    'CP2'    'FC5'    'FC6'    'CP5'    'CP6'    'TP9'    'TP10'    'POz'    'ECG'    'F1'    'F2'    'C1'    'C2'    'P1'    'P2'    'AF3'    'AF4'    'FC3'    'FC4'    'CP3'    'CP4'    'PO3'    'PO4'    'F5'    'F6'    'C5'    'C6'    'P5'    'P6'    'AF7'    'AF8'    'FT7'    'FT8'    'TP7'    'TP8'    'PO7'    'PO8' 'FT9'    'FT10'    'Fpz'    'CPz'}

% montage.labelnew={'HEOG F7-F8' 'EMG FT9-FT10' 'EOG RV Fp2-FT9' 'EOG LV Fp1-FT10' 'Fz-LM' 'Cz-LM' 'Oz-LM'    'Cz-C3' 'C4-Cz' 'Cz-Pz' 'Cz-Fz' 'F7-F3' 'F4-F8' 'Fp1-F7' 'Fp2-F8'   }
montage.labelnew={'EOG H F7-F8' 'EMG FT9-FT10' 'EOG RV Fp2-FT9' 'EOG LV Fp1-FT10' 'Fz-LM' 'Cz-LM' 'Oz-LM'    'Cz-C3' 'C4-Cz' 'Cz-Pz' 'Cz-Fz'  'Fz-F3' 'F4-Fz'}

montage.tra=zeros(length(montage.labelnew),length(montage.labelorg))

montage.tra(1,match_str(montage.labelorg,'F7'))=1;
montage.tra(1,match_str(montage.labelorg,'F8'))=-1;
montage.tra(2,match_str(montage.labelorg,'FT9'))=1;
montage.tra(2,match_str(montage.labelorg,'FT10'))=-1;
montage.tra(3,match_str(montage.labelorg,'Fp2'))=1;
montage.tra(3,match_str(montage.labelorg,'FT9'))=-1;
montage.tra(4,match_str(montage.labelorg,'Fp1'))=1;
montage.tra(4,match_str(montage.labelorg,'FT10'))=-1;

montage.tra(5,match_str(montage.labelorg,'Fz'))=1;
montage.tra(5,match_str(montage.labelorg,'TP9'))=-.5;
montage.tra(5,match_str(montage.labelorg,'TP10'))=-.5;
montage.tra(6,match_str(montage.labelorg,'Cz'))=1;
montage.tra(6,match_str(montage.labelorg,'TP9'))=-.5;
montage.tra(6,match_str(montage.labelorg,'TP10'))=-.5;
montage.tra(7,match_str(montage.labelorg,'Oz'))=1;
montage.tra(7,match_str(montage.labelorg,'TP9'))=-.5;
montage.tra(7,match_str(montage.labelorg,'TP10'))=-.5;

montage.tra(8,match_str(montage.labelorg,'Cz'))=1;
montage.tra(8,match_str(montage.labelorg,'C3'))=-1;
montage.tra(9,match_str(montage.labelorg,'C4'))=1;
montage.tra(9,match_str(montage.labelorg,'Cz'))=-1;
montage.tra(10,match_str(montage.labelorg,'Cz'))=1;
montage.tra(10,match_str(montage.labelorg,'Pz'))=-1;
montage.tra(11,match_str(montage.labelorg,'Cz'))=1;
montage.tra(11,match_str(montage.labelorg,'Fz'))=-1;

montage.tra(12,match_str(montage.labelorg,'Fz'))=1;
montage.tra(12,match_str(montage.labelorg,'F3'))=-1;
montage.tra(13,match_str(montage.labelorg,'F4'))=1;
montage.tra(13,match_str(montage.labelorg,'Fz'))=-1;

% montage.tra(12,match_str(montage.labelorg,'F7'))=1;
% montage.tra(12,match_str(montage.labelorg,'F3'))=-1;
% montage.tra(13,match_str(montage.labelorg,'F4'))=1;
% montage.tra(13,match_str(montage.labelorg,'F8'))=-1;
% montage.tra(14,match_str(montage.labelorg,'Fp1'))=1;
% montage.tra(14,match_str(montage.labelorg,'F7'))=-1;
% montage.tra(15,match_str(montage.labelorg,'Fp2'))=1;
% montage.tra(15,match_str(montage.labelorg,'F8'))=-1;


% montage.tra(2,match_str(montage.labelorg,'Fp1'))=1;
% montage.tra(2,match_str(montage.labelorg,'F3'))=-1;
% montage.tra(3,match_str(montage.labelorg,'Fp2'))=1;
% montage.tra(3,match_str(montage.labelorg,'F4'))=-1;
% montage.tra(5,match_str(montage.labelorg,'F7'))=1;
% montage.tra(5,match_str(montage.labelorg,'T7'))=-1;
% montage.tra(6,match_str(montage.labelorg,'F3'))=1;
% montage.tra(6,match_str(montage.labelorg,'C3'))=-1;
% montage.tra(7,match_str(montage.labelorg,'F4'))=1;
% montage.tra(7,match_str(montage.labelorg,'C4'))=-1;
% montage.tra(8,match_str(montage.labelorg,'F8'))=1;
% montage.tra(8,match_str(montage.labelorg,'T8'))=-1;
% montage.tra(9,match_str(montage.labelorg,'C3'))=1;
% montage.tra(9,match_str(montage.labelorg,'T7'))=-1;
% montage.tra(12,match_str(montage.labelorg,'T8'))=1;
% montage.tra(12,match_str(montage.labelorg,'C4'))=-1;
% montage.tra(13,match_str(montage.labelorg,'T7'))=1;
% montage.tra(13,match_str(montage.labelorg,'P7'))=-1;
% montage.tra(14,match_str(montage.labelorg,'C3'))=1;
% montage.tra(14,match_str(montage.labelorg,'P3'))=-1;
% montage.tra(15,match_str(montage.labelorg,'C4'))=1;
% montage.tra(15,match_str(montage.labelorg,'P4'))=-1;
% montage.tra(16,match_str(montage.labelorg,'T8'))=1;
% montage.tra(16,match_str(montage.labelorg,'P8'))=-1;
% montage.tra(17,match_str(montage.labelorg,'P7'))=1;
% montage.tra(17,match_str(montage.labelorg,'O1'))=-1;
% montage.tra(18,match_str(montage.labelorg,'P3'))=1;
% montage.tra(18,match_str(montage.labelorg,'O1'))=-1;
% montage.tra(19,match_str(montage.labelorg,'P4'))=1;
% montage.tra(19,match_str(montage.labelorg,'O2'))=-1;
% montage.tra(20,match_str(montage.labelorg,'P8'))=1;
% montage.tra(20,match_str(montage.labelorg,'O2'))=-1;

save('D:\audtac\eeg_data\sleepMontage_combo_noExG.mat','montage')
