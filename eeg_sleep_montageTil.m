
montage.labelorg={'Fp1'    'Fp2'    'F3'    'F4'    'C3'    'C4'    'P3'    'P4'    'O1'    'O2'    'F7'    'F8'    'T7'    'T8'    'P7'    'P8'    'Fz'    'Cz'    'Pz'    'Oz'    'FC1'    'FC2'    'CP1'    'CP2'    'FC5'    'FC6'    'CP5'    'CP6'    'TP9'    'TP10'    'POz'    'ECG'    'F1'    'F2'    'C1'    'C2'    'P1'    'P2'    'AF3'    'AF4'    'FC3'    'FC4'    'CP3'    'CP4'    'PO3'    'PO4'    'F5'    'F6'    'C5'    'C6'    'P5'    'P6'    'AF7'    'AF8'    'FT7'    'FT8'    'TP7'    'TP8'    'PO7'    'PO8' 'FT9'    'FT10'    'Fpz'    'CPz'}

montage.labelnew={'HEOG F7-F8' 'EMG FT9-FT10' 'EOG RV Fp2-FT9' 'EOG LV Fp1-FT10' 'Fz-LM' 'Cz-LM' 'Oz-LM'}

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
montage.tra(5,match_str(montage.labelorg,'FT9'))=-.5;
montage.tra(5,match_str(montage.labelorg,'FT10'))=-.5;
montage.tra(6,match_str(montage.labelorg,'Cz'))=1;
montage.tra(6,match_str(montage.labelorg,'FT9'))=-.5;
montage.tra(6,match_str(montage.labelorg,'FT10'))=-.5;
montage.tra(7,match_str(montage.labelorg,'Oz'))=1;
montage.tra(7,match_str(montage.labelorg,'FT9'))=-.5;
montage.tra(7,match_str(montage.labelorg,'FT10'))=-.5;

save('sleepmontageT.mat','montage')
