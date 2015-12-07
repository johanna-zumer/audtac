
montage.labelorg={'Fp1'    'Fp2'   'F7'    'F8'   'Fz'    'Cz'    'Pz'    'Oz'   'TP9'    'TP10'   'FT9'    'FT10' }

montage.labelnew={'EOG H F7-F8' 'EMG FT9-FT10' 'EOG RV Fp2-FT9' 'EOG LV Fp1-FT10' 'Fz-LM' 'Cz-LM' 'Oz-LM'}

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

save('sleepmontageTeog.mat','montage')
