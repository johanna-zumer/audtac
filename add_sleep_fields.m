function data_sum = add_sleep_fields(data1,data2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fn={'delta' 'kc' 'sw' 'sp_fast' 'sp_slow'};

data_sum={};
data_sum.trialinfo1=data1.trialinfo;
data_sum.trialinfo2=data2.trialinfo;
for ff=1:length(fn)
  for tr=1:length(data1.(fn{ff}))
    data_sum.(fn{ff}){tr}=[data1.(fn{ff}){tr} data2.(fn{ff}){tr}];
  end
end

end

