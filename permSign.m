function [permmedian,permmean] = permSign(data)

randomseed_ft(13);
nsub=length(data);
for np=1:100  % number of permutations
  permlabel=round(rand(1,nsub));
  for nt=1:nsub
    if permlabel(nt)
      permtrialvals(nt)=-data(nt);
    else
      permtrialvals(nt)=data(nt);
    end
  end
  permmedian(np)=median(permtrialvals);
  permmean(np)=mean(permtrialvals);
end


