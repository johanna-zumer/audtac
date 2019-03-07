num_sub = 5;
num_cond = 7;

y = rand(num_cond*num_sub, 1)  % fake data


Subject = nominal(reshape(repmat([1 2 3 4 5 6 7],num_sub,1),num_sub*num_cond,1));
G = {Subject};
Z = {ones(num_sub*num_cond,1)};

x_center_peri = [[1 2 3 4 3 2 1]-mean([1 2 3 4 3 2 1])]'; % pre-recal
x_Lin = [-3 -2 -1 0 1 2 3]'; % pre-recal

x_center_peri_all_sub = repmat(x_center_peri,num_sub,1);
x_Lin = repmat(x_Lin,num_sub,1);

m = 1; % model number

my_model(m).X = [x_Lin x_center_peri_all_sub  ones(num_sub*num_cond,1)]
my_model(m).num_par = size(my_model(m).X,2);
figure; imagesc(my_model(m).X);

lme = fitlmematrix(my_model(m).X,y,Z,G,'FitMethod','REML','FixedEffectPredictors',...
    {'Lin','CentrePeri','Intercept'},...
    'RandomEffectPredictors',{{'Intercept'}},'RandomEffectGroups',{'Subject'});
lme

my_model(m).lme = lme
my_model(m).BIC = my_model(m).lme.ModelCriterion.BIC;  % the better the model, the smaller BIC
my_model(m).AIC =  my_model(m).lme.ModelCriterion.AIC; % the better the model, the smaller AIC
my_model(m).loglike =  my_model(m).lme.ModelCriterion.LogLikelihood; % loglike




