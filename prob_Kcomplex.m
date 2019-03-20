function [pKpred,pKdiff]=prob_Kcomplex(m_n,m_a,m_t,m_at)
%-------------------------------------------------------------------------
% Boolean approach: Noisy-Or model
% -- correct!

% 
%  m_n = P(null event) % probability of spontaneous  k-complex for null condition
%  m_a = 1 - (1-m_n)*(1-sel_a) % sel_a is selective impact of auditory; m_a
%               % probability of k-complex for A condition (either spontaneously or A evoked, but you can't have two k-complexes in this window) 
%  m_t = 1 - (1-m_n)*(1-sel_t) % sel_t is the selective impact of v
% 
% Then, the null model is that 
% m_at = 1 - (1-m_n)*(1-sel_a)*(1-sel_t)
%      = 1 - (1-m_a)*(1-m_t)/(1-m_n)
% 

pKpred = 1 - ((1-m_a)*(1-m_t)/(1-m_n));

pKdiff = m_at-K_pred;


