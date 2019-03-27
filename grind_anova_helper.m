function [grind_TPA_MSPN,grind_tacPaud,grind_tacMSpN]=grind_anova_helper(grind_tacPaud_save,grind_tacMSpN_save,ll)

soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

grind_TPA_MSPN=ft_math(cfg,grind_tacPaud_save,grind_tacMSpN_save);
grind_tacPaud=grind_tacPaud_save;
grind_tacMSpN=grind_tacMSpN_save;
if ll>5
  grind_TPA_MSPN.time=grind_TPA_MSPN.time-soades(ll);
  grind_tacPaud.time=grind_tacPaud.time-soades(ll);
  grind_tacMSpN.time=grind_tacMSpN.time-soades(ll);
end
