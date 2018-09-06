% plot of sleep results

means=[25 24 58 20 2]
stds=[5 3 8 4 1]

figure(1)
bar(means);hold on;
errorbar(means,stds,'o')
ax=get(1);

set(ax.Children,'XTickLabel',{'W', 'N1', 'N2', 'N3', 'REM'})
ylabel('Minutes')
set(ax.Children,'FontSize',26)

