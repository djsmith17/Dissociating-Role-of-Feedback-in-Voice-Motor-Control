function drawJNDResults(allRunData)

figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',14,...
    'FontName','Arial');
box(axes1,'on');
hold(axes1,'all');

plot(UD.x,'LineWidth',3);
hold on;
plot(find(UD.response==1),UD.x(find(UD.response==1)),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',10);
plot(find(UD.response==0),UD.x(find(UD.response==0)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',10);
line([0 length(UD.response)], [meanJND meanJND],'LineStyle', '-.', 'LineWidth',3,'color',[1 0 1])
xlabel('Trials','FontSize',20,'FontName','Arial');
ylabel('Perturbations','FontSize',20,'FontName','Arial');

end