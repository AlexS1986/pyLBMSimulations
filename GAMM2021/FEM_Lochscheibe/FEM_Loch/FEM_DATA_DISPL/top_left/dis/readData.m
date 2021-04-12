clear;
close all;

%LBM
LBM_filename = 'displLBM.dis';
LBM_Data = txt2mat(LBM_filename);

LBM_t = LBM_Data(:,1);
%LBM_u = LBM_Data(:,2);
for i=1:length(LBM_t)
    LBM_u(i) = sqrt(LBM_Data(i,2)^2+LBM_Data(i,3)^2);
end

%FEM
FEM_filename = 'Plochscheibea.dis';
FEM_Data = txt2mat(FEM_filename);

FEM_t = FEM_Data(:,1);
% FEM_u = FEM_Data(:,2:3)
for i=1:length(FEM_t)
    FEM_u(i) = sqrt(FEM_Data(i,2)^2+FEM_Data(i,3)^2);
end
    

%Save 
save(['DataLBM','.mat'],'LBM_t','LBM_u');
save(['DataFEM','.mat'],'FEM_t','FEM_u');

%Plot
co = get(0,'defaultaxescolororder'); 
lso = get(0,'defaultaxeslinestyleorder'); 
%set(0,'defaultaxescolororder',[0 0 0]) %black and gray
%set(0,'defaultaxeslinestyleorder',{':','-','--','-.'}) %or whatever you want
figure('Position', [10 10 1200 1600])
plot(LBM_t, LBM_u, FEM_t,FEM_u,'LineWidth',7.5);
legend('{\it LBM}', '{\it FEM}','Location', 'southeast','FontSize', 50)
xlabel('time $t$','Interpreter','latex','FontSize',50);
ylabel(['displacement $\|\mathbf{u}\|$'],'Interpreter','latex','FontSize',50);
axis([0 2 0 0.0015]);
xticks([0 1 2])
yticks([0.000 0.00025 0.0005 0.00075 0.001 0.00125 0.0015])
%yticks([0.000 0.0005 0.005])
grid on 
set(gca,'FontSize',50)
set(legend,'FontSize',50,'Location','EastOutside','Interpreter','latex')
set(0,'defaultaxescolororder',co); 
set(0,'defaultaxeslinestyleorder',lso);
%hax =axes
