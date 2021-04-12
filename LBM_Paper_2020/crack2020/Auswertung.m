clear;
close all;

factor = 1.0 % 1.4962/0.8359;

rix_LBM=3;
rix_LBMnew_glob =3;
rix_LBMnew_loc =3;

rix_FEM_notch=3;
rix_FEM_crack=3;


% load('DataLBM1.mat');
% timeLBM100=timeLBM;
% KfactLBM100=KfactLBM;
% clear timeLBM KfactLBM;

load('DataLBM336_1.mat');
timeLBM336_100=timeLBM;
KfactLBM336_100=KfactLBM;
clear timeLBM KfactLBM;


load('DataAn.mat');

load('DataFEM_crack_fine.mat');
timeFEMcrack=timeFEM;
KfactFEMcrack=KfactFEM;
clear timeFEM KfactFEM;


load('DataGuangWuGlob.mat') % 'time_gw_glob','Kfact_gw_glob'
%load('DataChopardGlob.mat') % 'time_ch_glob','Kfact_ch_glob'
%load('DataChopardLoc.mat')  % 'time_ch_loc','Kfact_ch_loc'
load('DataGuangWuLoc.mat')  % 'time_gw_loc','Kfact_qw_loc'





%load('DataFEM.mat');


%co = get(0,'defaultaxescolororder'); 
%lso = get(0,'defaultaxeslinestyleorder'); 


% figure(4)
 %set(0,'defaultaxescolororder',[0 0 0]) %black and gray
 %set(0,'defaultaxeslinestyleorder',{':','-','--','-.'}) %or whatever you want

figure(1)
plot(timeAN,KfactAN, ...
    timeFEMcrack,KfactFEMcrack(:,rix_FEM_crack),...
    timeLBM336_100,KfactLBM336_100(:,rix_LBM), ... 
    time_gw_glob/2, Kfact_gw_glob(:,rix_LBMnew_glob)*factor, ...
    time_gw_loc/2, Kfact_gw_loc(:,rix_LBMnew_loc)*factor,'LineWidth',2.5); 
    % 'LineWidth',2.5) %, ...
   % time_gw_glob/2, Kfact_gw_glob(:,2)*factor)

%legend('analytic',['FEM crack ',num2str(rix_FEM_crack),'$\Delta{x}$'],'LBM $\tau=1.00$')
legend('analytic','FEM crack','LBM mesh-conforming', 'LBM non-mesh conforming', 'LBM non-mesh conforming local')


% 



xlabel('$t$ [2a/$c_s$]','Interpreter','latex','FontSize',45);
ylabel([' $K_{III}$ [$K_s$] '],'Interpreter','latex','FontSize',45);
%title(['Distance from crack tip ',num2str(rix_LBM),'$\Delta{x}$'],'Interpreter','latex','FontSize',45)
axis([0,2,0,2])
axis square
set(gca,'FontSize',30)
set(legend,'FontSize',30,'Location','EastOutside','Interpreter','latex')

% set(0,'defaultaxescolororder',co); 
%set(0,'defaultaxeslinestyleorder',lso);

% figure(2);
% plot(timeAN,KfactAN,timeFEMcrack,KfactFEMcrack(:,rix_FEM_crack),'LineWidth',2)
% legend('analytic',['FEM notch ',num2str(rix_FEM_notch),'$\Delta{x}$'],['FEM crack ',num2str(rix_FEM_crack),'$\Delta{x}$'],'LBM 1 coarse','LBM 2 coarse','LBM 3 coarse','LBM 2 finer','LBM 3 finer');
% axis([0,2,0,2])
% axis square
% set(legend,'FontSize',30,'Location','EastOutside','Interpreter','latex')
% grid on;
% xlabel('$t$ [2a/$c_s$]','Interpreter','latex','FontSize',30);
% ylabel(' $K$ [$K_s$]','Interpreter','latex','FontSize',30);
% set(gca,'FontSize',30)  
% 
% 
% figure(3)
% plot(timeAN,KfactAN, ...
%     timeFEMcrack,KfactFEMcrack,...
%     timeLBM100,KfactLBM100)
% legend('analytic','FEM crack 1$\Delta{x}$','FEM crack 2$\Delta{x}$','FEM crack 3$\Delta{x}$','FEM crack 4$\Delta{x}$','LBM 1$\Delta{x}$','LBM 2$\Delta{x}$','LBM 3$\Delta{x}$','LBM 4$\Delta{x}$','LBM 5$\Delta{x}$','LBM 6$\Delta{x}$','LBM 7$\Delta{x}$')
% xlabel('$t$ [2L/$c_s$]','Interpreter','latex','FontSize',30);
% ylabel([' $K_{III}$ [$K_s$] '],'Interpreter','latex','FontSize',30);
% axis([0,2,0,2])
% axis square
% set(legend,'FontSize',30,'Location','EastOutside','Interpreter','latex')
% title(['$\tau=1.0$ '],'Interpreter','latex','FontSize',30)
% 
% 
% co = get(0,'defaultaxescolororder'); 
% lso = get(0,'defaultaxeslinestyleorder'); 
% 
% figure(4)
% set(0,'defaultaxescolororder',[0 0 0]) %black and gray
% set(0,'defaultaxeslinestyleorder',{':','-','--','-.'}) %or whatever you want
% 
% % subplot(131)
% % plot(timeAN,KfactAN,'r', ...
% %      timeLBM336_100,KfactLBM336_100(:,1:4),'LineWidth',1.5)
% % set(gca,'FontSize',30)
% % xlabel('$t$ [$2L$/$c_s$]','Interpreter','latex','FontSize',30);
% % ylabel([' $K$ [$K_s$] '],'Interpreter','latex','FontSize',30);
% % %set(legend,'FontSize',24,'Location','EastOutside','Interpreter','latex');
% % %legend('analytic','1\Delta{x}','2\Delta{x}','3\Delta{x}','4\Delta{x}','5\Delta{x}','6\Delta{x}','7\Delta{x}','Interpreter','latex')
% % %title('LBM')
% % axis([0,2,0,2])
% 
% 
% subplot(132)
% plot(timeAN,KfactAN,'r','LineWidth',1.5)
% set(gca,'FontSize',30)
% xlabel('$t$ [$2L$/$c_s$]','Interpreter','latex','FontSize',30);
% ylabel([' $K$ [$K_s$] '],'Interpreter','latex','FontSize',30);
% %legend('analytic','1\Delta{x}','2\Delta{x}','3\Delta{x}','4\Delta{x}','Interpreter','latex')
% %title('FEM notch')
% %set(legend,'FontSize',24,'Location','EastOutside','Interpreter','latex');
% axis([0,2,0,2])
% 
% subplot(133)
% plot(timeAN,KfactAN)
% set(gca,'FontSize',30)
% xlabel('$t$ [$2L$/$c_s$]','Interpreter','latex','FontSize',30);
% ylabel([' $K$ [$K_s$] '],'Interpreter','latex','FontSize',30);
% legend('analytic','1\Delta{h}','2\Delta{h}','3\Delta{h}','4\Delta{h}','Interpreter','latex')
% %title('FEM crack')
% axis([0,2,0,2])
% %set(legend,'FontSize',24,'Location','EastOutside','Interpreter','latex');


% maxAN=max(KfactAN);
% maxLBM=max(KfactLBM336_100(:,4));
% maxFEMnotch=max(KfactFEMnotch(:,4));
% maxFEMcrack=max(KfactFEMcrack(:,4));
% relErrorLBM=(maxLBM-maxAN)/(maxAN);
% relErrorFEMnotch=(maxFEMnotch-maxAN)/(maxAN);
% relErrorFEMcrack=(maxFEMcrack-maxAN)/(maxAN);
% 
% tau=0:0.01:2;
% 
% 
% for i=1:length(tau)
%     R(i)=(tau(i)^2-tau(i)+1/6)/(tau(i)-0.5);
% end
% 
% figure(5)
% %subplot(121)
% %plot(tau,R,'k','LineWidth',2)
% %xlabel('$\tau$','Interpreter','latex','FontSize',30);
% %ylabel([' $\bar{R}$  '],'Interpreter','latex','FontSize',30);
% %axis square
% %set(gca,'FontSize',30)
% %grid on
% 
% 
% %subplot(122)
% 
% set(0,'defaultaxescolororder',[0 0 0; 0.5 0.5 0.5]) %black and gray
% plot(timeAN,KfactAN,'r', ...
%     timeLBM090,KfactLBM090(:,rix_LBM), ...
%     timeLBM095,KfactLBM095(:,rix_LBM), ...
%     timeLBM100,KfactLBM100(:,rix_LBM), ...
%     timeLBM120,KfactLBM120(:,rix_LBM), ...
%     timeLBM140,KfactLBM140(:,rix_LBM), ...
%     timeLBM160,KfactLBM160(:,rix_LBM), ...
%     timeLBM180,KfactLBM180(:,rix_LBM), ...
%     timeLBM200,KfactLBM200(:,rix_LBM),'LineWidth',2.5)
% 
% legend('analytic','$\hat{\tau}=0.90$','$\hat{\tau}=0.95$','$\hat{\tau}=1.00$','$\hat{\tau}=1.20$', ...
%        '$\hat{\tau}=1.40$','$\hat{\tau}=1.60$','$\hat{\tau}=1.80$','$\hat{\tau}=2.00$')
% 
% xlabel('$t$ [$2L$/$c_s$]','Interpreter','latex','FontSize',30);
% ylabel([' $K$ [$K_s$] '],'Interpreter','latex','FontSize',30);
% %title(['Distance from crack tip ',num2str(rix_LBM),'$\Delta{x}$'],'Interpreter','latex','FontSize',45)
% axis([0,2,0,2])
% axis square
% set(gca,'FontSize',30)
% set(legend,'FontSize',30,'Location','EastOutside','Interpreter','latex')
% 
% set(0,'defaultaxescolororder',co); 
% set(0,'defaultaxeslinestyleorder',lso);


