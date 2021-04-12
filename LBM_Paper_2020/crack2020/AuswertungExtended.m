clear;
close all;

rix_LBM=3;
rix_FEM_notch=3;
rix_FEM_crack=3;




load('DataLBM1.mat');
timeLBM100=timeLBM;
KfactLBM100=KfactLBM;
clear timeLBM KfactLBM;




load('DataLBM336_1.mat');
timeLBM336_100=timeLBM;
KfactLBM336_100=KfactLBM;
maxLBMnotch = max(KfactLBM336_100(:,3));
clear timeLBM KfactLBM;


load('DataAn.mat');
maxAn = max(KfactAN);



% load('DataFEM_notch.mat');
% timeFEMnotch=timeFEM;
% KfactFEMnotch=KfactFEM;
% clear timeFEM KfactFEM;

load('DataFEM_crack_fine.mat');
timeFEMcrack=timeFEM;
KfactFEMcrack=KfactFEM;
maxFEM = max(KfactFEMcrack(:,3));

clear timeFEM KfactFEM;



load('DataGuangWuGlob.mat') % 'time_gw_glob','Kfact_gw_glob''
maxGuangWuGlob = max(Kfact_gw_glob(:,3));
load('DataGuangWuLoc.mat')  % 'time_gw_loc','Kfact_qw_loc'
maxGuangWuLoc = max(Kfact_gw_loc(:,3));

%load('DataFEM.mat');


relFEM = abs(maxFEM-maxAn)/maxAn
relLBMnotch=abs(maxLBMnotch-maxAn)/maxAn
relLBMGuangWuGlob = abs(maxGuangWuGlob-maxAn)/maxAn
relLBMGuangWuLoc = abs(maxGuangWuLoc-maxAn)/maxAn


fontsize=20;  


co = get(0,'defaultaxescolororder'); 
lso = get(0,'defaultaxeslinestyleorder'); 

figure(4)
set(0,'defaultaxescolororder',[0 0 0]) %black and gray
set(0,'defaultaxeslinestyleorder',{'-',':','--','-'}) %or whatever you want

subplot(141)
plot(timeAN,KfactAN,'r', ...
     timeLBM336_100,KfactLBM336_100(:,1:3),'LineWidth',1.5)
set(gca,'FontSize',fontsize)
xlabel('$t$ [$2L$/$c_s$]','Interpreter','latex','FontSize',fontsize);
ylabel([' $K$ [$K_s$] '],'Interpreter','latex','FontSize',fontsize);
%set(legend,'FontSize',24,'Location','EastOutside','Interpreter','latex');
%legend('analytic','1\Delta{x}','2\Delta{x}','3\Delta{x}','4\Delta{x}','5\Delta{x}','6\Delta{x}','7\Delta{x}','Interpreter','latex')
%title('LBM')
axis([0,2,0,2])


subplot(142)
plot(timeAN,KfactAN,'r', ...
     time_gw_glob/2, Kfact_gw_glob,'LineWidth',1.5)
set(gca,'FontSize',fontsize)
xlabel('$t$ [$2L$/$c_s$]','Interpreter','latex','FontSize',fontsize);
ylabel([' $K$ [$K_s$] '],'Interpreter','latex','FontSize',fontsize);
%legend('analytic','1\Delta{x}','2\Delta{x}','3\Delta{x}','4\Delta{x}','Interpreter','latex')
%title('FEM notch')
%set(legend,'FontSize',24,'Location','EastOutside','Interpreter','latex');
axis([0,2,0,2])

subplot(143)
plot(timeAN,KfactAN,'r', ...
     time_gw_loc/2, Kfact_gw_loc,'LineWidth',1.5)
set(gca,'FontSize',fontsize)
xlabel('$t$ [$2L$/$c_s$]','Interpreter','latex','FontSize',fontsize);
ylabel([' $K$ [$K_s$] '],'Interpreter','latex','FontSize',fontsize);
%legend('analytic','1\Delta{x}','2\Delta{x}','3\Delta{x}','4\Delta{x}','Interpreter','latex')
%title('FEM notch')
%set(legend,'FontSize',24,'Location','EastOutside','Interpreter','latex');
axis([0,2,0,2])


subplot(144)
plot(timeAN,KfactAN,'r', ...
     timeFEMcrack,KfactFEMcrack(:,1:3),'LineWidth',1.5)
set(gca,'FontSize',fontsize)
xlabel('$t$ [$2L$/$c_s$]','Interpreter','latex','FontSize',fontsize);
ylabel([' $K$ [$K_s$] '],'Interpreter','latex','FontSize',fontsize);
legend('analytic','1\Delta{h}','2\Delta{h}','3\Delta{h}','Interpreter','latex')
%title('FEM crack')
axis([0,2,0,2])
%set(legend,'FontSize',24,'Location','EastOutside','Interpreter','latex');





