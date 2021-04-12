%
clear;
%close all;
L = 1.0 
%l_div_2 = 1.0 + 10.0 * 1.0/126.0
%l_div_2 = 1.0 %- 0.5*1.0/42.0
%l_div_2 = 0.98999999999999;
%l_div_2 = 1.00000000000001;
%r = [l_div_2-1.00000000e+00, l_div_2-9.76190476e-01, l_div_2-9.52380952e-01, l_div_2-9.28571429e-01];
r = [L-9.76190476e-01, L-9.52380952e-01, L-9.28571429e-01, L-9.04761905e-01];
r_ch_fine = r;
%r_ch_fine = [L-9.92063492e-01, L-9.84126984e-01, L-9.76190476e-01];
delta_local = 0.0; %1.0/42.0;
r_local_gw = [L-9.88095238e-01+delta_local, L-9.64285714e-01+delta_local, L-9.40476190e-01+delta_local , L-9.16666667e-01+delta_local];

index_macr = 3;
index_local_gw =1;

mu = 1.0;
Ks = 1.0*0.0025*sqrt(pi*L);


head = 'Computations/';
folder = [head,'GuangWu_Global/'];
filename = [folder,'output/','disc_with_crack_w.outloc'];
GuangWuGlobalData = txt2mat(filename);

time_gw_glob = GuangWuGlobalData(:,1)*1.0;
delta_gw_glob = [abs(GuangWuGlobalData(:,2)-GuangWuGlobalData(:,3)), ...
    abs(GuangWuGlobalData(:,4)-GuangWuGlobalData(:,5)), abs(GuangWuGlobalData(:,6)-GuangWuGlobalData(:,7))];
Kfact_gw_glob = 1.0/Ks*mu/4.0*sqrt(2*pi) *[delta_gw_glob(:,1)*sqrt(1.0/r(1)), delta_gw_glob(:,2)*sqrt(1.0/r(2)), delta_gw_glob(:,3)*sqrt(1.0/r(3))];

folder = [head,'Classic/'];
filename = [folder,'output/','disc_with_crack_w.outloc'];
ClassicData = txt2mat(filename);

time_cl_glob = ClassicData(:,1)*1.0;
delta_cl_glob = [abs(ClassicData(:,2)-ClassicData(:,3)), abs(ClassicData(:,4)-ClassicData(:,5)), abs(ClassicData(:,6)-ClassicData(:,7))];
Kfact_cl_glob = 1.0/Ks*mu/4.0*sqrt(2*pi) *[delta_cl_glob(:,1)*sqrt(1.0/r_ch_fine(1)), delta_cl_glob(:,2)*sqrt(1.0/r_ch_fine(2)), delta_cl_glob(:,3)*sqrt(1.0/r_ch_fine(3))];


folder = [head,'GuangWu_Local/'];
filename = [folder,'output/','disc_with_crack_w.outloc'];
GuangWuLocalData = txt2mat(filename);

time_gw_loc = GuangWuLocalData(:,1)*1.0;
delta_gw_loc = [abs(GuangWuLocalData(:,2)-GuangWuLocalData(:,3)), ...
    abs(GuangWuLocalData(:,4)-GuangWuLocalData(:,5)), abs(GuangWuLocalData(:,6)-GuangWuLocalData(:,7))];
sqrt(1.0/r_local_gw(1))
sqrt(1.0/r_local_gw(2))
sqrt(1.0/r_local_gw(3))
Kfact_gw_loc = 1.0/Ks*mu/4.0*sqrt(2*pi) *[delta_gw_loc(:,1)*sqrt(1.0/r_local_gw(1)), delta_gw_loc(:,2)*sqrt(1.0/r_local_gw(2)), delta_gw_loc(:,3)*sqrt(1.0/r_local_gw(3))];


save(['DataGuangWuGlob','.mat'],'time_gw_glob','Kfact_gw_glob');
% save(['DataChopardGlob','.mat'],'time_ch_glob','Kfact_ch_glob');
% save(['DataChopardLoc','.mat'],'time_ch_loc','Kfact_ch_loc');
save(['DataGuangWuLoc','.mat'],'time_gw_loc','Kfact_gw_loc');

figure(100)
%hax =axes
plot(time_gw_glob,Kfact_gw_glob(:,index_macr), time_gw_loc, Kfact_gw_loc(:,index_local_gw), time_cl_glob, Kfact_cl_glob(:,index_macr));
legend('gw glob', 'gw loc', 'classic')
line([2 2],[ 0 2])

