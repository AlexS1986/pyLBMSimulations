DataD2Q5 = txt2mat('outputD2Q5/disc_with_hole_w.outloc');
DataD2Q9 = txt2mat('output/disc_with_hole_w.outloc');

figure(1)
plot(DataD2Q5(:,1),DataD2Q5(:,2),DataD2Q9(:,1),DataD2Q9(:,2))
legend('D2Q5', 'D2Q9')
grid on