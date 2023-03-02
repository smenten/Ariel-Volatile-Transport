[ Tsurf ] = thermal_ariel_f(90) ;

%% plot surface temps together
fonts = 15;
fontn = 'Helvetica';
% fontw = 'bold';
% 
% set(gca, 'fontname', fontn, 'fontsize', fonts, 'fontweight', fontw, 'linewidth', 1);
% set(gcf, 'Color', 'w');
% set(gcf,'Units','centimeters');
% set(gcf,'outerposition',[2 20 25 17]); % left bottom width height

load('n90LatSurfaceTemp.mat', 'lsWrapped', 'Tsurf')   ;
whereCrossOver360to0 = find(floor(lsWrapped)==0,1); % First index after Ls wrapped from 360 back to 0
subplot(2,1,1)
plot(lsWrapped(1:whereCrossOver360to0-1), Tsurf(1:whereCrossOver360to0-1), 'color', [0 0.4470 0.7410] );
hold on
plot(lsWrapped(whereCrossOver360to0:end), Tsurf(whereCrossOver360to0:end), 'color', [0 0.4470 0.7410]);
ylabel('Temperature (K)', 'fontname', fontn, 'fontsize', fonts, 'fontweight', 'normal');


load('0LatSurfaceTemp.mat', 'lsWrapped', 'Tsurf')   ;
whereCrossOver360to0 = find(floor(lsWrapped)==0,1); % First index after Ls wrapped from 360 back to 0
subplot(2,1,1)
plot(lsWrapped(1:whereCrossOver360to0-1), Tsurf(1:whereCrossOver360to0-1), 'color', [0.8500 0.3250 0.0980]);

plot(lsWrapped(whereCrossOver360to0:end), Tsurf(whereCrossOver360to0:end), 'color', [0.8500 0.3250 0.0980]);
ylabel('Temperature (K)', 'fontname', fontn, 'fontsize', fonts, 'fontweight', 'normal');


load('90LatSurfaceTemp.mat', 'lsWrapped', 'Tsurf')   ;
whereCrossOver360to0 = find(floor(lsWrapped)==0,1); % First index after Ls wrapped from 360 back to 0
subplot(2,1,1)
plot(lsWrapped(1:whereCrossOver360to0-1), Tsurf(1:whereCrossOver360to0-1), 'color', [0.9290 0.6940 0.1250]);

plot(lsWrapped(whereCrossOver360to0:end), Tsurf(whereCrossOver360to0:end), 'color', [0.9290 0.6940 0.1250]);
ylabel('Temperature (K)', 'fontname', fontn, 'fontsize', fonts, 'fontweight', 'normal');
hold off 



xlim([0 360]);


%%

figure
x = 1:5;
s1 = scatter(x,[6 3 9 10 7],100,"filled");
hold on
s2 = scatter(x,[16 13 19 20 17],100,"filled");
s3 = scatter(x,[26 23 29 33 27],100,"filled");
hold off


