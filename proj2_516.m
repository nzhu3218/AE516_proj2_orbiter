function proj2_516

close all;

aoa = [10 20 30];    % AoAs collected by CFD
M = 3.8;                % mach number

cl = NaN(size(aoa));    % preallocation
cd = NaN(size(aoa));
% cm = NaN(size(aoa));

for i = 1:length(aoa)
%     [cl(i),cd(i),cm(i)] = readAoA(aoa(i));
    [cl(i),cd(i)] = readAoA(aoa(i));
end

[aoa2,cl2,cd2] = Dilao_and_Fonseca(aoa,M);  % fit of experimental data

graph_comparison(aoa,cl,cd,aoa2,cl2,cd2); % plotting

end

function [avg_cl,avg_cd] = readAoA(AoA)

cd = importdata([num2str(AoA), '/report-drag-rfile.out']);
cl = importdata([num2str(AoA), '/report-lift-rfile.out']);
% cm = importdata([num2str(AoA), '/report-mom-rfile.out']);

cd = cd.data(:,2);
cl = cl.data(:,2);
% cm = cm.data(:,2);

n = 300; % averaging over the last 300 points or however many is available

if numel(cd) < n
    avg_cd = sum(cd)/numel(cd);
    avg_cl = sum(cl)/numel(cl);
%     avg_cm = -sum(cm)/numel(cm);
else
    avg_cd = sum(cd(end-n:end))/(n+1);
    avg_cl = sum(cl(end-n:end))/(n+1);
%     avg_cm = -sum(cm(end-n:end))/(n+1);
end

% % shifting moment about leading edge to quarter chord
% avg_cm = avg_cm + avg_cl/4;

% plotting the convergence history

iter = 1:numel(cd);

figure;
plot(iter,cd,iter,cl);
title('Convergence History'); grid on; grid minor;
xlabel('Iterations '); ylabel('Aerodynamic Coefficient [-]');
legend('Cd','Cl');
set(gca,'fontname','times','fontsize',16);
axis([0 max(iter) min([cd;cl])/1.2 max([cd;cl])*1.2])

end

% Journal of aerospace enigneering. ISSN 0893-1321/04015012(12)
function [aoa_deg,cl,cd] = Dilao_and_Fonseca(aoa, M)

maxAoA = max(aoa);
minAoA = min(aoa);

aoa_deg = linspace(minAoA,maxAoA,20);
aoa = aoa_deg*pi/180;

a1 = -0.053; a2 = 2.73; a3 = -1.55;
b1 = -1.01; b2 = 1.1;
d3 = 1.79;
e1 = -1.4; e2 = 1.5;
f1 = 0.028; f2 = 1.4;
Mc = 1.25;

K = 1/2*(1 + sqrt(abs(1-(M/Mc)^2)));

cl = (a1 + a2*aoa + a3*aoa.^2).*K.^(b1 + aoa*b2);

cd = (0.01 + f1 * M^f2 + d3 * aoa.^2).*K.^(e1 + aoa*e2);

end

function graph_comparison(aoa,cl,cd,aoa2,cl2,cd2)

% error_plot = 0; % 1 shows error bar
% drag_error_plot = 0; % 1 to plot just the drag
% forces_plot = 0;

% Cl and Cd comparison
coeff = 1;
if coeff
    figure;
    plot(aoa2,cl2,'r',aoa2,cd2,'b'); hold on;
    plot(aoa,cl,'r+', aoa, cd,'b+','markersize',12); hold off;
    title('Aerodynamic Coefficients of the Space Shutte');
    grid on; grid minor;
    xlabel('AoA [deg]'); ylabel('Aerodynamic Coefficients [-]');
    set(gca,'fontname','times','fontsize',16);
    axis([0 50, 0 1.2]);
    legend('Cl (Experimental)','Cd (Experimental)', ...
        'Cl (CFD)','Cd (CFD)','location','nw');
end

% lift to drag ratio

L2D = 1;
if L2D
    
    l2d_anal = cl2./cd2;
    l2d_cfd = cl./cd;
    
    % matlab built-in polyfit
    p = polyfit(aoa,l2d_cfd,2);
    l2d_cfd_fit = polyval(p,aoa2);

    figure;
    plot(aoa2,l2d_anal,'r'); hold on;
    plot(aoa,l2d_cfd,'b*','markersize',12);
    plot(aoa2,l2d_cfd_fit,'b--'); hold off;
    
    title('Lift to drag ratio of the Space Shutte');
    grid on; grid minor;
    xlabel('AoA [deg]'); ylabel('L/D [-]');
    set(gca,'fontname','times','fontsize',16);
    axis([min(aoa)-5, max(aoa)+5, min([l2d_cfd,l2d_anal])/1.1,...
        max([l2d_cfd,l2d_anal])*1.1]);
    legend('L/D (Experimental)','L/D (CFD)','2nd-order fit of L/D (CFD)');
end


% if error_plot
%     
%     % 10% error bar plots for Cl and Cd
%     err_cl = cl2 * 0.1;
%     err_cd = cd2 * 0.1;
%     errorbar(aoa2,cl2,err_cl); hold on;
%     errorbar(aoa2,cd2,err_cd)
%     plot(aoa,cl,'-o',aoa,cd,'-o','linewidth',3); hold off;
%     grid on; grid minor;
%     xlabel('AoA [deg]'); ylabel('Aerodynamic Coefficients [-]');
%     title('C_l and C_d vs. AoA with error bars');
%     legend('C_l Xfoil','C_d Xfoil','C_l CFD','C_d CFD','Location','nw');
%     set(gca,'fontname','times','fontsize',18)
%     
% elseif drag_error_plot
% 
%     figure;
%     err_cd = cd2 * 0.1;
%     errorbar(aoa2,cd2,err_cd); hold on;
%     plot(aoa,cd,'-o','linewidth',3); hold off;
%     xlabel('AoA [deg]'); grid on; grid minor;   
%     ylabel('Coefficients of drag, C_d [-]');
%     title('Coefficients of drag vs. AoA');
%     
%     legend('C_d CFD','C_d Xfoil','Location','nw');
%     set(gca,'fontname','times','fontsize',18)
% 
%     
% elseif forces_plot
%     
%     rho = 1.225;    % [kg/m^3]
%     V = 17.88;      % [m/s]
%     q = 1/2*rho*V^2;    % [Pa]
%     c = 1;              % chord length
%     
%     % dimensionalization
%     L = cl*q*c;
%     D = cd*q*c;
%     M = cm*q*c^2;
%     
%     figure;
%     plot(aoa,L,'-o',aoa,D,'-o',aoa,M,'-o','linewidth',3);
%     
%     xlabel('AoA [deg]'); grid on; grid minor;
%     ylabel('Aerodynamic Forces and Moments [N, N.m]');
%     title('Aerodynamic Forces and Moments vs. AoA');
%     
%     legend('Lift','Drag','Moment, c/4','Location','nw');
%     set(gca,'fontname','times','fontsize',18)
% 
%     
%     
% else
%     
%     % aerodynamic coefficients plotting
%     figure;
%     
%     subplot(1,2,1);
%     plot(aoa,cl,'-o',aoa2,cl2,aoa,cd,'-o',aoa2,cd2, ...
%                     aoa,cm,'-o',aoa2,cm2,'linewidth',3);
%     
%     xlabel('AoA [deg]'); grid on; grid minor;
%     ylabel('Aerodynamic Coefficients [-]');
%     title('Aerodynamic Coefficients vs. AoA');
%     
%     legend('C_l CFD','C_l Xfoil','C_d CFD','C_d Xfoil','C_{m,c/4} CFD',...
%                                         'C_{m,c/4} Xfoil','Location','nw');
%     set(gca,'fontname','times','fontsize',18)
% 
%     subplot(1,2,2)
%     plot(aoa,cl./cd,'-o',aoa2,cl2./cd2,'linewidth',3);
%     xlabel('AoA [deg]'); ylabel('L/D ratio [-]');
%     grid on; grid minor;
%     title('L/D ratio vs. AoA');
%     legend('L/D CFD','L/D Xfoil','location','nw');
%     set(gca,'fontname','times','fontsize',18)
% 
% end

end