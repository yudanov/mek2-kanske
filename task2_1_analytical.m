clc; clear all;% clf;
%% Define constants
R1 = 1.5e3;%m
R2 = 2.5e3;%m
m = 5.0e13;%kg http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=1989AC
omega_0 = [20,32,98].*pi/(180*24*60*60); % initial values, rad/s
prec_time_period = 7.35;%precessions periodtid dygn

%% calculated values
r_bar = (R2^4-R1^4)/(R1^3+R2^3);%dist from sphere tangent to total body CoM
m1 = R1^3*m/(R1^3+R2^3);
m2 = R2^3*m/(R1^3+R2^3);
% Moments of inertia
I_zeta = 2/5*m1*R1^2 + 2/5*m2*R2^2;
I_0 = I_zeta + m1*(R1+r_bar)^2 + m2*(R2 - r_bar)^2;

A = (I_0-I_zeta)/I_0; % gamma_xi
B = (I_zeta-I_0)/I_0; % gamma_eta
gamma_zeta = 0; % 


% ODE solution constants
%C = omega_0(3); % solution to omega_zeta
lambda_1 = sqrt(A*B)*omega_0(3);
lambda_2 = -sqrt(A*B)*omega_0(3);
C1 = (omega_0(2) + omega_0(1)*sqrt(B/A))/(2*sqrt(B/A));
C2 = omega_0(1) - C1;

omega_xi  = @(t) C1*exp(lambda_1.*t) + C2*exp(lambda_2.*t);
omega_eta = @(t) C1*exp(lambda_1.*t)*1i - C2*exp(lambda_2.*t)*1i;
omega_zeta = @(t) omega_0(3)+0.*t; % t vektor, utan +0*t returnerar skalär :-(


t = linspace(0,100,1000).*24*60^2;
omega = [omega_xi(t);
        omega_eta(t);
        omega_zeta(t)];
omegaAbs = sqrt(sum(omega.^2));
    
I = [I_0 0 0;
    0 I_0 0;
    0 0 I_zeta];

L = I*omega;
LAbs = sqrt(sum(L.^2));

T = 1/2 * omegaAbs .* LAbs;
%T2 = 1/2 * omega'*(I*omega);
T2 = zeros(1,length(L));
T3 = zeros(1,length(L));
for i = 1:length(L)
    T2(i) = 1/2*omega(:,i)'*I*omega(:,i);
   % T3(i) = 1/2 *(I_0*omega(1,i)^2 + I_0*omega(2,i)^2 + I_zeta*omega(3,i)^2);
end

% nånting om precession
% calc_temp = sqrt(sum((I*omega).^2))./((I_0-I_zeta)*omega(3,:));
% psi_temp = [calc_temp;calc_temp;calc_temp];
% psi_dot = omega.*psi_temp;

%% Plot that fucker!
t_korr = 1/(3600*24); %s->dygn
omega_korr = 24*3600*180/pi; %rad/s -> grader/dygn
omegaPlot = omega*omega_korr;
omegaAbsPlot = omegaAbs*omega_korr;
tPlot = t*t_korr;
figure(1)
plot(tPlot,omegaPlot(1,:),'-k', ...
    tPlot,omegaPlot(2,:),'--k',...
    tPlot,omegaPlot(3,:),'-.k',...
    tPlot,omegaAbsPlot,':k',...
    'LineWidth',2);
    legend('\omega_\xi','\omega_\eta','\omega_\zeta','\omega')
    title('Analytisk lösning, symmetrisk kropp. Rotationsvektorns komponenter \omega_{\xi,\eta,\zeta} för Toutatis 4179')
    xlabel('Tid [Dygn]')
    ylabel('Vinkelfrekvens [ ^\circ / dygn]')
    axis([t(1)*t_korr t(end)*t_korr -40 120])


figure(2)
LPlot = L;
LAbsPlot = LAbs;
plot(tPlot,LPlot,...
    tPlot,LAbsPlot,...
    'LineWidth',2);
    legend('L_\xi','L_\eta','L_\zeta','|\bf{L}|')
    title('Analytisk lösning, symmetrisk kropp. Rörelsemängdsmomentets komponenter L_{\xi,\eta,\zeta} för Toutatis 4179')
    xlabel('Tid [Dygn]')
    ylabel('Rörelsemängdsmoment [ kg \cdot m^2 \cdot rad/s]')

figure(2)

plot(tPlot,T,...
    tPlot,T2,...
    'LineWidth',2);
    legend('T','T2')
    title('Analytisk lösning, symmetrisk kropp. Röreleseenergi')
    xlabel('Tid [Dygn]')
    ylabel('Energi [J]')
    
% Save plots
%saveas(1,'analytic_omega_task2_1','png')