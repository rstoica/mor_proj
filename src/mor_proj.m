% Main Project Script
% Date: 30.11.2014
% Author: Andrei Stoica

% system parameters:
% R1 = R13 = 1/100
% R3 = R5 = ... = R11 = 100
% R2 = R4 = ... = R14 = 2
% R15 = 3
% R16 = R17 = 4

clc;
clear all;
% inputs of system considered
m = 8;
% states of system considered
n = 17;
% outputs of system considered
p = 7;

sys = init_dyn_sys(m,n,p,'Original');

% generate the A matrix (see project report for details)
for i = 1:sys.n
    if (i==1 || i==13)
        sys.A(i,i) = -100;
    elseif (mod(i,2) == 0)
        sys.A(i,i) = -2;
    elseif (mod(i,2) == 1)
        sys.A(i,i) = -1/100;
    end
    if (i+1 <= sys.n-3)
        sys.A(i,i+1) = -1;
    end
    if (i+1 <= sys.n-3)
        sys.A(i+1,i) = 1;
    end
end
sys.A(15,15) = -3;
sys.A(16,16) = -4;
sys.A(17,17) = -4;
sys.A(15,1) = -1;
sys.A(16,3) = 1;
sys.A(17,5) = 1;
sys.A(17,11) = -1;
sys.A(16,13) = -1;
sys.A(1,15) = 1;
sys.A(3,16) = -1;
sys.A(13,16) = 1;
sys.A(5,17) = -1;
sys.A(11,17) = 1;

% generate the B matrix (see project report for details)
for i = 1:2:13
    sys.B(i,floor(i/2)+1) = -1;
    sys.B(i,end) = 2;
end

% generate the C matrix (see project report for details)
sys.C = sys.B(:,1:7)';
% get poles of the system
sys = anal_sys(sys);

tol1 = 10^-2;
tol1 = 0.5*10^-2;
tol1 = 10^-3;

% model reduction using balancing transformation (using Grammians) 
[sys_bal, res_bal] = bal_reduce(sys, tol1);

% model reduction using modal approximation (diagonalization of A matrix)
[sys_mod, res_mod] = mod_reduce(sys, tol1);

% plot of the normalized Hankel singular values (normalization done with
% respect to the highest sv) and of the modal coefficients (normalization 
% done with respect to the higher one)
svfig = figure(1);
hold on
semilogy(res_bal.hsv_norm,'b*-','MarkerSize',5);
semilogy(res_mod.coeff_norm,'rD-','MarkerSize',5);
set(gca,'yscale','log');
title('Semilog plot of normalized Hankel SV');
xlabel('dimension k');
% ylabel('$\displaystyle\frac{\sigma_{k+1}}{\sigma_1}$','interpreter','latex');
ylabel('reduction normalized indicators');
legend('Hankel normaliz. sv','Modal coeff.');
xlim([1 sys.n]);
grid on

% poles of the systems
poles_plot = [real(sys.poles),imag(sys.poles)];
bal_poles_plot = [real(sys_bal.poles),imag(sys_bal.poles)];
mod_poles_plot = [real(sys_mod.poles),imag(sys_mod.poles)];
pfig = figure(2);
hold on
grid on 
plot(poles_plot(:,1),poles_plot(:,2),'ko','MarkerSize',15,'LineWidth',2,'MarkerFaceColor','b');
plot(bal_poles_plot(:,1),bal_poles_plot(:,2),'ko','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','r');
plot(mod_poles_plot(:,1),mod_poles_plot(:,2),'ko','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','y');
legend('Original Sys','Bal. Trunc. System','Mod. Trunc. System','Location','NorthWest');
plot(zeros(1,50),linspace(-2,2,50),'k--');
plot(linspace(-100,0.5,50),zeros(1,50),'k--');
xlim([-100 0.5]);
ylim([-2 2]);
title('Poles of the original system vs. poles of approx. systems');
xlabel('real axis');
ylabel('imaginary axis');
pause()
xlim([-4 0.5]);

t = linspace(0,40,10000)';
f = logspace(-3,3,60);
% step response and Bode plot of the systems
% original system
[sys, sol_sys] = solve_sys(sys,t,f);
[sys_bal, sol_sys_bal] = solve_sys(sys_bal,t,f);
[sys_mod, sol_sys_mod] = solve_sys(sys_mod,t,f);

error_norms.bal_hinf = norm(sol_sys.hinf-sol_sys_bal.hinf);
error_norms.bal_h2 = norm(sol_sys.h2-sol_sys_bal.h2);
error_norms.mod_hinf = norm(sol_sys.hinf-sol_sys_mod.hinf);
error_norms.mod_h2 = norm(sol_sys.h2-sol_sys_mod.h2);
error_norms.bal_hinf_bound = 2*sum(res_bal.hsv(res_bal.k+1:end));

disp(error_norms)

% % error systems plots: step and Bode plot
% stfig_err = figure(6);
% for i = 1:in
%     subplot(2,4,i), plot(t,y_o(:,:,i)-y_m(:,:,i),'b-');
%     hold on
%     plot(t,y_o(:,:,i)-y_b(:,:,i),'r-');
%     title(['Error in step response for input ',inputs(i),' at outputs']);
%     xlabel('seconds')
%     ylabel('amplitude (V)')
% end
% 
% bfig_err = figure(10);
% hold on
% semilogx(bf,20*log10(bamp(1:size(bamp,1),:))-20*log10(bamp_bal(1:size(bamp_bal,1),:)),'b');
% semilogx(bf,20*log10(bamp(1:size(bamp,1),:))-20*log10(bamp_mod(1:size(bamp_mod,1),:)),'r');
% grid on
% set(gca,'xscale','log');
% title('Bode plot of error systems')
% xlabel('Frequency (rad/s)');
% ylabel('Amplitude (dB)');