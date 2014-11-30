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

% model reduction using balancing transformation (using Grammians) 
tol = 10^-2;
[sys_bal, res_bal] = bal_reduce(sys, tol);

% model reduction using modal approximation (diagonalization of A matrix)
[sys_mod, res_mod] = mod_reduce(sys, tol);

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

% poles of the system
poles_plot = [real(sys.poles),imag(sys.poles)];
mins = [min(poles_plot(:,1)),min(poles_plot(:,2))];
maxs = [max(poles_plot(:,1)),max(poles_plot(:,2))];
bal_poles_plot = [real(sys_bal.poles),imag(sys_bal.poles)];
mins = [mins;min(bal_poles_plot(:,1)),min(bal_poles_plot(:,2))];
maxs = [maxs;max(bal_poles_plot(:,1)),max(bal_poles_plot(:,2))];
mod_poles_plot = [real(sys_mod.poles),imag(sys_mod.poles)];
mins = [mins;min(mod_poles_plot(:,1)),min(mod_poles_plot(:,2))];
maxs = [maxs;max(mod_poles_plot(:,1)),max(mod_poles_plot(:,2))];
pfig = figure(2);
hold on
grid on 
plot(poles_plot(:,1),poles_plot(:,2),'ko','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','b');
plot(bal_poles_plot(:,1),bal_poles_plot(:,2),'ko','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','r');
plot(mod_poles_plot(:,1),mod_poles_plot(:,2),'ko','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','y');
if (max(maxs(:,1)) > 0.5)
    max_x = max(maxs(:,1));
else
    max_x = 0.5;
end
legend('Original Sys','Bal. Trunc. System','Mod. Trunc. System','Location','NorthWest');
plot(zeros(1,50),linspace(min(mins(:,2)),max(maxs(:,2)),50),'k--');
plot(linspace(min(mins(:,1)),max_x,50),zeros(1,50),'k--');
xlim([min(mins(:,1)) max_x]);
ylim([min(mins(:,2)) max(maxs(:,2))]);
pause()
xlim([-4 max_x]);

% step response of the original system
t = linspace(0,40,12500)';
[y_o,~,x_o] = step_sys(sys,t);
in = size(y_o,3);
inputs = {['1=I1'],['2=I3'],['3=I5'],['4=I7'],['5=I9'],['6=I11'],['7=I13'],['8=VDD']};
stfig_o = figure(3);
title('Step response of the original system');
hold on
for i = 1:in
    subplot(2,4,i), plot(t,y_o(:,:,i));
    legend('VC1','VC3','VC5','VC7','VC9','VC11','VC13')
    title(['Unit step response for input ',inputs(i)]);
    xlabel('seconds')
    ylabel('amplitude (V)')
end

% step response of the balanced truncated system
[y_b,t_b,x_b] = step_sys(sys_bal,t);
in = size(y_b,3);
stfig_b = figure(4);
title('Step response of the balanced reduced system');
hold on
for i = 1:in
    subplot(2,4,i), plot(t,y_b(:,:,i));
    legend('VC1','VC3','VC5','VC7','VC9','VC11','VC13')
    title(['Unit step response for input ',inputs(i)]);
    xlabel('seconds')
    ylabel('amplitude (V)')
end

% step response of the modal approximated system
[y_m,~,x_m] = step_sys(sys_mod,t);
in = size(y_m,3);
inputs = {['1=I1'],['2=I3'],['3=I5'],['4=I7'],['5=I9'],['6=I11'],['7=I13'],['8=VDD']};
stfig_c = figure(5);
title('Step response of the modal approximated system');
hold on
for i = 1:in
    subplot(2,4,i), plot(t,y_m(:,:,i));
    legend('VC1','VC3','VC5','VC7','VC9','VC11','VC13')
    title(['Unit step response for input ',inputs(i)]);
    xlabel('seconds')
    ylabel('amplitude (V)')
end

% error systems plots
stfig_err = figure(6);
for i = 1:in
    subplot(2,4,i), plot(t,y_o(:,:,i)-y_m(:,:,i),'b-');
    hold on
    plot(t,y_o(:,:,i)-y_b(:,:,i),'r-');
    title(['Error in step response for input ',inputs(i),' at outputs']);
    xlabel('seconds')
    ylabel('amplitude (V)')
end

% Bode plot for the outputs
bf = logspace(-3,3,60);
[bamp, ~] = sigma_sys(sys,bf);
[bamp_bal, ~] = sigma_sys(sys_bal,bf);
[bamp_mod, ~] = sigma_sys(sys_mod,bf);
bfig = figure(7);
semilogx(bf,20*log10(bamp(1:size(bamp,1),:)),'*-');
legend('VCap1','VCap3','VCap5','Vcap7','VCap9','Vcap11','Vcap13')
grid on
set(gca,'xscale','log');
title('Bode plot of system outputs - Original system')
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');

bfig_b = figure(8);
semilogx(bf,20*log10(bamp_bal(1:size(bamp_bal,1),:)),'*-');
legend('VCap1','VCap3','VCap5','Vcap7','VCap9','Vcap11','Vcap13')
grid on
set(gca,'xscale','log');
title('Bode plot of system outputs - Balanced red.')
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');

bfig_m = figure(9);
semilogx(bf,20*log10(bamp_mod(1:size(bamp_mod,1),:)),'*-');
legend('VCap1','VCap3','VCap5','Vcap7','VCap9','Vcap11','Vcap13')
grid on
set(gca,'xscale','log');
title('Bode plot of system outputs - Modal approx.')
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');

bfig_err = figure(10);
hold on
semilogx(bf,20*log10(bamp(1:size(bamp,1),:))-20*log10(bamp_bal(1:size(bamp_bal,1),:)),'b');
semilogx(bf,20*log10(bamp(1:size(bamp,1),:))-20*log10(bamp_mod(1:size(bamp_mod,1),:)),'r');
grid on
set(gca,'xscale','log');
title('Bode plot of error systems')
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');

% error norms 
% original system
sys_hinf = normhinf(sys.A,sys.B,sys.C,sys.D);
sys_b_hinf = normhinf(sys_bal.A,sys_bal.B,sys_bal.C,sys_bal.D);
sys_m_hinf = normhinf(sys_mod.A,sys_mod.B,sys_mod.C,sys_mod.D);

sys_h2 = normh2(sys.A,sys.B,sys.C,sys.D);
sys_b_h2 = normh2(sys_bal.A,sys_bal.B,sys_bal.C,sys_bal.D);
sys_m_h2 = normh2(sys_mod.A,sys_mod.B,sys_mod.C,sys_mod.D);