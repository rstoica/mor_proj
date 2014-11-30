function [ sys,res ] = solve_sys( sys,t,f )
% solv_sys Solves different problems and reduces a given system 
%   
%   INPUT: sys = system to solve
%   OUTPUT: sys = same system unmodified
%           res = results struct containing the step response and the bode
%           plot of the system

    % plot the step response of the system
    [res.y_step,~,res.x_step] = step_sys(sys,t);
    in = size(res.y_step,3);
    inputs = {['1=I1'],['2=I3'],['3=I5'],['4=I7'],['5=I9'],['6=I11'],['7=I13'],['8=VDD']};
    res.fig = zeros(3,1);
    res.fig(1) = figure();
    hold on
    for i = 1:in
        subplot(2,4,i), plot(t,res.y_step(:,:,i));
        legend('VC1','VC3','VC5','VC7','VC9','VC11','VC13')
        title(['Unit step response: input ',inputs(i),' sys: ',sys.name]);
        xlabel('seconds')
        ylabel('amplitude (V)')
    end
    res.t_step = t;

    % plot the Bode plot of the system
    [res.y_bode, ~] = sigma_sys(sys,f);
    res.f_bode = f;
    res.fig(2) = figure();
    semilogx(f,20*log10(res.y_bode(1:size(res.y_bode,1),:)),'*-');
    legend('VCap1','VCap3','VCap5','Vcap7','VCap9','Vcap11','Vcap13')
    grid on
    set(gca,'xscale','log');
    title(['Bode plot of ',sys.name,' system outputs']);
    xlabel('Frequency (rad/s)');
    ylabel('Amplitude (dB)');
    res.hinf = normhinf(sys.A,sys.B,sys.C,sys.D);
    fprintf('normhinf sys: %s = %.4f\n',sys.name,res.hinf);
    res.h2 = normh2(sys.A,sys.B,sys.C,sys.D);
    fprintf('normh2 sys: %s = %.4f\n',sys.name,res.h2);
end

