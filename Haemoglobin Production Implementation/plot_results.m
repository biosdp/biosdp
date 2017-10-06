%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function used to generate the figures, and compute three additionnal
% biological meaningful fits (step-function, piece-wise polynomial, sigmoid)
% of the generated control
%
% The input is a result file produced by the run_iterative_optimal_control
% script
%
% This implementation is dedicated solely to the Haemoglobin production model.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ done ] = plot_results(save_file)

close all;

load(save_file) 
xinit = save_results.xinit;
exp_table = save_results.exp_table;
result_table = save_results.result_table;
infos_control = save_results.infos_control;
sqsum = save_results.sqsum;

options = odeset('RelTol',1e-8,'AbsTol',1e-10);


%%  synth control plot
result_scaling_factor = 481/0.005;
tsim = [];    
utt  = [];
Fe = []; % Iron (not radioactive)
H = [];  % Heme (not radioactive)
G = [];  % Globin (not radioactive)
Hb = []; % Hemoglobin (not radioactive)
Fe59 = []; % Iron (not radioactive)
H59 = [];  % Heme (not radioactive)
Gtot = [];  % Globin (not radioactive)
Hb59 = []; % Hemoglobin (not radioactive)

for i = 1:length(infos_control)
    tsim = [tsim;infos_control{i}.tsim1;infos_control{i}.tsim2];
    utt = [utt;infos_control{i}.usim1;infos_control{i}.usim2];
    Fe = [Fe;infos_control{i}.xsim1(:,1);infos_control{i}.xsim2(:,1)]; % Iron (not radioactive)
    H = [H;infos_control{i}.xsim1(:,2);infos_control{i}.xsim2(:,2)];  % Heme (not radioactive)
    G = [G;infos_control{i}.xsim1(:,3);infos_control{i}.xsim2(:,3)];  % Globin (not radioactive)
    Hb = [Hb;infos_control{i}.xsim1(:,4);infos_control{i}.xsim2(:,4);]; % Hemoglobin (not radioactive)
%     Fe59 = [Fe59;infos_control{i}.xsim2(:,5)]; % Iron (not radioactive)
%     H59 = [H59;infos_control{i}.xsim2(:,6)];  % Heme (not radioactive)
%     Gtot = [Gtot;infos_control{i}.xsim2(:,7)];  % Globin (not radioactive)
%     Hb59 = [Hb59;infos_control{i}.xsim2(:,8)]; % Hemoglobin (not radioactive)

end

tim_plot = tsim/3600;

figure(1);
subplot(2,2,1)       
plot(tim_plot,Fe,'b');
title('Fe(t)');
subplot(2,2,2)      
plot(tim_plot,H,'b');
title('H(t)');
subplot(2,2,3)   
plot(tim_plot,G,'b');
title('G(t)');
subplot(2,2,4) 
plot(tim_plot,Hb,'b');
title('Hb(t)');

figure(2);
plot(tim_plot,utt);
title('Generated Optimal Control (t) on [0 1]');

figure(3);
merr = 0;
for i = 1:length(infos_control)
    xref_plot =  exp_table{i}(3);
    subH59 = infos_control{i}.xsim2(:,6);
    subHb59 = infos_control{i}.xsim2(:,8);
    plot(infos_control{i}.tsim2/3600,sqrt((subH59+4*subHb59 - 3*xref_plot).^2)/sqsum,'b');
    lerr = sqrt((subH59(end)+4*subHb59(end) - 3*xref_plot)^2)/sqsum
    title('Cost Function');
    hold on;
    merr =  merr + lerr;
end

figure(4);
for i = 1:length(infos_control)
subplot(2,2,1)  
hold on;
subFe59 = infos_control{i}.xsim2(:,5);
plot(infos_control{i}.tsim2/3600,subFe59,'b');
title('Fe59(t)');   
subplot(2,2,2)   
hold on;
subGtot = infos_control{i}.xsim2(:,7);
plot(infos_control{i}.tsim2/3600,subGtot,'b');
title('Gtot(t)');
subplot(2,2,[3 4]) 
hold on;
subH59 = infos_control{i}.xsim2(:,6);
subHb59 = infos_control{i}.xsim2(:,8);
% plot(infos_control{i}.tsim2/3600,(subH59+4*subHb59)*result_scaling_factor);
plot(infos_control{i}.tsim2/3600,(subH59+4*subHb59)*result_scaling_factor/3,'b');
title('Observable variable = H59(t)+4*Hb59(t)');
hold on;
end

%% some not plotted polynomial fit (same degree than the polynomial proposed in the previous paper describing the haemoglobin model )
%p = polyfit(tsim,utt,4); 
%f1 = polyval(p,tsim);
%figure(2);
%hold on;
%plot(tim_plot,f1,'--');

%% Piece-wise polynomial control
mean_infos = {};
lastTf = 0; 
xinit = [ 0.0664; 0; 0; 0; 0];
control_scaling = 1e-3;
mc_meansq = 0;
%pmean = polyfit(tsim,utt,4);
% Fit of each sub section using polyfit matlab function
p1 = polyfit(tsim(1:2004),utt(1:2004),2);
f11 = polyval(p1,tsim(1:2004));
p2 = polyfit(tsim(2005:7014),utt(2005:7014),3);
f12 = polyval(p2,tsim(2005:7014));

figure(2);
hold on;
plot([tsim(1:2004)/3600;tsim(2005:7014)/3600],[f11;f12],'r--');


for i = 1:length(exp_table)
        current_exp = exp_table{i};
        T1 = current_exp(1)-lastTf; % starting experiment time
        Tf = current_exp(2)-lastTf; % total time % intial point (always in temoin experiment)
        xref = current_exp(3); % Target value (minimization in cost function )            
        
        if i <= 2
            control_scaling = 1e-3;
            % Mode 1  
            controller1 =@(tt,xx) polyval(p1,(tt.*Tf)+lastTf);
            mode = 1;
            [ tval1, xval1 ] = ode15s( @(tt,xx) simu_scaledode( tt, xx, controller1, mode,Tf,control_scaling), [0:0.001:T1/Tf], xinit,options);%, ode_options ); % put precision options ...
            % Mode 2
            controller2 =@(tt,xx) polyval(p1,(tt*Tf)+lastTf);
            mode = 2; 
            [ tval2, xval2 ] = ode15s( @(tt,xx) simu_scaledode( tt, xx, controller2,mode,Tf,control_scaling ), ...
                                  [tval1(end):0.001:1], [xval1(end,1:4) 0 0 xval1(end,3) 0 xval1(end,5)],options);%, ode_options );                      
        else
            control_scaling = 8e-3;
            % Mode 1  
            controller1 =@(tt,xx) polyval(p2,(tt.*Tf)+lastTf);
            mode = 1; 
            [ tval1, xval1 ] = ode15s( @(tt,xx) simu_scaledode( tt, xx, controller1, mode,Tf,control_scaling), [0:0.001:T1/Tf], xinit,options);%, ode_options ); % put precision options ...
            % Mode 2
            controller2 =@(tt,xx) polyval(p2,(tt*Tf)+lastTf);
            mode = 2; 
            [ tval2, xval2 ] = ode15s( @(tt,xx) simu_scaledode( tt, xx, controller2,mode,Tf,control_scaling ), ...
                                  [tval1(end):0.001:1], [xval1(end,1:4) 0 0 xval1(end,3) 0 xval1(end,5)],options);%, ode_options );
        end
                        
%
        %Quality of results
        err = (xval2(end,6)+4*xval2(end,8) - 3*xref).^2/sqsum;
        infos.e = err;  
        infos.tsim1 = (tval1*Tf) + lastTf;
        infos.tsim2 = (tval2*Tf) + lastTf;
        infos.xsim1 = xval1;
        infos.xsim2 = xval2;
        mean_infos = [mean_infos; infos];
        
        xinit = [xval2(end,1:4)'; 0];
        lastTf = Tf+lastTf;
        mc_meansq = mc_meansq + err;
end


mean_tsim = [];   
mean_Fe = []; % Iron (not radioactive)
mean_H = [];  % Heme (not radioactive)
mean_G = [];  % Globin (not radioactive)
mean_Hb = []; % Hemoglobin (not radioactive)

for i = 1:length(mean_infos)
    mean_tsim = [mean_tsim;mean_infos{i}.tsim1;mean_infos{i}.tsim2];
    mean_Fe = [mean_Fe;mean_infos{i}.xsim1(:,1);mean_infos{i}.xsim2(:,1)]; % Iron (not radioactive)
    mean_H = [mean_H;mean_infos{i}.xsim1(:,2);mean_infos{i}.xsim2(:,2)];  % Heme (not radioactive)
    mean_G = [mean_G;mean_infos{i}.xsim1(:,3);mean_infos{i}.xsim2(:,3)];  % Globin (not radioactive)
    mean_Hb = [mean_Hb;mean_infos{i}.xsim1(:,4);mean_infos{i}.xsim2(:,4);]; % Hemoglobin (not radioactive)

end

mean_tsim_plot = mean_tsim/3600;
figure(1);
hold on;
subplot(2,2,1) 
hold on;
plot(mean_tsim_plot,mean_Fe,'r:');
title('Fe(t)');
subplot(2,2,2)  
hold on;
plot(mean_tsim_plot,mean_H,'r:');
title('H(t)');
subplot(2,2,3)   
hold on;
plot(mean_tsim_plot,mean_G,'r:');
title('G(t)');
subplot(2,2,4)
hold on;
plot(mean_tsim_plot,mean_Hb,'r:');
title('Hb(t)');


figure(3);
hold on;
merr2 = 0;
for i = 1:length(mean_infos)
    xref_plot =  exp_table{i}(3);
    subH59 = mean_infos{i}.xsim2(:,6);
    subHb59 = mean_infos{i}.xsim2(:,8);
    plot(mean_infos{i}.tsim2/3600,sqrt((subH59+4*subHb59 - 3*xref_plot).^2)/sqsum,'r:' );
    lerr2 = sqrt((subH59(end)+4*subHb59(end) - 3*xref_plot)^2)/sqsum
    title('Cost Function');
    hold on;
    merr2 =  merr2 + lerr2;
end

figure(4);
hold on;
for i = 1:length(mean_infos)
subplot(2,2,1)
hold on;
subFe59 = mean_infos{i}.xsim2(:,5);
plot(mean_infos{i}.tsim2/3600,subFe59,'r:');
title('Fe59(t)');   
subplot(2,2,2)
hold on;
subGtot = mean_infos{i}.xsim2(:,7);
plot(mean_infos{i}.tsim2/3600,subGtot,'r:');
title('Gtot(t)');
subplot(2,2,[3 4])
hold on;
subH59 = mean_infos{i}.xsim2(:,6);
subHb59 = mean_infos{i}.xsim2(:,8);
plot(mean_infos{i}.tsim2/3600,(subH59+4*subHb59)*result_scaling_factor/3,'r:' );
title('Observable variable = H59(t)+4*Hb59(t)');
end

%% Sigmoid Control
mean_infos = {};
lastTf = 0; 
xinit = [ 0.0664; 0; 0; 0; 0];
control_scaling = 1e-3;
mc_meansq = 0;

figure(2);
hold on;
psig = @(tt) 1.8e-04.*((tt).^7./((14*3600)^7+(tt).^7)) +0.1e-5;
fsig = psig(tsim);
plot(tim_plot,fsig,'k');
for i = 1:length(exp_table)
         current_exp = exp_table{i};
        T1 = current_exp(1)-lastTf; % starting experiment time
        Tf = current_exp(2)-lastTf; % total time % intial point (always in temoin experiment)
        xref = current_exp(3); % Target value (minimization in cost function )  
      
%%%%% SIG smoothing
      
        
        controller1 =@(tt,xx) psig((tt.*Tf)+lastTf);
        mode = 1; 
        [ tval1, xval1 ] = ode15s( @(tt,xx) simu_scaledode( tt, xx, controller1, mode,Tf,control_scaling), [0:0.001:T1/Tf], xinit,options);%, ode_options ); % put precision options ...
        % Mode 2
        controller2 =@(tt,xx) psig((tt.*Tf)+lastTf);
        mode = 2; 
        [ tval2, xval2 ] = ode15s( @(tt,xx) simu_scaledode( tt, xx, controller2,mode,Tf,control_scaling ), ...
        [tval1(end):0.001:1], [xval1(end,1:4) 0 0 xval1(end,3) 0 xval1(end,5)],options);%, ode_options );
        
        %Quality of results
        err = (xval2(end,6)+4*xval2(end,8) - 3*xref).^2/sqsum;
        infos.e = err;  
        infos.tsim1 = (tval1*Tf) + lastTf;
        infos.tsim2 = (tval2*Tf) + lastTf;
        infos.xsim1 = xval1;
        infos.xsim2 = xval2;
        mean_infos = [mean_infos; infos];
        
        xinit = [xval2(end,1:4)'; 0];
        lastTf = Tf+lastTf;
        mc_meansq = mc_meansq + err;
end

mean_tsim = [];   
mean_Fe = []; % Iron (not radioactive)
mean_H = [];  % Heme (not radioactive)
mean_G = [];  % Globin (not radioactive)
mean_Hb = []; % Hemoglobin (not radioactive)

for i = 1:length(mean_infos)
    mean_tsim = [mean_tsim;mean_infos{i}.tsim1;mean_infos{i}.tsim2];
    mean_Fe = [mean_Fe;mean_infos{i}.xsim1(:,1);mean_infos{i}.xsim2(:,1)]; % Iron (not radioactive)
    mean_H = [mean_H;mean_infos{i}.xsim1(:,2);mean_infos{i}.xsim2(:,2)];  % Heme (not radioactive)
    mean_G = [mean_G;mean_infos{i}.xsim1(:,3);mean_infos{i}.xsim2(:,3)];  % Globin (not radioactive)
    mean_Hb = [mean_Hb;mean_infos{i}.xsim1(:,4);mean_infos{i}.xsim2(:,4);]; % Hemoglobin (not radioactive)

end

mean_tsim_plot = mean_tsim/3600;
figure(1);
hold on;
subplot(2,2,1) 
hold on;
plot(mean_tsim_plot,mean_Fe,'k--');
title('Fe(t)');
subplot(2,2,2)  
hold on;
plot(mean_tsim_plot,mean_H,'k--');
title('H(t)');
subplot(2,2,3)   
hold on;
plot(mean_tsim_plot,mean_G,'k--');
title('G(t)');
subplot(2,2,4)
hold on;
plot(mean_tsim_plot,mean_Hb,'k--');
title('Hb(t)');


figure(3);
hold on;
merrsig = 0;
for i = 1:length(mean_infos)
    xref_plot =  exp_table{i}(3);
    subH59 = mean_infos{i}.xsim2(:,6);
    subHb59 = mean_infos{i}.xsim2(:,8);
    plot(mean_infos{i}.tsim2/3600,sqrt((subH59+4*subHb59 - 3*xref_plot).^2)/sqsum,'k--' );
    lerrsig = sqrt((subH59(end)+4*subHb59(end) - 3*xref_plot)^2)/sqsum
    title('Cost Function');
    hold on;
    merrsig =  merrsig + lerrsig;
end

figure(4);
hold on;
for i = 1:length(mean_infos)
subplot(2,2,1)
hold on;
subFe59 = mean_infos{i}.xsim2(:,5);
plot(mean_infos{i}.tsim2/3600,subFe59,'k--');
title('Fe59(t)');   
subplot(2,2,2)
hold on;
subGtot = mean_infos{i}.xsim2(:,7);
plot(mean_infos{i}.tsim2/3600,subGtot,'k--');
title('Gtot(t)');
subplot(2,2,[3 4])
hold on;
subH59 = mean_infos{i}.xsim2(:,6);
subHb59 = mean_infos{i}.xsim2(:,8);
plot(mean_infos{i}.tsim2/3600,(subH59+4*subHb59)*result_scaling_factor/3,'k--' );
title('Observable variable = H59(t)+4*Hb59(t)');
end

%% Step function control
mean_infos = {};
lastTf = 0; 
xinit = [ 0.0664; 0; 0; 0; 0];
mc_meansq = 0;
tswitch = [];
uplot_switch = [];
for i = 1:length(exp_table)
        
        %if i <= 2
        if i==1    
            control_scaling = 1e-3;
            controller1 =@(tt,xx) sum(utt(1:2004))/2004; 
            controller2 =@(tt,xx) sum(utt(1:2004))/2004;

        elseif i==2
            control_scaling = 1e-3;
            controller1 =@(tt,xx) sum(utt(1:2004))/2004; 
            controller2 =@(tt,xx) sum(utt(1:2004))/2004;

        else
            control_scaling = 8e-3;
            controller1 =@(tt,xx)  sum(utt(2005:7014))/5010;
            controller2 =@(tt,xx)  sum(utt(2005:7014))/5010;

        end
        current_exp = exp_table{i};
        T1 = current_exp(1)-lastTf; % starting experiment time
        Tf = current_exp(2)-lastTf; % total time % intial point (always in temoin experiment)
        xref = current_exp(3); % Target value (minimization in cost function )
       
        % Mode 1       
        mode = 1; 
        [ tval1, xval1 ] = ode15s( @(tt,xx) simu_sigmoidode( tt, xx, controller1, mode,Tf,control_scaling), [0:0.001:T1/Tf], xinit,options);%, ode_options ); % put precision options ...
       
        % Mode 2
        mode = 2; 
        [ tval2, xval2 ] = ode15s( @(tt,xx) simu_sigmoidode( tt, xx, controller2,mode,Tf,control_scaling ), ...
                                  [tval1(end):0.001:1], [xval1(end,1:4) 0 0 xval1(end,3) 0 xval1(end,5)],options);%, ode_options );
                  


        %% Quality of results
        err = (xval2(end,6)+4*xval2(end,8) - 3*xref).^2/sqsum;
        infos.e = err;  
        infos.tsim1 = (tval1*Tf) + lastTf;
        infos.tsim2 = (tval2*Tf) + lastTf;
        infos.xsim1 = xval1;
        infos.xsim2 = xval2;
        mean_infos = [mean_infos; infos];
        %% plot sigmoid control

        
        if i<=2 

            tswitch =  [tswitch;[infos.tsim1;infos.tsim2]/3600];
            uplot_switch = [uplot_switch;repmat(sum(utt(1:2004))/2004,size([infos.tsim1;infos.tsim2]))];
        else
            tswitch =  [tswitch;[infos.tsim1;infos.tsim2]/3600];
            uplot_switch = [uplot_switch;repmat(sum(utt(2005:7014))/5010,size([infos.tsim1;infos.tsim2]))];
        end

        
        xinit = [xval2(end,1:4)'; 0];
        lastTf = Tf+lastTf;
        mc_meansq = mc_meansq + err;
end


figure(2);
hold on;
plot(tswitch,uplot_switch,'g--');


mean_tsim = [];   
mean_Fe = []; % Iron (not radioactive)
mean_H = [];  % Heme (not radioactive)
mean_G = [];  % Globin (not radioactive)
mean_Hb = []; % Hemoglobin (not radioactive)

for i = 1:length(mean_infos)
    mean_tsim = [mean_tsim;mean_infos{i}.tsim1;mean_infos{i}.tsim2];
    mean_Fe = [mean_Fe;mean_infos{i}.xsim1(:,1);mean_infos{i}.xsim2(:,1)]; % Iron (not radioactive)
    mean_H = [mean_H;mean_infos{i}.xsim1(:,2);mean_infos{i}.xsim2(:,2)];  % Heme (not radioactive)
    mean_G = [mean_G;mean_infos{i}.xsim1(:,3);mean_infos{i}.xsim2(:,3)];  % Globin (not radioactive)
    mean_Hb = [mean_Hb;mean_infos{i}.xsim1(:,4);mean_infos{i}.xsim2(:,4);]; % Hemoglobin (not radioactive)

end

mean_tsim_plot  = mean_tsim/3600 ;

figure(1);
hold on;
subplot(2,2,1) 
hold on;
plot(mean_tsim_plot,mean_Fe,'g--');
title('Fe(t)');
subplot(2,2,2)  
hold on;
plot(mean_tsim_plot,mean_H,'g--');
title('H(t)');
subplot(2,2,3)   
hold on;
plot(mean_tsim_plot,mean_G,'g--');
title('G(t)');
subplot(2,2,4)
hold on;
plot(mean_tsim_plot,mean_Hb,'g--');
title('Hb(t)');


figure(3);
hold on;
merr3 = 0;
for i = 1:length(mean_infos)
    xref_plot =  exp_table{i}(3);
    subH59 = mean_infos{i}.xsim2(:,6);
    subHb59 = mean_infos{i}.xsim2(:,8);
    plot(mean_infos{i}.tsim2/3600,sqrt((subH59+4*subHb59 - 3*xref_plot).^2)/sqsum,'g--' );
    lerr3 = sqrt((subH59(end)+4*subHb59(end) - 3*xref_plot)^2)/sqsum
    title('Cost Function');
    hold on;
    merr3 =  merr3 + lerr3;
end

figure(4);
hold on;
for i = 1:length(mean_infos)
subplot(2,2,1)
hold on;
subFe59 = mean_infos{i}.xsim2(:,5);
plot(mean_infos{i}.tsim2/3600,subFe59,'g--');
title('Fe59(t)');   
subplot(2,2,2)
hold on;
subGtot = mean_infos{i}.xsim2(:,7);
plot(mean_infos{i}.tsim2/3600,subGtot,'g--');
title('Gtot(t)');
subplot(2,2,[3 4])
hold on;
subH59 = mean_infos{i}.xsim2(:,6);
subHb59 = mean_infos{i}.xsim2(:,8);
plot(mean_infos{i}.tsim2/3600,(subH59+4*subHb59)*result_scaling_factor/3,'g--' );
title('Observable variable = H59(t)+4*Hb59(t)');
end

figure(4);
hold on;
subplot(2,2,[3 4])
hold on;
plot(result_table(:,1)/3600,result_table(:,2)*result_scaling_factor,'+');

%% Legends
figure(2);
hold on;
set(gca,'FontSize',15);
set(gca,'XTick',0:5:60);
set(gca,'YTick',0:2.5e-4/10:2.5e-4);
title('Generated Optimal control & associated fitted functions for u(t)','Interpreter','tex', 'FontSize', 15);
legend('\fontsize{15} Generated Optimal control u_{gen}(t)','\fontsize{15}piece-wise polynomial fit u_{poly}(t)','\fontsize{15}Hill function fit u_{hill}(t)','\fontsize{15}piece-wise constant fit u_{pc}(t)','Location','southeast');
ylabel('k_3(t) [hours^{-1}]','Interpreter','tex', 'FontSize', 15);
xlabel('Time [hours]','Interpreter','tex', 'FontSize', 15);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(gcf,'../Result_paper/opt_ctrl_test','-dpdf','-r0');



end

