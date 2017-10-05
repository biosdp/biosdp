clear;
close all;
deg = 4; %(u is polynomial of degree deg/2)
control_scaling = 1e-3; % Put it as a variable
xinit = [ 0.0664; 0; 0; 0; 0];
exp_table = {[4*3600;7*3600; 1.6632e-04];[8*3600;11*3600;8.8358e-04];[16*3600;19*3600;0.0036];[24*3600;27*3600;0.0041];[32*3600;35*3600;0.0041];[42*3600;45*3600;0.005];[52*3600;55*3600;0.0041]};
result_table = [7*3600 1.6632e-04;11*3600 8.8358e-04;19*3600 0.0036;27*3600 0.0041;35*3600 0.0041;45*3600 0.005;55*3600 0.0041];
result_scaling_factor = 481/0.005;
lastTf = 0; 
tt_exectime = 0;
infos_control = {};
epsilon = 0.005; %~0.5%err max
sqsum = sum((3.*result_table(:,2))); %% sum(x_meas[i])
%sqsum = sum((3.*result_table(:,2)).^2); %% sum(x_meas[i]^2)
opt_meansq_err = 0;
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
ulast = 0 ;
tic
for i = 1:length(exp_table)
    if i == 1
        control_scaling = 1e-3;
    elseif i == 2
        control_scaling = 1e-3;
    else
        control_scaling = 8e-3;
    end
    %% Temp optimization
    current_exp = exp_table{i};
    T1 = current_exp(1)-lastTf; % starting experiment time
    Tf = current_exp(2)-lastTf; % total time % intial point (always in temoin experiment)
    xref = current_exp(3); % Target value (minimization in cost function )
    optim_out = local_hybrid_optim( T1,Tf,xinit,deg,xref,control_scaling,ulast);
    x = optim_out.x;
    t = optim_out.t;
    uvar = optim_out.u;

    %% different dependencies
    v1 = {[t];  [t]};
    v2 = {[t;x{1}(1)];  [t;x{2}(1)]};
    v3 = {[t;x{1}(1)];  [t;x{2}(1);x{2}(5)]};
    v4 = {[t;x{1}];  [t;x{2}]};
    % different degrees : [0,2,deg] (/2)
    %% possible controls
    c1.v = v1;c1.d = 0;
    c2.v = v1;c2.d = 2;
    c3.v = v1;c3.d = deg;
    c4.v = v2;c4.d = 2;
    c5.v = v2;c5.d = deg;
    c6.v = v3;c6.d = 2;
    c7.v = v3;c7.d = deg;
    c8.v = v4;c8.d = 2;
    c9.v = v4;c9.d = deg;
%     control_test_table = {c1;c2;c3;c4;c5;c6;c7;c8;c9};  
%     control_test_table = {c1;c2;c3}  ;
     control_test_table = {c1};  
    min_err = 1e30;
    for j=1:length(control_test_table)
         %% Temp sub-degree control synth
         synth_vars = control_test_table{j}.v;
         out_csynth = control_synth(optim_out.out_solver,uvar,synth_vars,control_test_table{j}.d);
         u = out_csynth.u;
         disp('done');
         %% simu
         % Mode 1
        controller1 =@(tt,xx) double(subs(u{1},[t;x{1}],[tt;xx])); %% ode_options = odeset('Events', @simu_guard1);
        mode = 1; 
        [ tval1, xval1 ] = ode15s( @(tt,xx) simu_ode( tt, xx, controller1, mode,Tf,control_scaling), [0:0.001:T1/Tf], xinit,options);%, ode_options ); % put precision options ...

        % Mode 2
        controller2 =@(tt,xx) double(subs(u{2},[t;x{2}],[tt;[xx(1:2);xx(7);xx(4:6);xx(8:9)]])); %% ode_options = odeset('Events', @simu_guard2);
        mode = 2; 
        [ tval2, xval2 ] = ode15s( @(tt,xx) simu_ode( tt, xx, controller2,mode,Tf,control_scaling ), ...
                                  [tval1(end):0.001:1], [xval1(end,1:4) 0 0 xval1(end,3) 0 xval1(end,5)],options);%, ode_options );

        utt1 = tval1;
        utt2 = tval2;
        % u(t) in the mode 1
        for i =1:length(tval1)
             uval_t = controller1(tval1(i),xval1(i,:)');
             uval_t(uval_t>(1/control_scaling)) = 1/control_scaling;
             uval_t(uval_t<0) = 0;
             utt1(i) = control_scaling*uval_t;
        end

        % u(t) in the mode 2
        for i =1:length(tval2)
             uval_t = controller2(tval2(i),xval2(i,:)');
             uval_t(uval_t>(1/control_scaling)) = 1/control_scaling;
             uval_t(uval_t<0) = 0;
             utt2(i) = control_scaling*uval_t;
        end
                    
        %% Quality of results
        %%% Near precision error
        %err = (xval2(end,6)+4*xval2(end,8) - 3*xref).^2/((3*xref)^2)
        %%%% Residual error as in the previous paper
        err = sqrt((xval2(end,6)+4*xval2(end,8) - 3*xref).^2)/sqsum
        infos.c = control_test_table{j};
        infos.u = u;
        infos.e = err;  
        infos.tsim1 = (tval1*Tf) + lastTf;
        infos.tsim2 = (tval2*Tf) + lastTf;
        infos.xsim1 = xval1;
        infos.xsim2 = xval2;
        infos.usim1 = utt1;
        infos.usim2 = utt2;
        if err<=epsilon
            min_err = err;
            min_infos = infos;
            break;
        end

        if(err<=min_err)
            min_err = err;
            min_infos = infos;
        end
    end
    min_err
    opt_meansq_err = opt_meansq_err + min_err;
    infos_control = [infos_control; min_infos];
    ulast = utt2(end)/control_scaling % ;0
    xinit = [xval2(end,1:4)'; 0];
    lastTf = Tf+lastTf;
    lastTf/3600
end
opt_meansq_err = opt_meansq_err
tt_exectime = toc/60;

%% Saving results
save_results.xinit = xinit;
save_results.exp_table = exp_table;
save_results.result_table = result_table;
save_results.infos_control = infos_control;
save_results.sqsum = sqsum;
save('results.mat','save_results');


