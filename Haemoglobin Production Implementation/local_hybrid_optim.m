%%%%%%%%%%%%%%%%%%
%
% Function used to descibe the local hybrid dynamical system corresponding to the  Haemoglobin
% production model on a given time interval as a 2 mode controlled hybrid system. 
% Inputs are: 
% T1 the initial time (normaly set to 0 previously)
% Tf the final time for the local optimization
% xint the intial condition
% relaxdeg the relaxation degree for the future optimization
% xrefTf the target data point at Tf
% scaling the scaling locally applied on the control
% ulast the last value of the control
%%%%%%%%%%%%%%%%%%
function [ out ] = local_hybrid_optim( T1,Tf,xinit,relaxdeg,xrefTf,scaling,ulast )

    T = Tf;         % maximum time horizon
    d = relaxdeg;   % degree of relaxation
    nmodes = 2;     % number of hybrid modes


    % ========================= Define variables =============================
    t = msspoly( 't', 1 );
    xa = msspoly( 'x', 8 );
    ua = msspoly( 'u', 1 );

    x = cell( nmodes, 1 );
    u = cell( nmodes, 1 );
    f = cell( nmodes, 1 );
    g = cell( nmodes, 1 );
    x0 = cell( nmodes, 1 );
    hX = cell( nmodes, 1 );
    hU = cell( nmodes, 1 );
    hXT = cell( nmodes, 1 );
    sX = cell( nmodes, nmodes );
    R = cell( nmodes, nmodes );
    h = cell( nmodes, 1 );
    H = cell( nmodes, 1 );
        
    x0{1} = xinit; % INTIAL CONDITION

    xref = xrefTf; % Target point

    trans_time = T1/Tf;

    % ===== Fixed parameters (taken from Table 1 in paper, column pmean)
    % We will search to synthesis of k3(t)

    %PARAMETERS NOT SCALED (scaled in the equations)
    k2 = 0; %3.78e-10;
    k4 = 4.47e-4;
    k5 = 7.27e-6;
    k6 = 4.47e-4;
    k7 = 0; %4.97e-10;
    k8 = 1.14e-5;
    control_scaling = scaling;
    %Feext input: 3.35e-5/3600 % inline
    %Fe59ext input: 0.0251/3600 

    % =============================== Dynamics ===============================
    % -------  Mode 1 --------
    x{1} = xa([1:4,8]); % 5 variables (4 temoin + clock)
    u{1} = ua;
    f{1} = T * [  3.35e-5/3600 - k2*x{1}(1);
                  -k4*x{1}(2) - 4*k5*x{1}(2)*x{1}(3) ;
                  k6*x{1}(2) - 4*k5*x{1}(2)*x{1}(3) - k7*x{1}(3) ;
                  k5*x{1}(2)*x{1}(3) - k8*x{1}(4) ;
                  1/T ];
    g{1} = T * [-control_scaling*x{1}(1);control_scaling*x{1}(1);0;0;0];
    % -------  Mode 2 --------
    x{2} =  xa; % 8 variables (4 temoin + 3 radio + clock)
    u{2} = ua;              
    f{2} = T * [  3.35e-5/3600 - k2*x{2}(1); % Fe
                  -k4*x{2}(2) - 4*k5*x{2}(2)*x{2}(3) ; % H
                  k6*(x{2}(2)+x{2}(6)) - 4*k5*(x{2}(2)+x{2}(6))*x{2}(3) - k7*x{2}(3) ; %G (depend of H and H59)
                  k5*x{2}(2)*x{2}(3) - k8*x{2}(4) ; %Hb
                  0.0251/3600  - k2*x{2}(5); %Fe59
                  -k4*x{2}(6) - 4*k5*x{2}(6)*x{2}(3) ; %H59
                  k5*x{2}(6)*x{2}(3) - k8*x{2}(7) ; %Hb59
                  1/T ]; %Clock          
    g{2} = T * [-control_scaling*x{2}(1);control_scaling*x{2}(1);0;0;-control_scaling*x{2}(5);control_scaling*x{2}(5);0;0];


    % ======================== Domains and transitions =======================
    % ----------- Mode 1 -----------
    y = x{1};
    hX{1} = [ y(1) * (10-y(1));           % Fe_intern = [0,10] 
              y(2) * (1-y(2));            % Heme = [0,1] 
              y(3) * (1-y(3));           % Globin = [0,1] 
              y(4) * (1-y(4));          % HemoGlobin = [0,1] 
              y(5) * (trans_time-y(5))];         % clock = [0,trans_time] (4h on 8)


    hU{1} = [ua * ((1/control_scaling)-ua)];                 % U_1 = [0,1]

    % Transition 1->2
    sX{1,2} = [ y(1) * (1-y(1));    % Fe_intern = [0,1] 
                y(2) * (1-y(2));     % Heme = [0,1] 
                y(3) * (1-y(3));     % Globin = [0,1] 
                y(4) * (1-y(4));     % HemoGlobin = [0,1] 
                 -(trans_time - y(5))^2]; % Clock = trans_time

    R{1,2} = [y(1:4);0;0;0;y(5)];       % R_12 = Identity for existing variables and  0 for radioactive variable initial concentration


    % ----------- Mode 2 -----------
    y = x{2};
    hX{2} = [ y(1) * (10-y(1));            % Fe_intern = [0,1] (mol/m^3)
              y(2) * (1-y(2));            % Heme = [0,1] 
              y(3) * (1-y(3));            % Globin = [0,1] 
              y(4) * (1-y(4));            % HemoGlobin = [0,1] 
              y(5) * (1-y(5));            % Fe_59_intern = [0,1]
              y(6) * (1-y(6));            % Heme_59 = [0,1]
              y(7) * (1-y(7));            % HemoGlobin_59 = [0,1]
              y(8) - trans_time ];            % clock = [trans_time,Tf] (~ hours scaled on [0,1])

    hU{2} = [ua * ((1/control_scaling)-ua)]; % U \in [0 1]

    % % ===================== Cost functions and target set ====================
    if(ulast ~= 0)
        h{1} = (scaling*(ulast-ua))^2; % cost to the last iteration control
        %%%%%%%%
        y = x{2};
        h{2} = (scaling*(ulast-ua))^2 + (0.01*ua)^2;
        H{2} = (y(6)+4*y(7) - 3*xref)^2; % Cost associated to the target
        hXT{2} = hX{2};
    else
        h{1} = 0;
        %%%%%%%%
        y = x{2};
        h{2} = (0.01*ua)^2 ; % Cost on the amplitude of the control
        H{2} = (y(6)+4*y(7) - 3*xref)^2;
        hXT{2} = hX{2}; 
    end


    % ============================== Options =================================
    options.freeFinalTime = 0;
    options.withInputs = 1;
    options.svd_eps = 1e12;

    % Call a function to build the spotless structure of the implmentation, and the mosek solver 
    [out_prog] = HybridOCPDualSolver_redo(t,x,u,f,g,hX,hU,sX,R,x0,hXT,h,H,d,options);

    pval = out_prog.pval;
    disp(['LMI deg = ' int2str(d) ' lower bound = ' num2str(pval)]);

    out.out_solver = out_prog;
    out.x = x;
    out.t = t;
    out.u = u;
end

