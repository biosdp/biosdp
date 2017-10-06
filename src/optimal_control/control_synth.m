%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function for Truncated control synthesis
%  This is based on the previous implementation of
%  https://github.com/pczhao/hybridOCP
%
% The input are :
% prog_sol : all the information to perform the control synthesis
% u : a precomputed form for the control
% vars : the variable of the desired control
% d : the degree of the desired control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ out ] = control_synth( prog_sol, u , vars, d )

u_infos =  prog_sol.u_infos;
nmodes = u_infos.nmodes;
uout = u_infos.uout;
u_real_basis = u_infos.u_real_basis;
mu_idx = u_infos.mu_idx;
svd_eps = u_infos.svd_eps;
dual_basis = u_infos.dual_basis;
y = u_infos.y;
sol = prog_sol.sol;

for i = 1 : nmodes
     urb = monomials( [ vars{i} ], 0 : d/2 );    
    u_real_basis{ i } = urb;
    % moments of mu
    mu_basis = dual_basis{ mu_idx(i) };
    y_mu = sol.dualEval( y{ mu_idx(i) } );
    
    % moment matrix of mu
    mypoly = mss_s2v( urb * urb' );
    coeff_mu = DecompMatch( mypoly, mu_basis );
    moment_mat = double( mss_v2s( coeff_mu * y_mu ) );
    
    [U,S,V] = svd(moment_mat);
    iS1 = S;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Default inverse Approx method %%%% 
%     startExp = 1e-10;   
%     while norm(pinv(iS1)) / norm(S) > svd_eps
%         iS1(iS1 < startExp) = 0;
%         startExp = startExp * 10;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    fprintf('norm of moment matrix %0.2f\n', norm(S));
    fprintf('norm of inverse moment matrix %0.2f\n', norm(pinv(iS1)));

    iMyt = V*pinv(iS1)*U';
   
    % yp and yn
    for j = 1 : length( u{ i } )
        fprintf('Processing mode %1.0f, input #%1.0f ...\n', i, j );
        mypoly = u{ i }( j ) * urb;
        coeff = DecompMatch( mypoly, mu_basis );
        y_u = sol.dualEval( coeff * y_mu );
        
        u_coeff = iMyt * y_u;
        uout{ i, j } = u_coeff' * urb;
        
    end
end
out.u = uout;
end

