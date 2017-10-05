function xdot = simu_sigmoidode( t, x, controller, mode,T,scaling )

%     disp('Inside of ODE')
    if mode ==  1
        xval = x(1:5);
    else
        xval = x(1:9);
    end 
    
%     disp('Not negative values')
%     for i=1:length(xval)
%         xloc = xval(i);
%         xloc(xloc<0) = 0;
%         xval(i) = xloc;
%     end
    uval = controller(t,xval); % control synthetiser

    uval(uval>1) = 1;
    uval(uval<0) = 0;
    
    k2 = 0; %3.78e-10;
    k4 = 4.47e-4;
    k5 = 7.27e-6;
    k6 = 4.47e-4;
    k7 = 0; %4.97e-10;
    k8 = 1.14e-5;
    
    if mode ==  1
        xdot = T*[ 3.35e-5/3600 - k2*xval(1) - (uval)*xval(1) ;
                  -k4*xval(2) - 4*k5*xval(2)*xval(3) + (uval)*xval(1);
                  k6*xval(2) - 4*k5*xval(2)*xval(3) - k7*xval(3) ;
                  k5*xval(2)*xval(3) - k8*xval(4) ;
                  1/T ];

    elseif mode == 2
        xdot = T*[  3.35e-5/3600 - k2*xval(1) - (uval)*xval(1); % Fe
                  -k4*xval(2) - 4*k5*xval(2)*xval(3) + (uval)*xval(1); % H
                  k6*xval(2) - 4*k5*xval(2)*xval(3) - k7*xval(3) ; % G without radioactive influence
                  k5*xval(2)*xval(3) - k8*xval(4) ; %Hb
                  0.0251/3600  - k2*xval(5) - (uval)*xval(5); %Fe59
                  -k4*xval(6) - 4*k5*xval(6)*xval(7) + (uval)*xval(5); %H59              
                  k6*(xval(2)+xval(6)) - 4*k5*(xval(2) + xval(6))*xval(7) - k7*xval(7) ; %G (depend of H and H59)
                  k5*xval(6)*xval(7) - k8*xval(8) ; %Hb59
                  1/T ]; %Clock    
    end
    

end
