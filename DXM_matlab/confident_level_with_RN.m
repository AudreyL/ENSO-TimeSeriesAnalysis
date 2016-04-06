function [EAC_cl_down,EAC_cl_up]=confident_level_with_RN(foac,VM,repet,N,umax,vmax,pmin,pmax,level)
    Tab_EAC_value=[];
    for i=1:repet
        disp(i)
        a=10*rand(1,1);                         
        randn('seed',floor(sum(10*clock)/a));   % changes the seed
        G(1,1) = randn(1,1);
        for k=1:(N-1)
            G(1,k+1) = (foac*G(1,k)) + ((1-foac)*randn(1,1)); 
        end
        G = G - repmat(mean(G),1,1);
        C = G/std(G);
        if (~isnan(VM))
        C(VM)=NaN;
        end
        EAC_profile_WN=EAC_value(C,N,umax,vmax,pmin,pmax); % calculate EAC value for white noise
        Tab_EAC_value=[Tab_EAC_value ; EAC_profile_WN]; % array of repet EAC values for white noise
    end
    EAC_cl_up=quantile(Tab_EAC_value,level);              % calculates the standard deviation of EAC values for each period
    EAC_cl_down=quantile(Tab_EAC_value,1-level);  