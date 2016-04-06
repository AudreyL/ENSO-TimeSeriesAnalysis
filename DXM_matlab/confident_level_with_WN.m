function [EAC_cl_down,EAC_cl_up]=confident_level_with_WN(variance,VM,repet,N,umax,vmax,pmin,pmax,level)
    Tab_EAC_value=[];
    for i=1:repet
        disp(i)
        a=10*rand(1,1);                         
        randn('seed',floor(sum(10*clock)/a));                   % changes the seed
        Simul_WN=randn(1,N)*sqrt(variance);                    % simulates wihte noise
        if (~isnan(VM))
        Simul_WN(VM)=NaN;
        end
        EAC_profile_WN=EAC_value(Simul_WN,N,umax,vmax,pmin,pmax); % calculate EAC value for white noise
        Tab_EAC_value=[Tab_EAC_value ; EAC_profile_WN];     % array of repet EAC values for white noise
    end
    EAC_cl_up=quantile(Tab_EAC_value,level);                % calculates the standard deviation of EAC values for each period
    EAC_cl_down=quantile(Tab_EAC_value,1-level);            % calculates the mean of EAC values for each period
    
    
    


    
