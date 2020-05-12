function [theta_c] = unit2comp(theta_u)
   
if isnan(theta_u) == 0
    
    theta_u = rad2deg(theta_u);
    
    for i = 1:length(theta_u)
        if isnan(theta_u)
            theta_c(i) = 0;
        elseif theta_u(i) <= 0
            theta_c(i) = 90-theta_u(i);
        elseif theta_u(i) > 0 & theta_u(i) < 90
            theta_c(i) = 90-theta_u(i);
        elseif theta_u(i) == 90
            theta_c(i) = 0;
        elseif theta_u(i) > 90 & theta_u(i) < 180
            theta_c(i) = 450-theta_u(i);
        elseif theta_u(i) == 180
            theta_c(i) = 270;
        elseif theta_u(i) > 180 & theta_u(i) < 270
            theta_c(i) = 450-theta_u(i);
        elseif theta_u(i) == 270
            theta_c(i) = 180;
        elseif theta_u(i) > 270 & theta_u(i) < 360
            theta_c(i) = 450-theta_u(i);
        elseif theta_u(i) == 360
            theta_c(i) = 90;
        end
    end
else
theta_c = 0;
end 

