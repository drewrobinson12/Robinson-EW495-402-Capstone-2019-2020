function [c,ceq] = nonlconZ(vars,x)
N = vars(1);
V = vars(2);
Wx = vars(3);
Wy = vars(4);
x_nm = vars(5);
y_nm = vars(6);

y = zeros(3,N);
% reformat x from single vector to 3-rows [dx;dy;theta]
    for p = 1:N
        y(:,p) = x((p-1)*3+1:(p-1)*3+3,1);
    end

c = [];

ceq = y(1:2,1:N)*([zeros(1,(N-1));eye((N-1),(N-1))]-[eye((N-1),(N-1));...
    zeros(1,(N-1))])-x(end)/(N-1)*([V*ones(1,(N-1)).*cos(y(3,1:(N-1))) + Wx*ones(1,(N-1));...
                     V*ones(1,(N-1)).*sin(y(3,1:(N-1))) + Wy*ones(1,(N-1))]);

 end

