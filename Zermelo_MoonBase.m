%Constants
N = 20; %step size
V = 2; %assumed boat velocity (kts)
Wx = 1.94*u_curr; %current speed x, m/s in --> kts out
Wy = 1.94*v_curr; %current speed y, m/s in --> kts out


% x = [dx; dy; theta] where dx - distance covered in x (east/west, longitude)
%                           dy - distance covered in y (north/south, latitude)
%                           theta - boat heading
% x is a 3*N x 1 + 1 length vector with each time step comprised of 3 rows

% **Think units. Distance covered is in meters or conversion, whereas lat
% and long are not. Grab math from trig solution script. 

%Variable Declarations


y0 = lat(end); %latitude in NM, 1 degree = 60 NM
yf = lat(1); %latitude in NM, 1 degree = 60 NM

x0 = long(1);
xf = long(end);
%long1 = 0.0034987*lat(end)^2 - 1.3302*lat(end) + 96.413; %longitude in NM of a degree at given latitude
%long2 = 0.0034987*lat(1)^2 - 1.3302*lat(1) + 96.413; %longitude in NM of a degree at given latitude

y_nm = (yf-y0)*60; %dist in NM from y0 to yf
x_low = 0.0034987*y0^2 - 1.3302*y0 + 96.413; %dist in NM at y0
x_hi = 0.0034987*yf^2 - 1.3302*yf + 96.413; %dist in NM at yf
x_nm = (xf-x0)*((x_low + x_hi)/2); %dist in NM from x0 to xf

vars = [N V Wx Wy x_nm y_nm];

slh = deg2rad(bth); %straight line heading = boat true heading from trig script, RADIANS
sldx = x_nm2; %NM
sldy = y_nm2; %NM
t_ub = 72; %guess for time upper bound, IN HOURS

%Bounds (1 x N+1)
lb = [0 -y_nm 0]; %[NM NM RAD]
ub = [x_nm y_nm 2*pi]; %[NM NM RAD]
lb_N = zeros(1,3*N+1);
ub_N = zeros(1,3*N+1);

for i = 1:N
    lb_N(1,1+3*(i-1):3+3*(i-1)) = lb;
    ub_N(1,1+3*(i-1):3+3*(i-1)) = ub;
end
lb_N(1,3*N+1) = 0; %time lower bound, HOURS
ub_N(1,3*N+1) = t_ub; %time upper bound, HOURS

%cost function
objective = @(x) objective(x);

%nonlinear constraints
nonlcon = @(x) nonlconZ(vars,x);

%linear constraints
A = [];
b = [];

% Aeq = [1,zeros(1,3*N);...
%        0,1,zeros(1,3*N-1);...
%        0,0,1,zeros(1,3*N-2);...
%        zeros(1,3*N+1);...
%        zeros(3*N-8,3*N+1);... %I bet this is wrong
%        zeros(1,3*N-3),1,0,0,0;...
%        zeros(1,3*N-2),1,0,0;...
%        zeros(3,3*N+1)]; %I bet this is wrong too

Aeq = [1,zeros(1,3*N);...
    0,1,zeros(1,3*N-1);...
    0,0,1,zeros(1,3*N-2);
    zeros(3*N-6,3*N+1);...
    zeros(1,3*N-3),1,0,0,0;...
    zeros(1,3*N-2),1,0,0;...
    zeros(1,3*N-1),1,0;...
    zeros(1,3*N+1)];

beq = [0;0;slh;zeros(3*N-6,1);x_nm;y_nm;slh;0]; %slh - straight line heading (i.e. if no wind/current) RADIANS
%                                                   %sldx/y = straight line
%                                                   dist in x/y **NM**
%                                                   %last position:
%                                                   estimated final time
%                                                   HOURS
%                                              

x0 = [0;0; slh; zeros(3*(N-2),1);x_nm;y_nm; slh; 48]; %initial guess

%fmincon options
options = optimoptions('fmincon','algorithm','sqp','TolFun',1e-8,'TolCon',1e-8);

%run fmincon
[x,fval,exitflag,output] = fmincon(objective,x0,A,b,Aeq,beq,lb_N,ub_N,nonlcon,options);

z = zeros(3,N);
    for q = 1:N
        z(:,q) = x((q-1)*3+1:(q-1)*3+3,1);
    end
%%
%data for plotting
t = zeros(N,1);
dx = zeros(N,1);
dy = zeros(N,1);
theta = zeros(N,1);

dt = x(end)/(N-1);

for i = 1:N
    dx(i) = x(1+3*(i-1));
    dy(i) = x(2+3*(i-1));
    theta(i) = x(3+3*(i-1));
    t(i) = dt*(i-1);
end 
t(N) = x(end);
theta = rad2deg(theta);
%Plot Results
plot(t,dx,'-o',t,dy,'-*',t,theta,'-p')
grid on

figure(200);clf
plot(dx,dy)
