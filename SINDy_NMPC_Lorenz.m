%%% Nonlinear System Identification and Nonlienar Model Predictive Control
%%% Allen Lee

%%% In this task, Lorenz system is used for its nonlinear and chaotic dynamcis,
%%% and SINDy and NMPC are used for system identification and control.
clc
clear all
 
global N n m x xr ur Eta timestep
global alpha Rho beta
%%% System Parameters
alpha = 10;Rho = 28;beta = 8/3;
n = 3;m = 2;

x0 = [-8 7 27]';% IC
timestep = 0.001; % too large, blow up!
RunTime = 10;
Sample_num = round(RunTime/timestep); % #data sample collected in time

X_noisy = false; % if X is noisy, poor performance
X_dot_noisy = true;% this term is rather robust to noise
U_noisy = false; % if input is noisy

%%% Data Generation
X_dot_clean = zeros(Sample_num,n); % record x_dot
X_clean = zeros(Sample_num,n); % record x
X_clean(1,:) = x0;
U_clean = zeros(Sample_num,m);
for i = 1:Sample_num
    U_clean(i,:) = -5+10.*rand(m,1);
    X_dot_clean(i,:) = Lorentz(X_clean(i,:),U_clean(i,:));
    
    if(i<Sample_num)
        X_clean(i+1,:) = X_clean(i,:) + X_dot_clean(i,:) .*timestep;
    end

end
%%% Corrupted Measurement Data
if(X_dot_noisy)
    X_dot = X_dot_clean + wgn(size(X_dot_clean,1),size(X_dot_clean,2),0); %0 dBW
else
    X_dot = X_dot_clean;
end
if(X_noisy)
    X = X_clean + wgn(size(X_dot_clean,1),size(X_dot_clean,2),-10);
else
    X = X_clean;
end
if(U_noisy)
    U = U_clean + wgn(size(U_clean,1),size(U_clean,2),-10);
else
    U = U_clean;
end

%%
save("LorenzwithControlData.mat","X","X_clean","U","X_dot","X_dot_clean","x0")
%%
%%% System Identification using SINDy

% Candidates Functions
Theta = Candidate_Library(X,U);
thresthold = 0.01; % as you see fit
Eta = pinv(Theta)*X_dot;
Eta_new = Eta;

%%% Apply STLS to find best-fit
while (true)
    % find big enough coefficients
    biginds = (abs(Eta) >= thresthold);
    % Threshold out small coefficients to eliminate from library
    Eta_new(~biginds) = 0;

    % Find coefficients for reduced library
    for i = 1:n  % one dynamics by one
        dummy = Eta_new(:,i);
        dummy(biginds(:,i)) = pinv(Theta(:,biginds(:,i)))*X_dot(:,i);
        Eta_new(:,i) = dummy;
    end

    if(sum(sum(abs(Eta - Eta_new))) == 0 ) % converged
        break
    else
        Eta = Eta_new;
    end
end
Eta;
sum(Eta~=0);


%%% Data Restoration
%%%
X_dot_rcv = zeros(size(X_dot));
X_rcv = zeros(size(X));
X_rcv(1,:) = x0; % Assume we know IC
theta_rcv = zeros(size(Theta(1,:)));
for i = 1:Sample_num

    theta_rcv = Candidate_Library(X_rcv(i,:),U_clean(i,:));
    X_dot_rcv(i,:) = theta_rcv*Eta;

    if(i<Sample_num)
        X_rcv(i+1,:) = X_rcv(i,:) + X_dot_rcv(i,:).*timestep;
    end
    
end
rmse(X_dot_clean,X_dot_rcv)
rmse(X_clean,X_rcv)

%% Data Comparison
%%%
rmse(X_dot_clean,X_dot_rcv)
rmse(X_clean,X_rcv)


for i=1:n
    figure
    hold on
    plot((0:Sample_num-1)*timestep,X_clean(:,i));
    plot((0:Sample_num-1)*timestep,X_rcv(:,i));
    plot((0:Sample_num-1)*timestep,X(:,i));
    legend('true','recoverd','noisy')
    hold off
end

for i=1:n
    figure
    hold on
    plot((0:Sample_num-1)*timestep,X_dot_clean(:,i));
    plot((0:Sample_num-1)*timestep,X_dot_rcv(:,i));
    plot((0:Sample_num-1)*timestep,X_dot(:,i));
    legend('true','recoverd','noisy')
    hold off
end
%% 3D plot
figure
hold on
view(45,15)

plot3(X_clean(:,1),X_clean(:,2),X_clean(:,3),'b.')
plot3(X_rcv(:,1),X_rcv(:,2),X_rcv(:,3),'r.')
plot3(X(:,1),X(:,2),X(:,3),'y.')

hold off
legend('true','recovered','noisy')
xlabel('X')
ylabel('Y')
zlabel('Z')

%% NMPC 
N = 10; % Control/Prediction Horizen
xr = zeros(n,N+1); % State Reference 
ur = zeros(m,N); % Manual

Run_Sample = 1e4*0.5;
x0 = [-3.77 5.69 7.02]';
u0 = zeros(m,1);
Xall = zeros(n,Run_Sample+1); Uall = zeros(m,Run_Sample);
Xr = zeros(n,size(Xall,2)+N+1); Ur = zeros(m,size(Uall,2)+N);

%%% Reference Tracking %%%
dummytimeline = (0:size(Xr,2)-1)*timestep;
Xr(1,:) = 10.*sin(1.*dummytimeline);
Xr(2,:) = 10.*cos(1.*dummytimeline);
for i=1:size(Ur,2)
    Ur(:,i) =  EquilibriumPoint(Xr(:,i));
end

% Xall = zeros(n,Run_Sample+1); Uall = zeros(m,Run_Sample);
Xall(:,1) = x0; Uall(:,1) = u0; x=x0; u=u0;
ur = Ur(:,1:1+N);
xr = Xr(:,1:1+N+1);

Xk = zeros(n*(N+1),1); % all future x from step k
Uk = zeros(m*N,1); % all future u from step k
Xk(1:n,1) = x0;
Uk(1:m,1) = u0;
zek = [Xk;Uk];
zek(1:n,1) = x0-xr(:,1); % Difference from the references
zek((N+1)*n+1:(N+1)*n+m) = u0-ur(:,1);

% Input & Output Constraints
xmin = [-10.1;-10.1;-100];
xmax = [10.1;10.1;100];
umin = [-100;-100].*1e2;
umax = [100;100].*1e2;

Fx = [eye(n);-eye(n)];
Fu = [eye(m);-eye(m)];

R = eye(m);
Qx = 1e6.*diag([1 1 0]); % Don't care 3rd state
QX = zeros(n*(N+1),n*(N+1)); RU = zeros(m*N,m*N);
FX = zeros(2*n*(N+1),n*(N+1)); FU = zeros(2*m*N,m*N);
Gxe = zeros(2*n*(N+1),1); Gue = zeros(2*m*N,1);
for i = 1:N
    QX((i-1)*n+1:i*n,(i-1)*n+1:i*n) = Qx;
    RU((i-1)*m+1:i*m,(i-1)*m+1:i*m) = R;

    FX((i-1)*2*n+1:i*2*n,(i-1)*n+1:i*n) = Fx;
    FU((i-1)*2*m+1:i*2*m,(i-1)*m+1:i*m) = Fu;
end

QX(N*n+1:end,N*n+1:end) = Qx;
FX(N*2*n+1:end,N*n+1:end) = Fx;
H = blkdiag(QX,RU.*0);

% Simulation with NMPC
ObjFunc = @(ze) (ze'*H*ze);
F = blkdiag(FX,FU);
Feq = []; Geq = []; lb = []; ub = []; % Linear constraints

nonlcon = @NonlinearConstraints;
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
for i = 1:Run_Sample
    x = Xall(:,i); 
    ur = Ur(:,i:i+N);
    xr = Xr(:,i:i+N+1);
    
    Gxe = zeros(2*n*(N+1),1);
    Gue = zeros(2*m*N,1);
    for j = 1:N
        Gxe((j-1)*2*n+1:j*2*n,1) = [xmax;-xmin] - Fx*xr(:,j);
        Gue((j-1)*2*m+1:j*2*m,1) =  [umax;-umin] - Fu*ur(:,j);
    end
    Gxe(N*2*n+1:end,1) = [xmax;-xmin] - Fx*xr(:,end);
    G = [Gxe;Gue];

    zek = fmincon(ObjFunc,zek,F,G,Feq,Geq,lb,ub,nonlcon,options);
    Uall(:,i) = zek(n*(N+1)+1:n*(N+1)+m,1) + ur(:,1);
    x_dot = Lorentz(Xall(:,i),Uall(:,i));
    if(i<=Run_Sample)
        Xall(:,i+1) =  Xall(:,i) +  x_dot.*timestep;
    end

end
%% Save and Figures
save("SINDy_NMPC_Lorentz.mat","Xall","Xr","Uall","Ur")
rmse(Xall(1:2,:)',Xr(1:2,1:size(Xall,2))')
%% Visualization and Errors
timeframe = 0:timestep:(size(Xall,2)-1)*timestep;
figure
hold on
plot(timeframe,Xall(1,:),'r.')
plot(timeframe,Xr(1,1:size(Xall,2)),'b')
hold off
legend("Model","True","Noise")


figure
hold on
plot(timeframe,Xall(2,:),'r.')
plot(timeframe,Xr(2,1:size(Xall,2)),'b')
% plot(Xg(2,:),'g')
hold off
legend("Model","True","Noise")

figure
hold on
plot(timeframe,Xall(3,:),'r.')
plot(timeframe,Xr(3,1:size(Xall,2)),'b')
% plot(Xg(3,:),'g')
hold off
legend("Model","True","Noise")

figure
hold on
plot(Uall(1,:),'r.')
plot(Uall(2,:),'b.')
hold off
legend("u1","u2")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ur = EquilibriumPoint(xr)
global alpha Rho
ur = zeros(2,1);
ur(1) = nthroot(-21.6*alpha*(xr(2)-xr(1)),3);
ur(2) = -1*(Rho*xr(1)-xr(1)*xr(3)-xr(2));
end

function [c,ceq] =NonlinearConstraints(z)
global N n m x xr ur Eta timestep
c=[];
shfu = (N+1)*n;
for i=1:N
    x_dummy = z((i-1)*n+1:i*n,1)+xr(:,i);
    u_dummy = z(shfu+(i-1)*m+1:shfu+i*m)+ur(:,i);
    x_prime =(Candidate_Library(x_dummy',u_dummy')*Eta)';
    c1((i-1)*n+1:i*n,1)=z(i*n+1:(i+1)*n,1)+xr(:,i+1) - timestep.*x_prime - x_dummy;

end
ceq = [z(1:n)-x+xr(:,1);c1];
end

function x_dot = Lorentz(x,u)
global alpha Rho beta
x_dot = [alpha*(x(2)-x(1)) + 1/21.6*u(1)^3;
    Rho*x(1)-x(1)*x(3)-x(2) + u(2);
    x(1)*x(2)-beta*x(3)];
end