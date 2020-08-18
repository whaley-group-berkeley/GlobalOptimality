%Find Globally optimal protocols
% PARAMETERS
dt  = 0.01; % WAS 0.0001     % Time step
T   = 3;      % Total simulation time
dr  = 0.001;     % Spacing of r in grid (0.0001)
Rcut=0.999;
us  = [0,1];   % Set of u values to try
k   = 1;         % measurement rate
eta = 0.3;         % quantum efficiency
rand_seed=222;   % seed for random number generation
sigma_cutoff = dr;    % Smallest allowed standard dev. of Guassian (avoid undersampling)
cost_func = 'r';
% Initialize arrays for u, r and t
num_us = length(us);
rs = 0:dr:Rcut; num_rs = length(rs);
[Us, Rs] = meshgrid(us, rs);
timesteps = T/dt;
rng(rand_seed);

% Provide a stochastic equation for r: dr1 = A dt + B dW
% A = k.*(Rs-eta./Rs).*(Us.^2-1);  % Qubit purification with feedback
rnew = sqrt(exp(2*k*(-1)*dt).*(rs.^2-eta)+eta); %rnew for A

% Set up the simulation:
% Cost function for optimization
if strcmp(cost_func, 'r')

    cost = 1-rs;

elseif strcmp(cost_func, 'purity')
    cost = (1-rs.^2)/2; % 1-(1+r^2)/2
end
rtsmatrix=repmat(rs, [length(rs),1]);
rsB=-Rcut:dr:Rcut;
rtsma=repmat(rs, [length(rsB),1]);

num_rsB = length(rsB);
rtsmatrixB=repmat(rsB, [length(rs),1]);
% r0smatrix=rtsmatrix'; % Removed redundant variable
%% Computing P(r(t+dt)|r(t), u) for all us
normedPA=zeros(num_rs,num_rs); %u, r0, r'
normedPB=zeros(num_rs,num_rsB); 
%     distr = exp(-( rtsmatrix-A1.*dt-r0smatrix).^2 ./(dt*2.*B1.^2)) + ...
%             exp(-(-rtsmatrix-A1.*dt-r0smatrix).^2 ./(dt*2.*B1.^2));
    % replace r with |r|
    r1 = repmat(rnew, [num_rs,1]);
    distrA= exp(-((r1'-rtsmatrix).^2./2/(dr)^2))/(sqrt(2*pi)*dr);
    zm=(rtsmatrixB-rtsma')./(1-rtsmatrixB.*(rtsma'));
    W=atanh(zm)./sqrt(2*eta*k);
    De=(1-rtsma'.^2)./(sqrt(2*eta*k)*(1-zm.^2).*(rtsma'.*rtsmatrixB-1).^2);
    distrB= exp(-W.^2/(2*dt)-k*eta*dt).*(cosh(sqrt(2*k*eta)*W)+sinh(sqrt(2*eta*k)*W).*rtsma').*abs(De)/sqrt(2*pi*dt);


    weightA=sum(distrA,2);
    normalize = repmat(weightA,1,num_rs);
    normedPA( :, :) = (distrA./(normalize));
    weightB=sum(distrB,2);
    normalize = repmat(weightB,1,num_rsB);
    normedPB( :, :) = (distrB./(normalize));

%% Free up some memory and propogate backwards in time, finding the optimal protocol u at each time step and r(t)
% clear A1 B1 distr rtsmatrix normalize % Free up some memory if needed

optUs=zeros(timesteps,num_rs);
costs=zeros(timesteps,num_rs);
allcost=zeros(num_us,num_rs,timesteps);
for ti=1:1:timesteps
    comparecost=zeros(num_us, num_rs);
    % Compute the expected value of the cost function for each value of u and r
       % ip=cost * ((squeeze(normedP(ui,:,:)))'); % Inner product
        comparecost(1,:)=cost * (normedPA(:,:)');
        allcost(1,:,ti)=comparecost(1, :);
        comparecost(2,:)=cost * (normedPB(:,num_rs:end)')+cost(2:end) * (normedPB(:,num_rs-1:-1:1)');
        allcost(2,:,ti)=comparecost(2, :);
   
    [newcost, optUi] = min(comparecost,[],1); % find the optimal u for each r0
    cost=newcost;
    costs(ti,:)=newcost; % Uncomment this to save the cost, then compare
    %with trajectories
    optUs(ti,:)=us(optUi);
end
%%
imagesc( [0,T],[0,1],rot90(optUs,3)); xlabel('t','FontName','Times New Roman','FontSize',14,'FontWeight','bold'); ylabel('r','FontName','Times New Roman','FontSize',14,'FontWeight','bold')
set(gca, 'YDir', 'normal')
set(gca,'linewidth',1);
colormap(colormap('gray'))

