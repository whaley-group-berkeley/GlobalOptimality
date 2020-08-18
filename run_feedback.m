%Simulate trajectories/u=0,i.e. feedback always on
eta = 0.3; %efficiency
num_runs = 10000;%number of total runs
ndt = 0.001;%dt
T = 1;%total time
k = 1;  %measurement strength
ntimestep = T/ndt;
sigmaz = diag([1,-1]);
sigmax = [0 1;1 0];
rho = zeros(2,2,ntimestep,num_runs);
drho = zeros(2,2);
rho(:,:,1,:) = repmat([1/2,0;0, 1/2],[1,1,1,num_runs]);%initial value of rho
r = zeros(ntimestep,num_runs);
r(1,:) = repmat(0,[1,num_runs]);%initial value of r
dW = sqrt(ndt)*randn(ntimestep,num_runs);
X = sigmaz/2; %measurement operator
for ii = 1:num_runs
    for kk = 1:ntimestep
        
        %W=exp(-k*dt)/sqrt(2*pi*dt)*rho(1,1,kk,ii)
   
    
    drho = -k*(X^2*rho(:,:,kk,ii)+rho(:,:,kk,ii)*X^2-2*X*rho(:,:,kk,ii)*X)*ndt...
        +dW(kk,ii)*sqrt(2*eta*k)*(X*rho(:,:,kk,ii)+rho(:,:,kk,ii)*X-2*trace(X*rho(:,:,kk,ii))*rho(:,:,kk,ii));
    rho(:,:,kk+1,ii)=rho(:,:,kk,ii)+drho;
    r(kk+1,ii) = sqrt(trace(rho(:,:,kk+1,ii)*sigmaz)^2+trace(rho(:,:,kk+1,ii)*sigmax)^2);
    rho(:,:,kk+1,ii) = (eye(2)+r(kk+1,ii)*sigmax)/2;
    if r(kk+1,ii)>1
        r(kk+1:ntimestep,ii) = ones(ntimestep-kk,1);
        break
    end
    
    end
  ii
end

rav_f = mean(r(ntimestep,:));
rr_f = mean(r(1:end-1,:),2);




        