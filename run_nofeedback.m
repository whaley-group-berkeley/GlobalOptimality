%Simulate trajectories without feedback;u=1;
%parameters same with File run_feedback
eta = 0.3;
num_runs = 10000;
ndt = 0.001;
T = 3;
ntimestep = T/ndt;
sigmaz = diag([1,-1]);
sigmax = [0 1;1 0];
rho = zeros(2,2,ntimestep,num_runs);
drho = zeros(2,2);
rho(:,:,1,:) = repmat([1/2,0;0, 1/2],[1,1,1,num_runs]);%initial value of rho
r = zeros(ntimestep,num_runs);
r(1,:) = repmat(0,[1,num_runs]);%initial value of r
dW = sqrt(ndt)*randn(ntimestep,num_runs);
X = sigmaz/2;
for ii = 1:num_runs
    for kk = 1:ntimestep
        
        %W=exp(-k*dt)/sqrt(2*pi*dt)*rho(1,1,kk,ii)
   
    
    drho = -k*(X^2*rho(:,:,kk,ii)+rho(:,:,kk,ii)*X^2-2*X*rho(:,:,kk,ii)*X)*ndt...
        +dW(kk,ii)*sqrt(2*eta*k)*(X*rho(:,:,kk,ii)+rho(:,:,kk,ii)*X-2*trace(X*rho(:,:,kk,ii))*rho(:,:,kk,ii));
    rho(:,:,kk+1,ii) = rho(:,:,kk,ii)+drho;
    r(kk+1,ii) = sqrt(trace(rho(:,:,kk+1,ii)*sigmaz)^2);
    if r(kk+1,ii)>1
        r(kk+1:ntimestep,ii) = ones(ntimestep-kk,1);
        break
    end
    
    end
  ii
end

rav_nf = mean(r(ntimestep,:));
rr_nf = mean(r(1:end-1,:),2);



        