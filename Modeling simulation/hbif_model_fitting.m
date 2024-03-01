function [xs,fitting] = hbif_model_fitting(FC_emp, C, Tmax, we, nSubs, nNodes, a1, TR, omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FINAL SIMULATION
fprintf(1, 'SIMULATING OPTIMIZED MODEL.\n');
a = a1; %use those avalues that have been found to be optimal

%wC = we*C;
%sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
%xs=zeros(10000/2,nNodes);
xs=zeros(Tmax*nSubs,nNodes);


%FROM HERE ON SIMULATIONS AND FITTING
dt=0.1.*(TR/2);  %BEFORE dt = 0.1;
sig = 0.04; %was 0.04
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step
wC = we*C;
sumC = repmat(sum(wC,2),1,2); % for sum Cij*x
z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
nn=0;

for t=1:dt:3000 %JVS is it really necessary to swing in for 3000secs?
    suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
    zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
    z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
end

for t=1:dt:Tmax*TR*nSubs %JVS: was 15000, now faster
    suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
    zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
    z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
    
    if (abs(mod(t,TR))<0.01)%%%% BEFORE mod(t,TR)==0
        nn=nn+1;
        xs(nn,:)=z(:,1)';
    end
end

fprintf(1, 'COMPUTING MODEL FIT.\n');
FC_simul = corrcoef(xs(1:nn,:)); %Now one FC_simul per G
cc=corrcoef(atanh(squareform(tril(FC_emp,-1))),atanh(squareform(tril(FC_simul,-1))));%atanh(FC...
fitting=cc(2);

