clear, clc
%%%%%%%%%%%%%%

PFlag = 1;
N = 200; % phenotype space
B0 = 1; m = 0.15;
P0 = 0.5;
a = -0.1;
ppos = 0.05;
POP = 1e4;

[v_av,sig,sig2] = simmer(PFlag,N,B0,m,P0,a,ppos,POP);
v_temp = v_av;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_av,wid,sig2]= simmer(PFlag,N,B0,m,P0,a,ppos,POP)
% general parameters
M = 1e6; % number of timesteps
dt = 1e-3; % timestep
x1 = 50; % start recording
x2 = 150; % end recording
%%%%%
x = (1:N)'; % phen coord
t = dt*(0:M); % time coord
% mutations
pneg = 1-ppos;
%%%%%
out = zeros(N,M);
out(100,1) = 1;
x_av_vec = zeros(M,1);
x_av_vec(1) = 100;

wid = zeros(M-1,1);
sig2 = zeros(M-1,1);

for i = 1:(M-1)
    pop_in = out(:,i); % input vector
    x_av = sum(x.*pop_in)./sum(pop_in); % average phen
    x_av_vec(i+1) = x_av;
    B = B0 + m*(x-x_av);
    
    ird = find(pop_in > 1/POP,1,'last');
    wid(i) = x(ird)-x_av;
    sig2(i) = sum(pop_in.*(x-x_av).^2)./sum(pop_in);

    % only participate if density is high enough
    B(pop_in < 1/POP) = 0;
    Bav = sum(B.*pop_in)./sum(pop_in); % average birth rate
    P = P0 + a*(x-x_av);
    P(P>1) = 1; P(P<0) = 0;
    % fluxes
    dBirth = dt*(B - Bav).*pop_in;
    dmut = zeros(N,1); 
    Nmut = P.*B.*pop_in;
    dmut(2:end-1) = -Nmut(2:end-1) + ppos*Nmut(1:end-2) + pneg*Nmut(3:end);
    dmut = dt*dmut;
    out(:,i+1) = out(:,i) + (dBirth + dmut);

    if ~mod(i+1,1e3) && PFlag == 1
        figure(1);clf; 
        bar(out(:,i+1),'k')
        ylim([0,1.2*max(out(:,i+1))])
        drawnow;
    end

    if x_av < x1 || x_av > x2
        break
    end
end
ind = min([1e3,round(i/3)]);

wid = mean(wid(ind:i));
sig2 = mean(sig2(ind:i));

x_av_vec = x_av_vec(1:i+1);
t_vec = t(1:i+1);
v_vec = gradient(x_av_vec,t_vec);
v_av = mean(v_vec(ind:i+1));
%%%%%%%
end









