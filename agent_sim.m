clear; clc
% params
N = 1e4;
M = 1e5;
dt = 3e-3;
B0 = 1;
m = 0.15;
p0 = 0.5;
pb = 0.2;
ome = 0;
% initialize agents
% [x, N+, N-, lineage ID]
ags = zeros(N,4);
ags(:,4) = (1:N)';
x_av_vec = zeros(M,3);
frac = zeros(M,1);
N_nose = zeros(M,2);
x_nose = 0;
for i = 1:M
    x_nose_temp = max(ags(:,1));
    if x_nose_temp > x_nose
        x_nose = x_nose_temp;
    end
    N_nose(i,:) = [length(find(ags(:,1) == x_nose)), x_nose];
    x_av = mean(ags(:,1));
    x_av_vec(i,:) = [x_av, mean(ags(:,2)), mean(ags(:,3))];
    temp = zeros(2*N,4);
    cc = 0;
    for j = 1:N
        runner = diver(ags(j,:),pb,p0,ome,B0,m,x_av,dt);
        cc2 = cc + length(runner(:,1));
        cc = cc + 1;
        temp(cc:cc2,:) = runner;
        cc = cc2;
    end
    ags = datasample(temp(1:cc2,:),N);
    if x_av_vec(i,2) + x_av_vec(i,3) > 150
        break
    end
end
M2 = i;
ind = find(x_av_vec(:,2) + x_av_vec(:,3) > 30,1);
x_av_vec = x_av_vec(ind:M2,:);
N_nose = N_nose(ind:M2,:);
% velocity
v_x = (x_av_vec(end,1)-x_av_vec(1,1))/(dt*(M2-ind));
v_Npos = (x_av_vec(end,2)-x_av_vec(1,2))/(dt*(M2-ind));
v_Nneg = (x_av_vec(end,3)-x_av_vec(1,3))/(dt*(M2-ind));
v_all = [v_x, v_Npos, v_Nneg];
output = {v_all,N_nose};

%%%%%%%
function offsp = diver(ag_in,pb,p0,ome,B0,m,x_av,dt)
p = p0 + ome*(ag_in(1) - x_av);
B = B0 + m*(ag_in(1) - x_av);

offsp = ag_in;
if rand(1) < B*dt
    temp = ag_in;
    if rand(1) < p
        switch rand(1) < pb
            case 0 % deleterious
                temp(1) = temp(1) - 1;
                temp(3) = temp(3) + 1;
            case 1 % beneficial
                temp(1) = temp(1) + 1;
                temp(2) = temp(2) + 1;
        end
    end
    offsp = [offsp; temp];
end
end
