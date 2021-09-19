function [Aeq,beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
n_all_poly = n_seg*(n_order+1);
%#####################################################
% p,v,a,j constraint in start,
Aeq_start = zeros(4, n_all_poly);
beq_start = zeros(4, 1);
% STEP 2.1: write expression of Aeq_start and beq_start
%
%
%
%
for i = 1:size(Aeq_start,1)
    Aeq_start(i,1:n_order+1) = getDerivative(0,n_order,i-1);
end
beq_start = start_cond';

%#####################################################
% p,v,a constraint in end
Aeq_end = zeros(4, n_all_poly);
beq_end = zeros(4, 1);
% STEP 2.2: write expression of Aeq_end and beq_end
%
%
%
%
for i = 1:size(Aeq_end,1)
    Aeq_end(i,n_all_poly-n_order:end) = getDerivative(ts(end),n_order,i-1);
end
beq_end = end_cond';

%#####################################################
% position constrain in all middle waypoints
Aeq_wp = zeros(n_seg-1, n_all_poly);
beq_wp = zeros(n_seg-1, 1);
% STEP 2.3: write expression of Aeq_wp and beq_wp
%
%
%
%
for i = 1:n_seg-1
    Aeq_wp(i,(i)*(n_order+1)+1:(i+1)*(n_order+1)) = getDerivative(0,n_order,0);
end
beq_wp = waypoints(2:end-1);
%#####################################################
% position continuity constrain between each 2 segments
Aeq_con_p = zeros(n_seg-1, n_all_poly);
beq_con_p = zeros(n_seg-1, 1);
% STEP 2.4: write expression of Aeq_con_p and beq_con_p
%
%
%
%
for i = 1:n_seg-1
    p1 = getDerivative(ts(i),n_order,0);
    p2 = getDerivative(0,n_order,0);
    Aeq_con_p(i,(i-1)*(n_order+1)+1:(i+1)*(n_order+1)) = [p1 -p2];
end
%#####################################################
% velocity continuity constrain between each 2 segments
Aeq_con_v = zeros(n_seg-1, n_all_poly);
beq_con_v = zeros(n_seg-1, 1);
% STEP 2.5: write expression of Aeq_con_v and beq_con_v
%
%
%
%
for i = 1:n_seg-1
    p1 = getDerivative(ts(i),n_order,1);
    p2 = getDerivative(0,n_order,1);
    Aeq_con_v(i,(i-1)*(n_order+1)+1:(i+1)*(n_order+1)) = [p1 -p2];
end
%#####################################################
% acceleration continuity constrain between each 2 segments
Aeq_con_a = zeros(n_seg-1, n_all_poly);
beq_con_a = zeros(n_seg-1, 1);
% STEP 2.6: write expression of Aeq_con_a and beq_con_a
%
%
%
%
for i = 1:1:n_seg-1
    p1 = getDerivative(ts(i),n_order,2);
    p2 = getDerivative(0,n_order,2);
    Aeq_con_a(i,(i-1)*(n_order+1)+1:(i+1)*(n_order+1)) = [p1 -p2];
end
%#####################################################
% jerk continuity constrain between each 2 segments
Aeq_con_j = zeros(n_seg-1, n_all_poly);
beq_con_j = zeros(n_seg-1, 1);
% STEP 2.7: write expression of Aeq_con_j and beq_con_j
%
%
%
%
for i = 1:1:n_seg-1
    p1 = getDerivative(ts(i),n_order,3);
    p2 = getDerivative(0,n_order,3);
    Aeq_con_j(i,(i-1)*(n_order+1)+1:(i+1)*(n_order+1)) = [p1 -p2];
end
%#####################################################
% combine all components to form Aeq and beq
Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
beq = [beq_start; beq_end; beq_wp; beq_con];
end

function derivative = getDerivative(ts,n_order,k)
    derivative = zeros(1,n_order+1);
    for i = k+1:n_order+1
        derivative(1,i) = factorial(i-1)*ts^(i-1-k)/factorial(i-1-k);
    end
end