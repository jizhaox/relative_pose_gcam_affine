function [cay_sols, t_sols, R_sols] = solver_depth_intra_2ac_enhance(Image1, Image2, At, R_cam, t_cam)

% add coordinate transformation
% according to our evaluation, it has better numerical stability

[R0, tr0] = find_optimal_transformation(t_cam(:,1), t_cam(:,2));
T0 = [R0, tr0; 0, 0, 0, 1];
inv_T0 = [R0', -R0'*tr0; 0, 0, 0, 1];
Tf_cam_1 = T0 * [R_cam(:,:,1), t_cam(:,1); 0, 0, 0, 1];
Tf_cam_2 = T0 * [R_cam(:,:,2), t_cam(:,2); 0, 0, 0, 1];
R_cam_new = cat(3, Tf_cam_1(1:3,1:3), Tf_cam_2(1:3,1:3));
t_cam_new = [Tf_cam_1(1:3,4); Tf_cam_2(1:3,4)];

[cay_sols_tmp, t_sols_tmp, R_sols_tmp, ~] = solver_depth_intra_2ac(Image1, Image2, At, R_cam_new, t_cam_new);

n_sol = size(cay_sols_tmp, 2);
cay_sols = zeros(size(cay_sols_tmp));
t_sols = zeros(size(t_sols_tmp));
R_sols = zeros(size(R_sols_tmp));
for ii = 1:n_sol
    q_sol = cay_sols_tmp(:, ii);
    t_sol = t_sols_tmp(:, ii);
    R_sol = cayley_rotation(q_sol);
    T_sol = inv_T0*[R_sol, t_sol; 0, 0, 0, 1]*T0;
    
    R = T_sol(1:3, 1:3);
    q = rotm2quat(R);
    R_sols(:, :, ii) = R;
    cay_sols(:, ii) = q(2:4)/q(1);
    t_sols(:, ii) = T_sol(1:3, 4);
end