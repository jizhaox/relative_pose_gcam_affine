function sol = find_solution_t_singlecamera(sols, gt)

num_sol = size(sols, 2);
if (num_sol < 1)
    sol = [];
    return;
end

for i = 1:num_sol
    sols(:,i) = sols(:,i)/norm(sols(:,i))* norm(gt);
end

err = sols - repmat(gt, [1 num_sol]);
err_s = sum(abs(err), 1);
[~, idx] = min(err_s);
sol = sols(:, idx);