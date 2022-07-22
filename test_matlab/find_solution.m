function sol = find_solution(sols, gt)

num_sol = size(sols, 2);
if (num_sol < 1)
    sol = [];
    return;
end
err = sols - repmat(gt, [1 num_sol]);
err_s = sum(abs(err), 1);
[~, idx] = min(err_s);
sol = sols(:, idx);