function [Image1, Image2, At, R_gt, cay_gt, t_gt] = generate_2AC_single_data(n_point)

%% Two affine correspondences of a single camera
n_view = 2;

%% define relative pose
cay = rand(3, 1);
R_gt = cayley_rotation(cay);
q = rotm2quat(R_gt);
cay_gt = q(2:4)/q(1);
cay_gt = cay_gt(:);

t_gt = rand(3, 1);

% transformation from perspective camera reference at time i to time j
T_gt = [R_gt t_gt; 0 0 0 1];

% camera pose at time i
T_cami = eye(4,4);
R_cam{1} = T_cami(1:3,1:3);
t_cam{1} = T_cami(1:3,4);

% camera pose at time j
T_camj = T_gt;
R_cam{2} = T_camj(1:3,1:3);
t_cam{2} = T_camj(1:3,4);

%% generating affine correspondences
% generating random scene points
[PT, Distance, Nvector] = generate_3Dscenepoints(n_point);
points_all = cell(n_point, 1);
for ii = 1:n_point
    points_all{ii} = struct('point', PT(:,:,ii));
end

%% extract point observations
x_i = cell(n_point, n_view);
for ii = 1:n_point
    PT = points_all{ii}.point;
    for jj = 1:n_view
        tmp = R_cam{jj}*PT+t_cam{jj};
        x_i{ii,jj} = tmp./tmp(3,:);
    end
end

%% construct observations
Image1 = zeros(3, 2);
Image2 = zeros(3, 2);
At = zeros(2, 2, 2);

% compute an AC using the frist plane
Image1(:,1) = x_i{1,1}(:,1);
Image2(:,1) = x_i{1,2}(:,1);
H(:,:,1) = DLT_homography(x_i{1,1}(:,2:end)',x_i{1,2}(:,2:end)');
At(:,:,1) = get_affinetransformation( H(:,:,1), Image1(:,1), Image2(:,1));

% compute the other AC using the second plane
Image1(:,2) = x_i{2,1}(:,1);
Image2(:,2) = x_i{2,2}(:,1);
H(:,:,2) = DLT_homography(x_i{2,1}(:,2:end)',x_i{2,2}(:,2:end)');
At(:,:,2) = get_affinetransformation(H(:,:,2), Image1(:,2), Image2(:,2));
end