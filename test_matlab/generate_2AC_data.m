function [Image1, Image2, At, R_extrinsic_para, T_extrinsic_para, R_gt, t_gt] = generate_2AC_data(dof_type,match_type)

%% Two types of affine correspondence
%% inter-cam affine correspondences refer to correspondences which are seen by different cameras over two consecutive views, which contain
% (1) 1 affine correspondence: (cam1 at time i) <--> (cam2 at time j)
% (2) 1 affine correspondence: (cam2 at time i) <--> (cam1 at time j)
%% intra-cam affine correspondences refer to correspondences which are seen by the same camera over two consecutive views, which contain
% (1) 1 affine correspondence: (cam1 at time i) <--> (cam1 at time i)
% (2) 1 affine correspondence: (cam2 at time i) <--> (cam2 at time j)

n_cam = 2;
n_point = 2;

%% define random extrinsic parameters
cam_body_rotation = cell(n_cam, 1);
cam_body_offset = cell(n_cam, 1);
R_cam = cell(n_cam, 1);
t_cam = cell(n_cam, 1);
T_body_cam = cell(n_cam, 1);
for ii = 1:n_cam
    cay = rand(3, 1);
    cam_body_rotation{ii} = cayley_rotation(cay);
    cam_body_offset{ii} = rand(3, 1);
    
    R_cam{ii} = cam_body_rotation{ii}';
    t_cam{ii} = -R_cam{ii}*cam_body_offset{ii};
    % transformation from body reference to perspective camera references
    T_body_cam{ii} = [R_cam{ii} t_cam{ii}; 0 0 0 1];
end


%% define relative pose
if strcmp(dof_type, '3DOF')
    % planar motion
    Ay_degree = 10*rand(1);
    R_gt = [ cosd(Ay_degree)   0   -sind(Ay_degree);
        0      1     0;
        sind(Ay_degree)   0   cosd(Ay_degree)];
    
    At_degree = 10*rand(1);
    t_gt = [sind(At_degree); 0; -cosd(At_degree)];
    
elseif strcmp(dof_type, '4DOF')
    % motion with known vertical direction
    Ay_degree = 10*rand(1);
    R_gt = [ cosd(Ay_degree)   0   -sind(Ay_degree);
        0      1     0;
        sind(Ay_degree)   0   cosd(Ay_degree)];
    
    t_gt = rand(3, 1);
    
else
    error('unsupported dof type!')
end

% transformation from body reference at time i to time j
T_gt = [R_gt t_gt; 0 0 0 1];

%% generating affine correspondences
% generating random scene points
[PT, Distance, Nvector] = generate_3Dscenepoints(n_point);
points_all = cell(n_point, 1);
for ii = 1:n_point
    points_all{ii} = struct('point', PT(:,:,ii));
end

%% extract point observations
% images at time i
x_i = cell(n_point, n_cam);
for ii = 1:n_point
    PT = points_all{ii}.point;
    for jj = 1:n_cam
        tmp = R_cam{jj}*PT+t_cam{jj};
        x_i{ii,jj} = tmp./tmp(3,:);
    end
end
% images at time j
Rc_j = cell(n_cam, 1);
tc_j = cell(n_cam, 1);
for ii = 1:n_cam
    tmp = T_body_cam{ii}*T_gt;
    Rc_j{ii} = tmp(1:3,1:3);
    tc_j{ii} = tmp(1:3,4);
end
x_j = cell(n_point, n_cam);
for ii = 1:n_point
    PT = points_all{ii}.point;
    for jj = 1:n_cam
        tmp = Rc_j{jj}*PT+tc_j{jj};
        x_j{ii,jj} = tmp./tmp(3,:);
    end
end

%% construct observations
R_extrinsic_para = cat(3, cam_body_rotation{1}, cam_body_rotation{2});
T_extrinsic_para = [cam_body_offset{1}, cam_body_offset{2}];
Image1 = zeros(3, 2);
Image2 = zeros(3, 2);
At = zeros(2, 2, 2);
if strcmp(match_type, 'inter')
    % compute an AC using the frist plane
    Image1(:,1) = x_i{1,1}(:,1);
    Image2(:,1) = x_j{1,2}(:,1);
    H_c1i_c2j = DLT_homography(x_i{1,1}(:,2:end)',x_j{1,2}(:,2:end)');
    At(:,:,1) = get_affinetransformation(H_c1i_c2j, Image1(:,1), Image2(:,1));
    % compute the other AC using the second plane
    Image1(:,2) = x_i{2,2}(:,1);
    Image2(:,2) = x_j{2,1}(:,1);
    H_c2i_c1j = DLT_homography(x_i{2,2}(:,2:end)',x_j{2,1}(:,2:end)');
    At(:,:,2) = get_affinetransformation(H_c2i_c1j, Image1(:,2), Image2(:,2));
    
elseif strcmp(match_type, 'intra')
    % compute an AC using the frist plane
    Image1(:,1) = x_i{1,1}(:,1);
    Image2(:,1) = x_j{1,1}(:,1);
    H_c1i_c1j = DLT_homography(x_i{1,1}(:,2:end)',x_j{1,1}(:,2:end)');
    At(:,:,1) = get_affinetransformation(H_c1i_c1j, Image1(:,1), Image2(:,1));
    % compute the other AC using the second plane
    Image1(:,2) = x_i{2,2}(:,1);
    Image2(:,2) = x_j{2,2}(:,1);
    H_c2i_c2j = DLT_homography(x_i{2,2}(:,2:end)',x_j{2,2}(:,2:end)');
    At(:,:,2) = get_affinetransformation(H_c2i_c2j, Image1(:,2), Image2(:,2));
else
    error('unsupported match type!')
end



