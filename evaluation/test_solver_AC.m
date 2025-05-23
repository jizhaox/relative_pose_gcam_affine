clear;
%addpath('../mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1AC plane solver & inter-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------1AC plane solver & inter-camera case--------------')
%% generate synthetic data
[Image1, Image2, At_33, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('3DOF', 'inter');

%% run inter-1AC plane solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose3d_1AC_plane(Image1(:,1), Image2(:,1), At_33(:,:,1), R_cam, t_cam, 'inter');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC plane solver & inter-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC plane solver & inter-camera case--------------')
%% generate synthetic data
[Image1, Image2, At_33, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('3DOF', 'inter');

%% run inter-2AC plane solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose3d_2AC_plane(Image1, Image2, At_33, R_cam, t_cam, 'inter');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC plane solver & intra-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC plane solver & intra-camera case--------------')
%% generate synthetic data
[Image1, Image2, At_33, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('3DOF', 'intra');

%% run intra-2AC plane solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose3d_2AC_plane(Image1, Image2, At_33, R_cam, t_cam, 'intra');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC vertical solver & inter-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC vertical solver & inter-camera case--------------')
%% generate synthetic data
[Image1, Image2, At_33, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('4DOF', 'inter');

%% run 2AC vertical solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose4d(Image1, Image2, At_33, R_cam, t_cam, 'inter');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC vertical solver & intra-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC vertical solver & intra-camera case--------------')
%% generate synthetic data
[Image1, Image2, At_33, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('4DOF', 'intra');

%% run 2AC vertical solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose4d(Image1, Image2, At_33, R_cam, t_cam, 'intra');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt
