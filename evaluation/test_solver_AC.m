clear;
addpath('../mex')

%% Two types of affine correspondence
%% inter-cam affine correspondences refer to correspondences which are seen by different cameras over two consecutive views, which contain
% (1) 1 affine correspondence: (cam1 at time i) <--> (cam2 at time j)
% (2) 1 affine correspondence: (cam2 at time i) <--> (cam1 at time j)
%% intra-cam affine correspondences refer to correspondences which are seen by the same camera over two consecutive views, which contain
% (1) 1 affine correspondence: (cam1 at time i) <--> (cam1 at time j)
% (2) 1 affine correspondence: (cam2 at time i) <--> (cam2 at time j)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1AC plane solver & inter-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------1AC plane solver & inter-camera case--------------')
%% generate synthetic data
[Image1, Image2, At, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('3DOF', 'inter');

%% check data
% For noise-free data, the residual of epipolar constraint and affine transformation constraints should be small (Eqs.(4) and (9) in the paper)
err_epipolar = zeros(2, 1);
err_affinetransformation = zeros(2, 2);
P1All = zeros(3, 2);
P2All = zeros(3, 2);
Line_iAll = zeros(6, 2);
Line_jAll = zeros(6, 2);
for ii = 1:2
    if ii <= 1
        idx1 = 1;
        idx2 = 2;
    else
        idx1 = 2;
        idx2 = 1;
    end
    x1 = Image1(:, ii);
    x2 = Image2(:, ii);
    
    P1All(:, ii) = R_cam (:,:,idx1)*x1;
    P2All(:, ii) = R_cam (:,:,idx2)*x2;
    
    u1 = P1All(:, ii)/norm(P1All(:, ii));
    u2 = P2All(:, ii)/norm(P2All(:, ii));
    
    Line_iAll(:, ii) = [u1; cross(t_cam(:,idx1),u1)];
    Line_jAll(:, ii) = [u2; cross(t_cam(:,idx2),u2)];
    
    % Generalized Epipolar Constraint: Eq.(4) in the paper
    E_gem = zeros(6,6);
    E_gem(1:3,1:3) = skew(t_gt)*R_gt;
    E_gem(1:3,4:6) = R_gt;
    E_gem(4:6,1:3) = R_gt;
    err_epipolar(ii) = Line_jAll(:, ii)'*E_gem*Line_iAll(:, ii);
    
    % Affine Transformation Constraint: Eq.(9) in the paper
    At_33(1:2,1:2,ii) = At(:,:,ii);
    At_33(3,3,ii) = 0;
    equationerror = R_cam (:,:,idx1)'*(skew(t_cam(:,idx1))*R_gt'+ R_gt'*skew(t_gt) - R_gt'*skew(t_cam(:,idx2)))*P2All(:, ii) - At_33(:,:,ii)'*R_cam (:,:,idx2)'*(R_gt*skew(t_cam(:,idx1))+ skew(t_gt)*R_gt - skew(t_cam(:,idx2))*R_gt)*P1All(:, ii);
    err_affinetransformation(:,ii) = equationerror(1:2,1);
end
disp(['The maximum residual of epipolar constraint: ' num2str(max(abs(err_epipolar(:))))])
disp(['The maximum residual of affine transformation constraints: ' num2str(max(abs(err_affinetransformation(:))))])

%% run inter-1AC plane solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose3d_1AC_plane(Image1(:,1), Image2(:,1), At_33(:,:,1), R_cam, t_cam, 'inter');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC plane solver & inter-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC plane solver & inter-camera case--------------')
%% generate synthetic data
[Image1, Image2, At, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('3DOF', 'inter');

%% check data
% For noise-free data, the residual of epipolar constraint and affine transformation constraints should be small (Eqs.(4) and (9) in the paper)
err_epipolar = zeros(2, 1);
err_affinetransformation = zeros(2, 2);
P1All = zeros(3, 2);
P2All = zeros(3, 2);
Line_iAll = zeros(6, 2);
Line_jAll = zeros(6, 2);
for ii = 1:2
    if ii <= 1
        idx1 = 1;
        idx2 = 2;
    else
        idx1 = 2;
        idx2 = 1;
    end
    x1 = Image1(:, ii);
    x2 = Image2(:, ii);
    
    P1All(:, ii) = R_cam (:,:,idx1)*x1;
    P2All(:, ii) = R_cam (:,:,idx2)*x2;
    
    u1 = P1All(:, ii)/norm(P1All(:, ii));
    u2 = P2All(:, ii)/norm(P2All(:, ii));
    
    Line_iAll(:, ii) = [u1; cross(t_cam(:,idx1),u1)];
    Line_jAll(:, ii) = [u2; cross(t_cam(:,idx2),u2)];
    
    % Generalized Epipolar Constraint: Eq.(4) in the paper
    E_gem = zeros(6,6);
    E_gem(1:3,1:3) = skew(t_gt)*R_gt;
    E_gem(1:3,4:6) = R_gt;
    E_gem(4:6,1:3) = R_gt;
    err_epipolar(ii) = Line_jAll(:, ii)'*E_gem*Line_iAll(:, ii);
    
    % Affine Transformation Constraint: Eq.(9) in the paper
    At_33(1:2,1:2,ii) = At(:,:,ii);
    At_33(3,3,ii) = 0;
    equationerror = R_cam (:,:,idx1)'*(skew(t_cam(:,idx1))*R_gt'+ R_gt'*skew(t_gt) - R_gt'*skew(t_cam(:,idx2)))*P2All(:, ii) - At_33(:,:,ii)'*R_cam (:,:,idx2)'*(R_gt*skew(t_cam(:,idx1))+ skew(t_gt)*R_gt - skew(t_cam(:,idx2))*R_gt)*P1All(:, ii);
    err_affinetransformation(:,ii) = equationerror(1:2,1);
end
disp(['The maximum residual of epipolar constraint: ' num2str(max(abs(err_epipolar(:))))])
disp(['The maximum residual of affine transformation constraints: ' num2str(max(abs(err_affinetransformation(:))))])

%% run inter-2AC plane solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose3d_2AC_plane(Image1, Image2, At_33, R_cam, t_cam, 'inter');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC plane solver & intra-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC plane solver & intra-camera case--------------')
%% generate synthetic data
[Image1, Image2, At, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('3DOF', 'intra');

%% check data
% For noise-free data, the residual of epipolar constraint and affine transformation constraints should be small (Eqs.(4) and (9) in the paper)
err_epipolar = zeros(2, 1);
err_affinetransformation = zeros(2, 2);
P1All = zeros(3, 2);
P2All = zeros(3, 2);
Line_iAll = zeros(6, 2);
Line_jAll = zeros(6, 2);
for ii = 1:2   
    idx1 = ii;
    idx2 = ii;
    
    x1 = Image1(:, ii);
    x2 = Image2(:, ii);
    
    P1All(:, ii) = R_cam (:,:,idx1)*x1;
    P2All(:, ii) = R_cam (:,:,idx2)*x2;
    
    u1 = P1All(:, ii)/norm(P1All(:, ii));
    u2 = P2All(:, ii)/norm(P2All(:, ii));
    
    Line_iAll(:, ii) = [u1; cross(t_cam(:,idx1),u1)];
    Line_jAll(:, ii) = [u2; cross(t_cam(:,idx2),u2)];
    
    % Generalized Epipolar Constraint: Eq.(4) in the paper
    E_gem = zeros(6,6);
    E_gem(1:3,1:3) = skew(t_gt)*R_gt;
    E_gem(1:3,4:6) = R_gt;
    E_gem(4:6,1:3) = R_gt;
    err_epipolar(ii) = Line_jAll(:, ii)'*E_gem*Line_iAll(:, ii);
    
    % Affine Transformation Constraint: Eq.(9) in the paper
    At_33(1:2,1:2,ii) = At(:,:,ii);
    At_33(3,3,ii) = 0;
    AtempAll(:,:,ii) =  At_33(:,:,ii)'*R_cam (:,:,ii)';
    equationerror = R_cam (:,:,idx1)'*(skew(t_cam(:,idx1))*R_gt'+ R_gt'*skew(t_gt) - R_gt'*skew(t_cam(:,idx2)))*P2All(:, ii) - At_33(:,:,ii)'*R_cam (:,:,idx2)'*(R_gt*skew(t_cam(:,idx1))+ skew(t_gt)*R_gt - skew(t_cam(:,idx2))*R_gt)*P1All(:, ii);
    err_affinetransformation(:,ii) = equationerror(1:2,1);
end
disp(['The maximum residual of epipolar constraint: ' num2str(max(abs(err_epipolar(:))))])
disp(['The maximum residual of affine transformation constraints: ' num2str(max(abs(err_affinetransformation(:))))])

%% run intra-2AC plane solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose3d_2AC_plane(Image1, Image2, At_33, R_cam, t_cam, 'intra');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC vertical solver & inter-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC vertical solver & inter-camera case--------------')
%% generate synthetic data
[Image1, Image2, At, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('4DOF', 'inter');

%% check data
% For noise-free data, the residual of epipolar constraint and affine transformation constraints should be small (Eqs.(4) and (9) in the paper)
err_epipolar = zeros(2, 1);
err_affinetransformation = zeros(2, 2);
P1All = zeros(3, 2);
P2All = zeros(3, 2);
Line_iAll = zeros(6, 2);
Line_jAll = zeros(6, 2);
for ii = 1:2
    if ii <= 1
        idx1 = 1;
        idx2 = 2;
    else
        idx1 = 2;
        idx2 = 1;
    end
    
    x1 = Image1(:, ii);
    x2 = Image2(:, ii);

    P1All(:, ii) = R_cam (:,:,idx1)*x1;
    P2All(:, ii) = R_cam (:,:,idx2)*x2;

    u1 = P1All(:, ii)/norm(P1All(:, ii));
    u2 = P2All(:, ii)/norm(P2All(:, ii));

    Line_iAll(:, ii) = [u1; cross(t_cam(:,idx1),u1)];
    Line_jAll(:, ii) = [u2; cross(t_cam(:,idx2),u2)];

    % Generalized Epipolar Constraint: Eq.(4) in the paper
    E_gem = zeros(6,6);
    E_gem(1:3,1:3) = skew(t_gt)*R_gt;
    E_gem(1:3,4:6) = R_gt;
    E_gem(4:6,1:3) = R_gt;
    err_epipolar(ii) = Line_jAll(:, ii)'*E_gem*Line_iAll(:, ii);

    % Affine Transformation Constraint: Eq.(9) in the paper
    At_33(1:2,1:2,ii) = At(:,:,ii);
    At_33(3,3,ii) = 0;
    AtempAll(:,:,ii) =  At_33(:,:,ii)'*R_cam (:,:,ii)';
    equationerror = R_cam (:,:,idx1)'*(skew(t_cam(:,idx1))*R_gt'+ R_gt'*skew(t_gt) - R_gt'*skew(t_cam(:,idx2)))*P2All(:, ii) - At_33(:,:,ii)'*R_cam (:,:,idx2)'*(R_gt*skew(t_cam(:,idx1))+ skew(t_gt)*R_gt - skew(t_cam(:,idx2))*R_gt)*P1All(:, ii);
    err_affinetransformation(:,ii) = equationerror(1:2,1);
end
disp(['The maximum residual of epipolar constraint: ' num2str(max(abs(err_epipolar(:))))])
disp(['The maximum residual of affine transformation constraints: ' num2str(max(abs(err_affinetransformation(:))))])

%% run 2AC vertical solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose4d( Image1, Image2, At_33, R_cam, t_cam, 'inter');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2AC vertical solver & intra-camera case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------2AC vertical solver & intra-camera case--------------')
%% generate synthetic data
[Image1, Image2, At, R_cam , t_cam, R_gt, t_gt] = generate_2AC_data('4DOF', 'intra');

%% check data
% For noise-free data, the residual of epipolar constraint and affine transformation constraints should be small (Eqs.(4) and (9) in the paper)
err_epipolar = zeros(2, 1);
err_affinetransformation = zeros(2, 2);
P1All = zeros(3, 2);
P2All = zeros(3, 2);
Line_iAll = zeros(6, 2);
Line_jAll = zeros(6, 2);
for ii = 1:2
    idx1 = ii;
    idx2 = ii;
    
    x1 = Image1(:, ii);
    x2 = Image2(:, ii);

    P1All(:, ii) = R_cam (:,:,idx1)*x1;
    P2All(:, ii) = R_cam (:,:,idx2)*x2;

    u1 = P1All(:, ii)/norm(P1All(:, ii));
    u2 = P2All(:, ii)/norm(P2All(:, ii));

    Line_iAll(:, ii) = [u1; cross(t_cam(:,idx1),u1)];
    Line_jAll(:, ii) = [u2; cross(t_cam(:,idx2),u2)];

    % Generalized Epipolar Constraint: Eq.(4) in the paper
    E_gem = zeros(6,6);
    E_gem(1:3,1:3) = skew(t_gt)*R_gt;
    E_gem(1:3,4:6) = R_gt;
    E_gem(4:6,1:3) = R_gt;
    err_epipolar(ii) = Line_jAll(:, ii)'*E_gem*Line_iAll(:, ii);

    % Affine Transformation Constraint: Eq.(9) in the paper
    At_33(1:2,1:2,ii) = At(:,:,ii);
    At_33(3,3,ii) = 0;
    AtempAll(:,:,ii) =  At_33(:,:,ii)'*R_cam (:,:,ii)';
    equationerror = R_cam (:,:,idx1)'*(skew(t_cam(:,idx1))*R_gt'+ R_gt'*skew(t_gt) - R_gt'*skew(t_cam(:,idx2)))*P2All(:, ii) - At_33(:,:,ii)'*R_cam (:,:,idx2)'*(R_gt*skew(t_cam(:,idx1))+ skew(t_gt)*R_gt - skew(t_cam(:,idx2))*R_gt)*P1All(:, ii);
    err_affinetransformation(:,ii) = equationerror(1:2,1);
end
disp(['The maximum residual of epipolar constraint: ' num2str(max(abs(err_epipolar(:))))])
disp(['The maximum residual of affine transformation constraints: ' num2str(max(abs(err_affinetransformation(:))))])

%% run 2AC vertical solver
[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose4d( Image1, Image2, At_33, R_cam, t_cam, 'intra');

Rf1tof2_recover - R_gt
Tf1tof2_recover - t_gt
