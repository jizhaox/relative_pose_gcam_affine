function match_info = match_type_gcam_ac(match_type)

%% Two types of affine correspondence 
%% inter-cam affine correspondences contain
% (1) 1 affine correspondence: (cam1 at time i) <--> (cam2 at time j)
% (2) 1 affine correspondence: (cam2 at time i) <--> (cam1 at time j)
%% intra-cam affine correspondences contain
% (1) 1 affine correspondence: (cam1 at time i) <--> (cam1 at time j)
% (2) 1 affine correspondence: (cam2 at time i) <--> (cam2 at time j)

if strcmp(match_type, 'inter')
    match_info{1} = struct('idx1', 1, 'idx2', 2);
    match_info{2} = struct('idx1', 2, 'idx2', 1);
elseif strcmp(match_type, 'intra')
    match_info{1} = struct('idx1', 1, 'idx2', 1);
    match_info{2} = struct('idx1', 2, 'idx2', 2);
else
    error('unsupported match type!')
end