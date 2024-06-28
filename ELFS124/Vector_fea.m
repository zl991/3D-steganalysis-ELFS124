function [ F ] = Vector_fea( ori_file )
%YANG13_FEA Summary of this function goes here
%   Detailed explanation goes here
[vertex_o, face_o] = read_off(ori_file);

%%
%smoothing
vertex_o=preprocess(vertex_o);

options.symmetrize = 0;
options.normalize = 1;

% laplacian_type = 'distance';
laplacian_type = 'combinatorial';
% laplacian_type = 'conformal';

disp('--> Computing laplacian');
L = compute_mesh_laplacian(vertex_o,face_o,laplacian_type,options);

options.dt = 0.2;
options.Tmax = 0.3;
vertex_s = perform_mesh_heat_diffusion(vertex_o,face_o,L,options);
face_s=face_o;

global vring; %vring{i} is the set of vertices that are adjacent to vertex i.
% global e2f; %e2f(i,j) and e2f(j,i) are the number of the two faces adjacent to edge (i,j).
% global fring; %fring{i} is the set of faces that are adjacent to face i.

vring = compute_vertex_ring(face_o);
% e2f = compute_edge_face_ring(face_o);
% fring = compute_face_ring(face_o);

vertex_o=vertex_o';
vertex_s=vertex_s';
vertex_o_ave=vertex_o;
vertex_s_ave=vertex_s;
o_1ring_angle_mean=zeros(length(vertex_o),1);
s_1ring_angle_mean=zeros(length(vertex_s),1);
o_1ring_angle_var=zeros(length(vertex_o),1);
s_1ring_angle_var=zeros(length(vertex_s),1);
for i=1:length(vertex_o)
    vertex_o_ave(i,:)=mean(vertex_o(vring{i},:));
    vertex_s_ave(i,:)=mean(vertex_s(vring{i},:));
    
       
    o_1ring_angle =vectorAngle3d(vertex_o(i,:),vertex_o(vring{i},:))+eps;
    s_1ring_angle =vectorAngle3d(vertex_s(i,:),vertex_s(vring{i},:))+eps;
    o_1ring_angle(isnan(o_1ring_angle))=eps; 
    s_1ring_angle(isnan(s_1ring_angle))=eps;
    o_1ring_angle_mean(i,:)=mean(o_1ring_angle);
    o_1ring_angle_var(i,:)=var(o_1ring_angle);
    s_1ring_angle_mean(i,:)=mean(s_1ring_angle);
    s_1ring_angle_var(i,:)=var(s_1ring_angle);
      
end


Point_vector_angle=vectorAngle3d(vertex_o,vertex_s)+eps;

o_ave_angle =vectorAngle3d(vertex_o,vertex_o_ave);
s_ave_angle =vectorAngle3d(vertex_s,vertex_s_ave);
o_ave_angle(isnan(o_ave_angle))=eps; 
s_ave_angle(isnan(s_ave_angle))=eps;
Ave_angle_diff=abs(o_ave_angle-s_ave_angle)+eps;

Ring_angle_mean_diff=abs(o_1ring_angle_mean-s_1ring_angle_mean)+eps;
Ring_angle_var_diff=abs(o_1ring_angle_var-s_1ring_angle_var)+eps;


F=zeros(1,16);

%vertex position difference-Cartesian
F(1,1:4) = [mean(log(Point_vector_angle)),var(log(Point_vector_angle)),...
    skewness(log(Point_vector_angle)), kurtosis(log(Point_vector_angle))];
F(1,5:8) = [mean(log(Ave_angle_diff)),var(log(Ave_angle_diff)),...
    skewness(log(Ave_angle_diff)), kurtosis(log(Ave_angle_diff))];
F(1,9:12) = [mean(log(Ring_angle_mean_diff)),var(log(Ring_angle_mean_diff)),...
    skewness(log(Ring_angle_mean_diff)), kurtosis(log(Ring_angle_mean_diff))];
F(1,13:16) = [mean(log(Ring_angle_var_diff)),var(log(Ring_angle_var_diff)),...
    skewness(log(Ring_angle_var_diff)), kurtosis(log(Ring_angle_var_diff))];



% F(1,13:16) = [mean(log(edgediffN)),var(log(edgediffN)),...
%     skewness(log(edgediffN)), kurtosis(log(edgediffN))];
% F(1,17:20) = [mean(log(edgeAngle)),var(log(edgeAngle)),...
%     skewness(log(edgeAngle)), kurtosis(log(edgeAngle))];
% F(1,21:24) = [mean(log(edgeNormDiff)),var(log(edgeNormDiff)),...
%     skewness(log(edgeNormDiff)), kurtosis(log(edgeNormDiff))];
% F(1,25:36) = [mean(log(crossvector_diff)),var(log(crossvector_diff)),...
%     skewness(log(crossvector_diff)), kurtosis(log(crossvector_diff))];
% F(1,37:40) = [mean(log(crossvector_angle)),var(log(crossvector_angle)),...
%     skewness(log(crossvector_angle)), kurtosis(log(crossvector_angle))];
% F(1,41:44) = [mean(log(crossvector_normdiff)),var(log(crossvector_normdiff)),...
%     skewness(log(crossvector_normdiff)), kurtosis(log(crossvector_normdiff))];
% diff_EdgeLength=abs(EdgeLength_o-EdgeLength_s)+eps;
% 
% diff_FaceArea=abs(FaceArea_o-FaceArea_s)+eps;


%Calculate the stastics of the vertors
% F=zeros(1,8);
% 
% %vertex position difference-Cartesian
% F(1,1:4) = [mean(log(diff_EdgeLength')),var(log(diff_EdgeLength')),...
%     skewness(log(diff_EdgeLength')), kurtosis(log(diff_EdgeLength'))];
% F(1,5:8) = [mean(log(diff_FaceArea')),var(log(diff_FaceArea')),...
%     skewness(log(diff_FaceArea')), kurtosis(log(diff_FaceArea'))];
% F(1,9:12) = [mean(log(diff_rho')),var(log(diff_rho')),...
%     skewness(log(diff_rho')), kurtosis(log(diff_rho'))];
% F(1,13:16) = [mean(log(diff_BSang')),var(log(diff_BSang')),...
%     skewness(log(diff_BSang')), kurtosis(log(diff_BSang'))];
% F(1,17:20) = [mean(log(diff_EAaz')),var(log(diff_EAaz')),...
%     skewness(log(diff_EAaz')), kurtosis(log(diff_EAaz'))];
% F(1,21:24) = [mean(log(diff_EAel')),var(log(diff_EAel')),...
%     skewness(log(diff_EAel')), kurtosis(log(diff_EAel'))];
% F(1,25:28) = [mean(log(diff_EABSang')),var(log(diff_EABSang')),...
%     skewness(log(diff_EABSang')), kurtosis(log(diff_EABSang'))];
% F(1,29:32) = [mean(log(diff_Erho')),var(log(diff_Erho')),...
%     skewness(log(diff_Erho')), kurtosis(log(diff_Erho'))];
end

