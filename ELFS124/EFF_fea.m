function [ F ] = EFF_fea( ori_file ,dt,Tmax)

[vertex_o, face_o] = read_off(ori_file);

%%
%smoothing
options.symmetrize = 0;
options.normalize = 1;

% laplacian_type = 'distance';
laplacian_type = 'combinatorial';
% laplacian_type = 'conformal';

disp('--> Computing laplacian');
L = compute_mesh_laplacian(vertex_o,face_o,laplacian_type,options);

options.dt = dt;
options.Tmax = Tmax;
vertex_s = perform_mesh_heat_diffusion(vertex_o,face_o,L,options);
face_s=face_o;

vertex_o=preprocess(vertex_o);
vertex_s=preprocess(vertex_s);

% global vring; %vring{i} is the set of vertices that are adjacent to vertex i.
global e2f; %e2f(i,j) and e2f(j,i) are the number of the two faces adjacent to edge (i,j).
% global fring; %fring{i} is the set of faces that are adjacent to face i.

% vring = compute_vertex_ring(face_o);
e2f = compute_edge_face_ring(face_o);
% fring = compute_face_ring(face_o);

edge_o=meshEdges(face_o');

edge_s=edge_o;

crossvector_o=zeros(3,length(edge_o));
crossvector_s=zeros(3,length(edge_o));
for i=1:length(edge_o)
    %get the face index adjacent to the edge_i
    face_1=e2f(edge_o(i,1),edge_o(i,2));
    face_2=e2f(edge_o(i,2),edge_o(i,1));
    %get the vertex index in the two adjacent faces
    v_f1=my_setdiff(face_o(:,face_1),edge_o(i,:));
    v_f2=my_setdiff(face_o(:,face_2),edge_o(i,:));
    %vector form v_f1 to v_f2
    crossvector_o(:,i)=vertex_o(:,v_f2)-vertex_o(:,v_f1);
    crossvector_s(:,i)=vertex_s(:,v_f2)-vertex_s(:,v_f1);
end

edge_vector_o=vertex_o(:,(edge_o(:,1)))-vertex_o(:,(edge_o(:,2)));
edge_vector_s=vertex_s(:,(edge_s(:,1)))-vertex_s(:,(edge_s(:,2)));

edgediff=abs(edge_vector_o-edge_vector_s)+eps;
edgediffN=vectorNorm3d(edgediff');
edgeAngle=vectorAngle3d(edge_vector_o',edge_vector_s')+eps;
edgeNormDiff=abs(vectorNorm3d(edge_vector_o')-vectorNorm3d(edge_vector_s'))+eps;

crossvector_diff=abs(crossvector_o'-crossvector_s')+eps;
crossvector_diffN=vectorNorm3d(crossvector_diff');
crossvector_angle=vectorAngle3d(crossvector_o',crossvector_s')+eps;
crossvector_normdiff=abs(vectorNorm3d(crossvector_o')-vectorNorm3d(crossvector_s'))+eps;


F=zeros(1,48);

%vertex position difference-Cartesian
F(1,1:12) = [mean(log(edgediff')),var(log(edgediff')),...
    skewness(log(edgediff')), kurtosis(log(edgediff'))];
F(1,13:16) = [mean(log(edgediffN)),var(log(edgediffN)),...
    skewness(log(edgediffN)), kurtosis(log(edgediffN))];
F(1,17:20) = [mean(log(edgeAngle)),var(log(edgeAngle)),...
    skewness(log(edgeAngle)), kurtosis(log(edgeAngle))];
F(1,21:24) = [mean(log(edgeNormDiff)),var(log(edgeNormDiff)),...
    skewness(log(edgeNormDiff)), kurtosis(log(edgeNormDiff))];
F(1,25:36) = [mean(log(crossvector_diff)),var(log(crossvector_diff)),...
    skewness(log(crossvector_diff)), kurtosis(log(crossvector_diff))];
F(1,37:40) = [mean(log(crossvector_angle)),var(log(crossvector_angle)),...
    skewness(log(crossvector_angle)), kurtosis(log(crossvector_angle))];
F(1,41:44) = [mean(log(crossvector_normdiff)),var(log(crossvector_normdiff)),...
    skewness(log(crossvector_normdiff)), kurtosis(log(crossvector_normdiff))];
F(1,45:48) = [mean(log(crossvector_diffN)),var(log(crossvector_diffN)),...
    skewness(log(crossvector_diffN)), kurtosis(log(crossvector_diffN))];


end

