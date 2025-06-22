function faces  = calculate_faces(vertices,idx)
%Calculates hyper-plane coefficients (faces) passing through a set of
%vertices. 
% idx - rows in idx contain indices of vertices forming a face.

%Matrix faces does not contain offsets in its last column as offsets equal
%1 i.e. faces are already ready for calculating their duals. When written
%in this form they are already vertices in their dual space (original
%primal)i.e. vertices of the required half-space intersection. These
%vertices only need to be translated by the interior point. 

%Michal Mackiewicz, University of East Anglia, 2012-2021

dim = numel(idx(1,:));
faces = zeros(size(idx,1),dim);
discard=false(1,size(idx,1));
for i=1:size(idx,1)
    if rank(vertices(idx(i,:),:))<dim
        discard(i)=1;
        continue;
    end
    faces(i,:) = (vertices(idx(i,:),:)\-ones(dim,1))';  
end
faces(discard,:)=[];
