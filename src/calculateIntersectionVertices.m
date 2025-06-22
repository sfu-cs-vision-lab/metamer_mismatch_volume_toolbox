function vertices = calculateIntersectionVertices(M,intPoint)
%Calculates vertices of half-space intersection according to the algorithm
%described in [1].

%The rows of M contain hyper-plane coefficients
% in the form a1*x1+a2*x2+...+b<=0 i.e. [a1 a2 .... b] 
%intPoint - interior point inside half-space intersection.
%vertices - vertices of half-space intersection
%[1] Preparata & Shamos, Computational Geometry, Springer 1985

%Michal Mackiewicz, University of East Anglia, 2012-2021

if any(find(isinf(M)))||any(find(isinf(intPoint)))
    display('Non number input not allowed')
    vertices=[];
    faces=[];
    idx=[];
    return
end

%check interior point
tmp = M(:,end)+M(:,1:end-1)*intPoint;
if find(tmp>0)
    display('Interior point outside input half planes')
    vertices=[];
    faces=[];
    idx=[];
    return
end

%calculate hyper-planes in Hessian form
M = M./repmat(sqrt(sum(M(:,1:end-1).^2,2)),1,size(M,2));
%move interior point
M(:,end) = M(:,end)+M(:,1:end-1)*intPoint;
%calculate dual representation
Dual = M(:,1:end-1)./repmat(M(:,end),1,size(M,2)-1);
%convex hull in dual space
if size(Dual,2)>1
    idx = convhulln(Dual);
elseif size(Dual,2)==1
    idx = [find(Dual==max(Dual));find(Dual==min(Dual))];
end
%calculate faces in dual space
faces = calculate_faces(Dual,idx);
%the faces returned are already vertices (translated) in primal, see calculate_faces
%move back interior point
vertices = faces + repmat(intPoint',size(faces,1),1);
end