function [k,b]=objectColSol_sphericalSampling(sensors,num,rngset,orth_flag)
% returns half-space representation of 6-D Object Colour Solid
% sensors - N_Wavelengths x 6
% num - number of normals
% rngset - used to have repeateable results from rnd function, can be
% useful in some experiments.
% orth_flag is 0 or 1, 1 for "orthonormal sensors" - see Fig. 2 in [1]
%k - surface normals
%b - offsets

%[1] Michal Mackiewicz, Hans Jakob Rivertz, and Graham Finlayson, 
%  "Spherical sampling methods for the calculation of metamer mismatch volumes," 
%  J. Opt. Soc. Am. A 36, 96-104 (2019)

%Michal Mackiewicz, Graham Finlayson, University of East Anglia, 2016-2021

function r=sphsampn(N,nm,rngset)
    rng(rngset)
    r=randn(nm,N);
    r=r./sqrt(sum(r'.^2)'*ones(1,N));
end

[~, b] = size(sensors);
[u s v] = svd(sensors);
ort=u(:,1:6);
s = diag(diag(s));
 
k=sphsampn(b,num,rngset);

if orth_flag
    k = k*inv(s)*v';
    k = k./(sqrt(sum(k.^2,2))*ones(1,6));
end

%for a very large number of samples, using batches helps with memory
%problems.

responses=zeros(num,6);
batchs = 100;
numb = num/batchs;

for batch = 1:numb
        spec=sensors*k(batchs*(batch-1)+1:batchs*batch,:)';
        spec=spec>0;
        responses(batchs*(batch-1)+1:batchs*batch,:)=spec'*sensors;
end

b = sum(k.*responses,2);
%keyboard
end