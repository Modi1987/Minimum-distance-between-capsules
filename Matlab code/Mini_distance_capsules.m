%% Efficient method for calculating the minimum distance between capsules.
% Copy right: Mohammad SAFEEA
% 16th-September-2017

function [ capsulesMiniDistance ] = Mini_distance_capsules( u,p,n,row )

%% Arreguments:
% n: scalar representing the number of capsules.
% u: (3xn) array, each column vector with index (i), represents the biginning point of
% the line segment ((at the axes)) of capsule (i).
% p: (3xn) array, each column vector with index (i), represents the end
% point of the line segment ((at the axes)) of capsule (i).
% row: (1xn) array, each element (i), represnts radious of capsule (i).

%% Return value:
% capsulesMiniDistance: (nxn) upper triangular matrix,
% while the (i,j) element of this matrix represents the minimum
% distance between capsule (i) and capsule (j).

% First calculate the sequare of the minimum distance between line segments of
% the capsules,
% The return of the following method is an (nxn) upper triangular matrix,
% while the (i,j) element of this matrix equals the sequare of the minimum
% distance between the two line segments at the axes of (capsule i) and
% (capsule j).
[seqMiniDistance]=Mini_distance_qr5_28(u,p,n);

% calculating the minimum distance between capsules
capsulesMiniDistance=zeros(n,n);
for i=1:n
    for j=i+1:n
        capsulesMiniDistance(i,j)=seqMiniDistance(i,j)^0.5-row(i)-row(j);
    end
end


end

