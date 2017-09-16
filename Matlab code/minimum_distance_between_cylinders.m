%% Efficient method for calculating the minimum distance between capsules.
% Copy right: Mohammad SAFEEA
% 16th-September-2017

n=3000; % number of capsules

u=rand(3,n)*100; % biginning point of line segment at axes of capsule
p=rand(3,n)*100; % end point of line segment at axes of capsule
row=rand(1,n)*0.5; % radious of capsules 

% The following method calculates the minimum distance muttually between
% the capsules
[ capsulesMiniDistance ] = Mini_distance_capsules( u,p,n,row );
% Where: capsulesMiniDistance: (nxn) upper triangular matrix,
% while the (i,j) element of this matrix represents the minimum
% distance between capsule (i) and capsule (j).

% If the value of element (i,j) is negative then the two cylinders (i,j)
% are in collision state