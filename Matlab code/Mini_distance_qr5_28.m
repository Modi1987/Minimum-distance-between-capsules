%% Efficient method for calculating the minimum distance between line segments.
% Copy right: Mohammad SAFEEA
% 16th-September-2017

function [mini_distance_seq]=Mini_distance_qr5_28(u,p,num_of_cyl)
%% This function calculates all the sequare of minimum distance
%% between the line segments defined by starting points (p) and 
%% end points (u),

%% Arreguments
% (u) is 3xn matrix, each column of which corresponds to the end of the segment
% (p) is 3xn matrix, each column of which corresponds to the beginning point of the segment
% (num_of_cyl) is scalar representing the total number of line segments

%% Return value
% (mini_distance_seq) is nXn upper triangular matrix, each (i,j) entry of this matrix
% represents the sequare of the minimum distance between segments i and j.

mini_distance_seq=zeros(num_of_cyl,num_of_cyl);
S1= u(:,:)-p(:,:);
riiseq=sum(S1.^2);
rii=riiseq.^0.5;

q1(1,1:num_of_cyl)=S1(1,:)./rii;
q1(2,1:num_of_cyl)=S1(2,:)./rii;
q1(3,1:num_of_cyl)=S1(3,:)./rii;

for i_count=1:num_of_cyl
for j_count=i_count+1:num_of_cyl
        r1j=q1(1,i_count)*S1(1,j_count)+q1(2,i_count)*S1(2,j_count)+...
            q1(3,i_count)*S1(3,j_count); % 5 operations
       
        y=p(:,i_count)-p(:,j_count);% 3 operations

        node14=q1(1,i_count)*y(1)+...
            q1(2,i_count)*y(2)+q1(3,i_count)*y(3); % 5 operations
        
        node24r2j=(node14*r1j-S1(1,j_count)*y(1)-S1(2,j_count)*y(2)-...
            S1(3,j_count)*y(3)); % 7 operations       
        node13=node14+rii(i_count); % 1 operation
        node14seq=node14*node14; % 1 operation
        yTy=y(1)*y(1)+y(2)*y(2)+y(3)*y(3); % 5 operation
        % Up to here total of 27 obligatory O(n2) mathematical operations, 

%% The nodes of the paralello are as the following:
%                             1----------------------------2
%
%
%                4----------------------------3 
%
%         r4---------------------------r3===============X axes=============================>


%% The closest point is on the horizontal
if (  (node24r2j>0)&&(node14*node13<0) )% 1 operation
                % 7 operations
                mini_distance_seq(i_count,j_count)=yTy-node14seq;% 1 operation
                continue;
end

        r2jseq=(riiseq(j_count)-r1j*r1j);%2 operations
        node11=node14-r1j; % 1 operaiton
        node12=node11+rii(i_count); % 1 operation 
        node21r2j=r2jseq+node24r2j; % 1 operation
        
%% The two segments are parallel       
if(r2jseq==0) 
%% Up to here 33 mandatory opreations
    % in such case the parallelopipe degenerate into Four coincedent
    % segments, two cases (a) and (b):
    
    % CASE (a):
    % 4-----1-----3-----2
    % 4-----3-----1-----2
    if ( node14*node12<=0)
                mini_distance_seq(i_count,j_count)=yTy-node14seq;
                continue;
    end
    
    % CASE (b):
    % 1-----2-----4-----3
    % 1-----4-----2-----3
    if ( node11*node13<=0)
                mini_distance_seq(i_count,j_count)=yTy-node14seq;
                continue;
    end
    
%% Otherwise the closest points are on the extremeties of the segments
    if(r1j<0)
    %% r1j<0
        if(node14>0)
            mini_distance_seq(i_count,j_count)=yTy;         
        elseif(node12<0)
            mini_distance_seq(i_count,j_count)=node12*node12+yTy-node14seq;  
        end
    else
	%% r1j>0
        if(node11>0)
            mini_distance_seq(i_count,j_count)=node11*node11+yTy-node14seq;      
        elseif(node13<0)
            mini_distance_seq(i_count,j_count)=node13*node13+yTy-node14seq; 
        end        
    end    
end
%%%%%%%%%%%%%%%%%%%%%%
%% The two cylinders are not parallel%
%%%%%%%%%%%%%%%%%%%%%%     

if ( (node21r2j<0) &&(node11*node12<0))
                mini_distance_seq(i_count,j_count)=(node21r2j*node21r2j-node24r2j*node24r2j)/r2jseq+...
                    yTy-node14seq;
                continue;
end


% The two cylinders are not parallel
%                       rii       r1j  
% [q1  q2 ]                        
%                       0        r2j 
r4=node11+r1j*node21r2j/r2jseq;
r3=r4+rii(i_count);

%  If the origin is inside the parallelo
if ((node21r2j*node24r2j<0) && (r3*r4<0)) 
                mini_distance_seq(i_count,j_count)=yTy-...
                    node14seq-node24r2j*node24r2j/r2jseq;

elseif ((r4*r4<r3*r3))
 %% closest oblique is the segment (1,4)
        if ( r1j*node11>node21r2j )% then point one is the closest
            node12seq=node12*node12;
            node11seq=node11*node11;
            if ( node11seq<node12seq)
                mini_distance_seq(i_count,j_count)=node11seq+...
                    +yTy...
                    -node14seq...
                    +(node21r2j*node21r2j-node24r2j*node24r2j)/r2jseq;  
            else
                mini_distance_seq(i_count,j_count)=node12seq+...
                    yTy-node14seq...
                    +(node21r2j*node21r2j-node24r2j*node24r2j)/r2jseq; 
            end
        elseif ( r1j*node14<node24r2j)
            if ( node13*node13>node14seq)
                mini_distance_seq(i_count,j_count)=yTy;
            else
                mini_distance_seq(i_count,j_count)=node13*node13+...
                    yTy-node14seq;
            end
        else
            mini_distance_seq(i_count,j_count)=r4*r4*r2jseq/riiseq(j_count)+...
                yTy-...
                node14seq-node24r2j*node24r2j/r2jseq;
        end
else        
%% closest oblique is segment (2,3)     
        if ( r1j*node12>node21r2j)
            node11seq=node11*node11;
            node12seq=node12*node12;
            if ( node11seq<node12seq)
                mini_distance_seq(i_count,j_count)=node11seq+...
                    +yTy-node14seq...
                    +(node21r2j*node21r2j-node24r2j*node24r2j)/r2jseq;  
            else
                mini_distance_seq(i_count,j_count)=node12seq+...
                    yTy...
                    -node14seq...
                    +(node21r2j*node21r2j-node24r2j*node24r2j)/r2jseq;
            end
        elseif ( r1j*node13<node24r2j)
            if ( node13*node13<node14seq)
                mini_distance_seq(i_count,j_count)=node13*node13+...
                    yTy-node14seq;
            else
                mini_distance_seq(i_count,j_count)=yTy;
            end
        else
            mini_distance_seq(i_count,j_count)=r3*r3*r2jseq/riiseq(j_count)+...
                yTy-node14seq...
                -node24r2j*node24r2j/r2jseq;
        end
end
end %% for first loop
end %% for second loop



