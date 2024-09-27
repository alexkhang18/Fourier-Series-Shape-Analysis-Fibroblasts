function [orientation,coherence] = aSMA_orientation(I,sigma)
%aSMA_orientation - computes local orientation of aSMA image
%   This function computes the local orientation within an aSMA image by
%   computting the gradients in the x- and y- directions. A Guassian blurr
%   is then added to the gradients and the structure tensor is assembled.
%   The eigenvalues of the structure tensor is used to determine coherence
%   which is a measure of anisotropy ranging from [0,1]. The eigenvectors
%   of the structure tensor are used to determine the local fiber orientation.
%   The first eigenvector points along the direction of greatest change
%   (i.e., the largest gradient). The second eigvector points along the
%   local fiber direction. This analysis is done on a pixel-by-pixel basis.

h = 1/10*ones(10,1);
H = h*h';

% I = imgaussfilt(I);  % Apply median filter with a 3x3 kernel

% image gradient
[Gx,Gy] = imgradientxy(I);

% figure
% imshow(Gx)
% 
% figure
% imshow(Gy)


% filtered gradients
w_Gx = imgaussfilt(Gx.^2,sigma);
w_Gy = imgaussfilt(Gy.^2,sigma);
w_Gxy = imgaussfilt(Gx.*Gy,sigma);
% w_Gx = Gx.^2;
% w_Gy = Gy.^2;
% w_Gxy = Gx.*Gy;


% imwrite(Gx,"Gx.png")
% imwrite(Gy,"Gy.png")
% 
% imshow(w_Gx,[])
% exportgraphics(gcf,'Gx.png','Resolution',300)
% imshow(w_Gy,[])
% exportgraphics(gcf,'Gy.png','Resolution',300)
% imshow(w_Gxy,[])
% exportgraphics(gcf,'Gxy.png','Resolution',300)
% 
% imwrite(w_Gx,"w_Gx.png")
% imwrite(w_Gy,"w_Gy.png")
% imwrite(w_Gxy,"w_Gxy.png")


% pre-allocation for coherence and orientation
coherence = zeros(size(I));
orientation = zeros(size(I));

% loop solves for coherence and orientation on a pixel-by-pixel basis
for i = 1:length(coherence(:))

    % structure tensor 
    J = [w_Gx(i),w_Gxy(i);w_Gxy(i),w_Gy(i)];

    % eigendecomposition 
    [V,D] = eig(J);
    Eig = diag(D) ;
    [~,idx] = sort(Eig,'descend'); % sorts eigenvalues
    PV = Eig(idx); % principal values
    PC = V(:,idx); % principal directions
    coherence(i) = (PV(1) - PV(2))/(PV(1) + PV(2));
    orientation(i) = atan2d(PC(2,2),PC(2,1));

end

% rescales orientation range to [-90,90]
indx = find(orientation>90 & orientation<=180);
indx2 = find(orientation>=-180 & orientation<-90);
orientation(indx) = orientation(indx) - 180;
orientation(indx2) = orientation(indx2) + 180;

end