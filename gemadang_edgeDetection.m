%Canny edge detection algorithm
% Geri Madanguit, 02/24/17


% READ ME: In this code, the outputs of an image are as followed in this
% specific order: Original Image, Smoothed Imaged through Gaussian Filter,
% the Edge Strength of the Image, 1-Pixel Wide Image, and the Image with
% Hysteresis Thresh.
% To run through the code for each different image, you need to comment all
% the initializations of A except for the one initialization of the image
% you want to display. 

clear all;
clc;

flowers = imread('Flowers.jpg');
syr1 = imread('Syracuse_01.jpg');
syr2 = imread('Syracuse_02.jpg');
mypic = imread('fun.jpg');

%COMMENT ALL EXCEPT FOR ONE BELOW TO RUN CODE ON SPECIFIC IMAGE

%A = flowers;
%A = syr1;
%A = syr2;
A = syr1;

figure, imshow(A);


%1) noise smoothing 
%   supress as much noise as possible, without destroying true edges
%   assume Gaussian noise



%%GAUSSIAN SMOOTHING
kerG2r = 3; %gaussian mask row
kerG2c = 3; %gaussian mask column
std =1.5; %sigma
    
    centerx = round((kerG2r+1)/2);
    centery = round((kerG2c+1)/2);
    
    int_gauss_ker = ones(kerG2r,kerG2c); 
    float_ker = zeros(kerG2r,kerG2c);

    for k=1:kerG2c
        for h=1:kerG2r
            float_ker(h,k) = exp((-(((h-centerx)^2+(k-centery)^2)/2*(std)^2)));
        end
    end
    
    gmin = float_ker(centerx,centery);
    norm_factor = 1/gmin;
    int_gauss_ker = round(norm_factor*float_ker);
    gauss_filtered_im = conv2(int_gauss_ker,A,'full');

    max_gauss_filtered = max(max(gauss_filtered_im));

    gauss_filtered_im = uint8((gauss_filtered_im/max_gauss_filtered)*255);

    figure, imshow(gauss_filtered_im);

    
    
% CANNY_ENHANCER

% Edge normal: The direction (unit vector) of the maximum intensity 
%   variation at the edge point. 
%   The direction that is perpendicular to the edge.
% Edge direction: the direction that is perpendicular to the edge normal, 
%   and thus tangent to the contour of which the edge is part.
% Compute the gradient components, Jx and Jy.

jx_vector = [-1,0,1];
jy_vector = jx_vector';

[row, col] = size(gauss_filtered_im);

for r=1:row
    jx(r,:) = conv(jx_vector, gauss_filtered_im(r,:));
end

for c=1:col
    jy(:,c) = conv(jy_vector, gauss_filtered_im(:,c));
end

% Estimate the edge strength, es(i,j).
% Estimate the orientation of the edge normal, eo(i,j).

for r=1:row
    for c=1:col
        edge_s(r,c) = sqrt((jx(r,c)^2)+(jy(r,c)^2));
        edge_o(r,c) = radtodeg(atan(jy(r,c)/jx(r,c)));
    end
end

% The output is a strength image, Es , and an orientation image Eo.
%figure, imshow(edge_s);
%figure, imshow(edge_o);



% NONMAX SUPPRESSION
%Find the direction d?k (out of 0, 45, 90 and 135 degrees), w
%   which best approximates the direction Eo(i,j).

%2) edge enhancement 
%   design filter responding to edges
%   Filter output: large at edge pixels, low elsewhere
%   edges can be located at local maxima

%modeling edges: - Step edges - Ridge edges - Roof edges

[edgeR, edgeC] = size(edge_o);

for r=1:edgeR
    for c=1:edgeC
        num = edge_o(r,c);
        if (num < 0)
            num = num + 180;
        end

        if ((num >= 0 && num <= 22.5) || (num >= 157.5 && num <= 180))   
            edge_d(r,c) =  0;
        elseif (num >= 22.5 && num <= 67.5)
           edge_d(r,c) =  45;
        elseif (num >= 67.5 && num <= 112.5)
           edge_d(r,c) =  90;
        elseif (num >= 112.5 && num <= 157.5)
           edge_d(r,c) =  135;
        else
            edge_d(r,c) = NaN;
        end
    end    
end

%figure, imshow(edge_d);
edge_s = uint8(edge_s);
figure, imshow(edge_s);

[dkR, dkC] = size(edge_d);

for r=1:dkR
    for c=1:dkC
        
        one = 0;zero = 0;
        
        top_left = 0;top_right = 0;top_edge = 0;
        bottom_left = 0;bottom_right = 0;bottom_edge = 0;
        left_edge = 0;right_edge = 0;
        
        num = edge_d(r,c);
        
        if (r==1)   %topEdge
            top_edge = 1;
            if(c==1)    %top_leftcorner
                top_left = 1;
                top_edge = 0;
            elseif(c==dkC) %top_rightcorner
                top_right = 1;
                top_edge = 0;
            end
        elseif (r==dkR) %bottomEdge
            bottom_edge = 1;
            if (c==1)   %bottom_leftcorner
                bottom_left = 1;
                bottom_edge = 0;
            elseif (c==dkC) %bottom_rightcorner
                bottom_right = 1;
                bottom_edge = 0;
            end
        elseif (c==1) %left_edge
            left_edge = 1;
        elseif (c==dkC) %right_edge
            right_edge = 1;
        end

        if (num==0)
           if (top_left || bottom_left || left_edge)
               t0 = edge_s(r,c+1); %rightEdge
               one = 1;
           elseif (top_right || bottom_right || right_edge)
               t0 = edge_s(r,c-1); %leftEdge
               one = 1;
           else %everything in between (includes top_edge and bottom_edge)
               t1 = edge_s(r,c-1); %leftEdge
               t2 = edge_s(r,c+1); %rightEdge
           end
        elseif (num==45)
           if(top_edge || top_right || right_edge)
               t0=edge_s(r+1,c-1); %leftBottomCorner
               one = 1;
           elseif(bottom_edge || bottom_left || left_edge)
               t0 = edge_s(r-1,c+1); %rightTopCorner
               one=1;
           elseif(top_left || bottom_right)
               zero=1;
           else %everything in between
               t1 = edge_s(r-1,c+1); %rightTopCorner
               t2 = edge_s(r+1,c-1); %leftBottomCorner
           end
        elseif (num==90)
           if(top_edge || top_left || top_right)
               t0 = edge_s(r+1,c); %bottomEdge
               one=1;
           elseif(bottom_edge || bottom_left || bottom_right)
               t0 = edge_s(r-1,c); %topEdge
               one=1;
           else %everything in between (includes left_edge and right_edge)
               t1 = edge_s(r-1,c); %topEdge
               t2 = edge_s(r+1,c); %bottomEdge
           end
        elseif (num==135)
           if(top_edge || top_left || left_edge)
               t0 = edge_s(r+1,c+1); %rightBottomCorner
               one=1;
           elseif(bottom_edge || bottom_right || right_edge)
               t0 = edge_s(r-1,c-1); %leftTopCorner
               one=1;
           elseif(top_right || bottom_left)
               zero=1;
           else
               t1 = edge_s(r-1,c-1); %leftTopCorner
               t2 = edge_s(r+1,c+1); %rightBottomCorner
           end
        else
            edge_d(r,c) = NaN;
            zero = 1;
        end
        
% If Es (i,j) is smaller than at least one of its two neighbors along dk? , 
%   assign IN(i,j) = 0 (suppression); otherwise assign IN(i,j) = Es (i,j).

        i_n(r,c) = edge_s(r,c);
        
        if(one)
            if (edge_s(r,c) < t0)
                i_n(r,c) = 0;
            end
        elseif(zero)
            %nothing
            %i_n(r,c) = 0;
        else %both
            if (edge_s(r,c) < t1 || edge_s(r,c) < t2)
                i_n(r,c) = 0;
            end
        end
 
    end
end

i_n = uint8(i_n);
figure, imshow(i_n);



% HYSTERESIS THRESH
%3) edge localization
%   thinning wide edges to 1-pixel width, thresholding
%   Decide which local maxima in the filter?s output are edges
%   which are just caused by noise
%   thinning wide edges to 1-pixel width (nonmaximum suppression)
%   establishing the minimum value to declare a local maxima an
%   edge (thresholding)

%input is IN ,(the output of NONMAX_SUPPRESSION), 
% the edge orientation image, and two thresholds (?l and ?h ) such
% that ?l < ?h . For all edge points in IN:
% 1. Locate the next unvisited edge pixel IN (i,j) such that IN(i,j)> ?h;
% 2. Starting from IN (i,j), follow the chains of connected local
%   maxima, in both directions perpendicular to the edge normal, 
%   as long as IN(i,j)> ?l. Mark all visited points, 
%   and save a list of the locations of all points 
%   in the connected contour found.

%Note that hysteresis thresholding performs edge tracking. 
%It finds chains of connected contours.

t_l = 10;
t_h = 25;

[in_r,in_c] = size(i_n);

for r=1:in_r
    for c=1:in_c
        if (i_n(r,c) > t_h)
            bin_mtx(r,c) = 1;
        else
            bin_mtx(r,c) = 0;
        end
    end
end

save = zeros(size(i_n));
mat1 = [];
for r=1:in_r
    for c=1:in_c
        if (bin_mtx(r,c) == 1)
            [mat, save] = hyst(edge_d, i_n, r, c, t_l, save);
            mat1 = [mat1; mat];
        end        
    end
end

Out_edge = ones(size(i_n))*255;
[c1, ~] = size(mat1);
for i = 1:c1
    Out_edge(mat1(i, 1), mat1(i,2)) = 0;
end

figure, imshow(Out_edge);

function [edge_xy, search] = hyst(e_dk, out, x, y, low, search)
edge_xy = [x,y];
[out_r,out_c] = size(out);

r = x;
c = y;

search(r,c) = 1;

if (r==1)   %topEdge
    if(c==1)    %top_leftcorner
    elseif(c==out_c) %top_rightcorner
    end
elseif (r==out_r) %bottomEdge
    if (c==1)   %bottom_leftcorner
    elseif (c==out_c) %bottom_rightcorner
    end
elseif (c==1) %left_edge
elseif (c==out_c) %right_edge
    
else
    
    a = e_dk(r,c);
    
    if (a==0)
        t1 = out(r,c-1); %leftEdge
        t2 = out(r,c+1); %rightEdge
        if (search(r,c-1)==0)
            if (t1 > low)
                [edge_xy1, search] = hyst(e_dk, out, r, (c-1), low, search);
                edge_xy = [edge_xy1; [r, c-1]];
            end
        end
        if (search(r,c+1)==0)
            if (t2 > low)
                [edge_xy1, search] = hyst(e_dk, out, r, (c+1), low, search);
                edge_xy = [edge_xy1; [r, c+1]];
            end
        end
    elseif (a==45)
        t1 = out(r-1,c+1); %rightTopCorner
        t2 = out(r+1,c-1); %leftBottomCorner
        if (search(r-1,c+1)==0)
            if (t1 > low)
                [edge_xy1, search] = hyst(e_dk, out, r-1, (c+1), low, search);
                edge_xy = [edge_xy1; [r-1, c+1]];
            end
        end
        if (search(r+1,c-1)==0)
            if (t2 > low)
                [edge_xy1, search] = hyst(e_dk, out, r+1, (c-1), low, search);
                edge_xy = [edge_xy1; [r+1, c-1]];
            end
        end
    elseif (a==90)
        t1 = out(r-1,c); %topEdge
        t2 = out(r+1,c); %bottomEdge
        if (search(r-1,c)==0)
            if (t1 > low)
                [edge_xy1, search] = hyst(e_dk, out, r-1, (c), low, search);
                edge_xy = [edge_xy1; [r-1, c]];
            end
        end
        if (search(r+1,c)==0)
            if (t2 > low)
                [edge_xy1, search] = hyst(e_dk, out, r+1, c, low, search);
                edge_xy = [edge_xy1; [r+1, c]];
            end
        end
    elseif (a==135)
        t1 = out(r-1,c-1); %leftTopCorner
        t2 = out(r+1,c+1); %rightBottomCorner
        if (search(r-1,c-1)==0)
            if (t1 > low)
                [edge_xy1, search] = hyst(e_dk, out, r-1, (c-1), low, search);
                edge_xy = [edge_xy1; [r-1, c-1]];
            end
        end
        if (search(r+1,c+1)==0)
            if (t2 > low)
                [edge_xy1, search] = hyst(e_dk, out, r+1, (c+1), low, search);
                edge_xy = [edge_xy1; [r+1, c+1]];
            end
        end
    else %NaN
    end
end
end

% Test your algorithm on images ?Flowers.jpg?, ?Syracuse 01.jpg? 
%   and ?Syracuse 02.jpg?
% Try different ? (the standard deviation of the Gaussian), 
%   ?l (lower threshold) and ?h (higher threshold)
% Compare and evaluate your results.
% Test your algorithm on a favorite image of yours.