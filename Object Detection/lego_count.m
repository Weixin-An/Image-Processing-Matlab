% 使用方法
% 读取图像： I=imread('train01.jpg');
% 目标检测： [NumA, NumB] = lego_count( I )
function [NumA, NumB] = lego_count( I )
%% Procedure
% 1. Transform the image into lab color space, and extract red pixel block 
%    and blue pixel block respectively;
% 2. Extract edge, morphological operation (use to clean up the results 
%    of thresholding);
% 3. Use Watershed algorithm segment image;
% 4. Classify shapes of this image as square, rectangle, triangle and circle 
%    according to pixel block properties (area and perimeter). 
%    And count the object A, object B pixel blocks respectively 
figure, imshow(I),title('Original Image')
%% 1.Transform the image into lab color space, and extract red pixel block and blue pixel block respectively;
% Convert image from RGB color space to lab color space
cform = makecform('srgb2lab');
lab_he = applycform(I, cform);
% Extract the red and blue blocks of the image
segmented_blue = im2bw(lab_he(:,:,3), 0.45);
segmented_red = im2bw(lab_he(:,:,2), 0.60);
segmented_blue = double(lab_he) .* cat(3,cat(3, ~segmented_blue, ~segmented_blue), ~segmented_blue);
segmented_red = double(lab_he) .* cat(3,cat(3, segmented_red, segmented_red), segmented_red);
% figure, subplot(1, 2, 1), imshow(segmented_blue),  title('Blue block image')
% subplot(1, 2, 2), imshow(segmented_red), title('Red block image')

%% (Blue)2. Extract edge, morphological operation for blue blocks (use to clean up the results of thresholding);
% Extract edge by 'canny' operator
grey_blue = rgb2gray(uint8(segmented_blue));
[~, threshold] = edge(grey_blue, 'canny');
result_edge = edge(grey_blue,'canny', threshold * 0.2);
% Morphological operation (dilate and fill)
se = strel('square',7);
result_dilate = imdilate(result_edge, se);
result_fill = imfill(result_dilate, 'holes');
% figure, imagesc(result_fill), title('Results of morphological operation-blue')

%% (Blue)3. Use Watershed algorithm segment blue blocks;
% Compute the distance transform of the complement of the binary image and 
% complement the distance transform
D1 = -bwdist(~result_fill);% figure, imshow(D, [])
% Filter out tiny local minima using imextendedmin and then modify the 
% distance transform so that no minima occur at the filtered-out locations. 
mask = imextendedmin(D1, 2);
D2 = imimposemin(D1, mask);
% Compute the watershed transform
Ld2 = watershed(D2);
segmented_img_blue = result_fill;
% Set the color of watershed pixels in the original image as background-color
segmented_img_blue(Ld2 == 0) = 0;
% figure(4), imagesc(segmented_img_blue), title('Blue objects in the image');

%% (Red)2. Extract edge, morphological operation for red blocks (use to clean up the results of thresholding);
% Extract edge by 'canny' operator
grey_red = rgb2gray(uint8(segmented_red));
[~, threshold] = edge(grey_red, 'canny');
result_edge = edge(grey_red,'canny', threshold * 0.4);
% figure, imagesc(BWs), title('binary gradient mask');
% Morphological operation (dilate and fill)
se = strel('square',5);
result_dilate = imdilate(result_edge, se);
result_fill = imfill(result_dilate, 'holes');
% figure, imagesc(BWdfill), title('Results of morphological operation-red')

%% (Red)3. Use Watershed algorithm segment red blocks;
% Compute the distance transform of the complement of the binary image and 
% complement the distance transform
D = -bwdist(~result_fill);
% filter out tiny local minima using imextendedmin and then modify the 
% distance transform so that no minima occur at the filtered-out locations. 
mask = imextendedmin(D,2);
D2 = imimposemin(D,mask);
% Compute the watershed transform
Ld2 = watershed(D2);
segmented_img_red = result_fill;
% Set the color of watershed pixels in the original image as background-color
segmented_img_red(Ld2 == 0) = 0;
% figure(4), imagesc(segmented_img_red), title('Red objects in the image');

%% 4. Classify shapes of this image as square, rectangle, triangle and circle 
%%    according to pixel block properties (area and perimeter). Count the 
%%    object A, object B pixel blocks respectively.
% Blue square lego perimeter range: [blue_thre_min, blue_thre_max]
% Red square lego perimeter range: [red_thre_min, red_thre_max]
blue_thre_min = 672;
blue_thre_max = 1400;
red_thre_min = 430;
red_thre_max = 650;
% Color = 1, 2. 1 -- blue(ObjectA), 2 -- red(ObjectB)
NumA = Classify_and_Count(segmented_img_blue,blue_thre_min, blue_thre_max, 1);
NumB = Classify_and_Count(segmented_img_red, red_thre_min, red_thre_max, 2);
end


%% There are two subfunctions
%% The first one
% Classify shapes of this image as square, rectangle, triangle and circle
% Count the number of ObjectA and ObjectB
function [ number_of_legos ] = Classify_and_Count( segmented_img, thre1, thre2, Color)
%% Classify shapes of this image as square, rectangle, triangle and circle
% Reference website:
% https://ww2.mathworks.cn/matlabcentral/answers/116793-how-to-classify-shapes-of-this-image-as-square-rectangle-triangle-and-circle
% Get pixel block bounding box, connected region and serial number of connected region
[B,L,N] = bwboundaries(segmented_img);
% Get stats including perimeter, area and metrics for each shape
stats=  regionprops(L, 'Centroid', 'Area', 'Perimeter');
Centroid = cat(1, stats.Centroid);
Perimeter = cat(1,stats.Perimeter);
Area = cat(1,stats.Area);
% The closer to 1, the closer CircleMetric is to a circle
CircleMetric = (Perimeter.^2)./(4*pi*Area);  %circularity metric
% 1 -- a square, 0 -- not a square
SquareMetric = NaN(N,1);
% 1 -- a rectangle, 0 -- not a rectangle
TriangleMetric = NaN(N,1);
% For each boundary, fit to bounding box, and calculate parameters
for k=1:N
   boundary = B{k};
   [rx,ry,boxArea] = minboundrect( boundary(:,2), boundary(:,1));
   % Get width and height of bounding box
   width = sqrt( sum( (rx(2)-rx(1)).^2 + (ry(2)-ry(1)).^2));
   height = sqrt( sum( (rx(2)-rx(3)).^2+ (ry(2)-ry(3)).^2));
   aspectRatio = width/height;
   if aspectRatio > 1  
       % Make aspect ratio less than unity
       aspectRatio = height/width;  
   end
   % Aspect ratio of box sides
   SquareMetric(k) = aspectRatio;  
   % Filled area vs box area
   TriangleMetric(k) = Area(k)/boxArea;  
end
% Define thresholds for each metric
% Do in order of circle, triangle, square, rectangle to avoid assigning the
% same shape to multiple objects
isCircle =   (CircleMetric < 1.1);
isTriangle = ~isCircle & (TriangleMetric < 0.6);
isSquare =   ~isCircle & ~isTriangle & (SquareMetric > 0.75);
isRectangle= ~isCircle & ~isTriangle & ~isSquare;  %rectangle isn't any of these
%% Count the number of ObjectA and ObjectB
RGB = label2rgb(L);
if Color == 2
    figure, imshow(RGB), title('Red objects'); hold on;
else
    figure, imshow(RGB), title('Blue objects'); hold on;
end
number_of_legos = 0;
for k=1:N
   % Display metric values and shape name next to an object
    if Color == 2
        if isSquare(k) == 1 && Perimeter(k) > thre1 && Perimeter(k) < thre2
           text( Centroid(k,1)-20, Centroid(k,2)+20, 'Square');
           number_of_legos = number_of_legos + 1;
        end
    else
        if isRectangle(k) == 1 && Perimeter(k) > thre1 && Perimeter(k) < thre2
           text( Centroid(k,1)-20, Centroid(k,2)+20, 'Rectangle');
           number_of_legos = number_of_legos + 1;
        end
    end
end
end

%% The second one
% minboundrect: Compute the minimal bounding rectangle of points in the plane
% Reference:
% John D'Errico (2020). A suite of minimal bounding objects 
% (https://www.mathworks.com/matlabcentral/fileexchange/34767-a-suite-of-minimal-bounding-objects), 
% MATLAB Central File Exchange. Retrieved December 17, 2020.
function [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
% usage: [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
%
% arguments: (input)
%  x,y - vectors of points, describing points in the plane as
%        (x,y) pairs. x and y must be the same lengths.
%
%  metric - (OPTIONAL) - single letter character flag which
%        denotes the use of minimal area or perimeter as the
%        metric to be minimized. metric may be either 'a' or 'p',
%        capitalization is ignored. Any other contraction of 'area'
%        or 'perimeter' is also accepted.
%        DEFAULT: 'a'    ('area')
%
% arguments: (output)
%  rectx,recty - 5x1 vectors of points that define the minimal
%        bounding rectangle.
%
%  area - (scalar) area of the minimal rect itself.
%
%  perimeter - (scalar) perimeter of the minimal rect as found
% 
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com


% default for metric
if (nargin<3) || isempty(metric)
  metric = 'a';
elseif ~ischar(metric)
  error 'metric must be a character flag if it is supplied.'
else
  % check for 'a' or 'p'
  metric = lower(metric(:)');
  ind = strmatch(metric,{'area','perimeter'});
  if isempty(ind)
    error 'metric does not match either ''area'' or ''perimeter'''
  end
  % just keep the first letter.
  metric = metric(1);
end

% preprocess data
x=x(:);
y=y(:);

% not many error checks to worry about
n = length(x);
if n~=length(y)
  error 'x and y must be the same sizes'
end

% start out with the convex hull of the points to
% reduce the problem dramatically. Note that any
% points in the interior of the convex hull are
% never needed, so we drop them.
if n>3
  edges = convhull(x,y);

  % exclude those points inside the hull as not relevant
  % also sorts the points into their convex hull as a
  % closed polygon
  
  x = x(edges);
  y = y(edges);
  
  % probably fewer points now, unless the points are fully convex
  nedges = length(x) - 1;
elseif n>1
  % n must be 2 or 3
  nedges = n;
  x(end+1) = x(1);
  y(end+1) = y(1);
else
  % n must be 0 or 1
  nedges = n;
end

% now we must find the bounding rectangle of those
% that remain.

% special case small numbers of points. If we trip any
% of these cases, then we are done, so return.
switch nedges
  case 0
    % empty begets empty
    rectx = [];
    recty = [];
    area = [];
    perimeter = [];
    return
  case 1
    % with one point, the rect is simple.
    rectx = repmat(x,1,5);
    recty = repmat(y,1,5);
    area = 0;
    perimeter = 0;
    return
  case 2
    % only two points. also simple.
    rectx = x([1 2 2 1 1]);
    recty = y([1 2 2 1 1]);
    area = 0;
    perimeter = 2*sqrt(diff(x).^2 + diff(y).^2);
    return
end
% 3 or more points.

% will need a 2x2 rotation matrix through an angle theta
Rmat = @(theta) [cos(theta) sin(theta);-sin(theta) cos(theta)];
% get the angle of each edge of the hull polygon.
ind = 1:(length(x)-1);
edgeangles = atan2(y(ind+1) - y(ind),x(ind+1) - x(ind));
% move the angle into the first quadrant.
edgeangles = unique(mod(edgeangles,pi/2));
% now just check each edge of the hull
nang = length(edgeangles);
area = inf;
perimeter = inf;
met = inf;
xy = [x,y];
for i = 1:nang
  % rotate the data through -theta 
  rot = Rmat(-edgeangles(i));
  xyr = xy*rot;
  xymin = min(xyr,[],1);
  xymax = max(xyr,[],1);
  
  % The area is simple, as is the perimeter
  A_i = prod(xymax - xymin);
  P_i = 2*sum(xymax-xymin);
  
  if metric=='a'
    M_i = A_i;
  else
    M_i = P_i;
  end
  % new metric value for the current interval. Is it better?
  if M_i<met
    % keep this one
    met = M_i;
    area = A_i;
    perimeter = P_i;
    
    rect = [xymin;[xymax(1),xymin(2)];xymax;[xymin(1),xymax(2)];xymin];
    rect = rect*rot';
    rectx = rect(:,1);
    recty = rect(:,2);
  end
end
end 

