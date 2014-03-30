function batch_crop( dirName, positivesDir, negativesDir )
%\brief
% This function is used to CROP rectangles from FERET faces.

% Coordinates are read from .gnd files and images are read from .tif files.
% 	Coordinates are eye, nose and mouth coordinates rather than direct
% 	face rectangles.
%
%	Key function is rect_from_coords() where the face rectangle is formed
%	using face feature (eye,nose,mouth) coordinates.
%
%\params
%	dirName = directory where images and coordinate files do exist
%	positivesDir = where to save positives
%	negativesDir = where to save negatives/cropped from around the face
%
    global IMAGE_FILE_EXT;
    IMAGE_FILE_EXT    ='.tif';
    textFileExt     ='.gnd';
    
    allFileNames = dir([dirName '/*' textFileExt ]);
    
    if ~isdir(positivesDir)
        mkdir(positivesDir);
    end
    
    if ~isdir(negativesDir)
        mkdir(negativesDir)
    end
    
    for i=1:length(allFileNames)
        % First get file name and then filter out the extension
        fileName = allFileNames(i).name;
        fileName = fileName(1:(length(fileName)-length(textFileExt)));
        
        [positive negatives] = draw_face_in_image([dirName '/' fileName]);
        
        if max(positive(:)) > 0
            imwrite( positive, [positivesDir '/' fileName '_cropped.png' ] );
            
            for k=1:size(negatives,3)
                curNeg = uint8(negatives(:,:,k));
                if max(curNeg(:)) == 0
                    continue;
                end
                
                imwrite( curNeg, sprintf( '%s/%s%d_cropped.png', negativesDir, fileName, k) );
            end
            
        end
        
    end
    
end

function [positive negatives] = draw_face_in_image( imName )
%DRAW_FACE_IN_IMAGE Summary of this function goes here

    global IMAGE_FILE_EXT;

    % First read image
    im = imread([ imName IMAGE_FILE_EXT]);
    [leftEye rightEye noseTip mouth] = get_coords(imName,0);
    
    % If all coordinates do not exist, 
    %   image must be passed by...
    if leftEye == -1
        positive = 0;
        negatives = 0;
        return;
    end
    
    [positive negatives] = rect_from_coords(im, leftEye, rightEye, noseTip, mouth);
end

function [ positive negatives faceTopLeft faceBottomRight ] = rect_from_coords( im, leftEye, rightEye, noseTip, mouth )
%RECT_FROM_COORDS Summary of this function goes here
%   Detailed explanation goes here
    
    % Image must be padded, because some rotations may cause the face point
    % to move out of the frame
    padSize = 120;
    im = padarray(im,[padSize padSize],'symmetric', 'both');
    
    % Coordinates must be padded too
    leftEye     = leftEye + padSize;
    rightEye    = rightEye + padSize;
    noseTip     = noseTip + padSize;
    mouth       = mouth + padSize;
    
    % First find the HORIZONTAL axis of the FACE
    % theta: The angle of slope of eye axis (y=mx+b)
    theta = atan((leftEye(2)-rightEye(2))/(leftEye(1)-rightEye(1)));
    
    % distance bw eyes
    eyeDist = sqrt( (leftEye(1)-rightEye(1)).^2+(leftEye(2)-rightEye(2)).^2 );
    
    % The width of the face is 5 eyes. One eye is 2*a.
    %   eyeDistance is 4*a
    %   faceWidth is 10*a
    t = eyeDist/4;
    
    % Coordinates of eye centre is necessery for the next bunch of steps..
    eyeCenter = [(leftEye(1)+rightEye(1))/2 (leftEye(2)+rightEye(2))/2];
    const = 2.4;
    
    % Left and right ends of eye axis.
    eyeLeft     = [eyeCenter(1)-const*t*cos(theta) eyeCenter(2)-const*t*sin(theta)];
    eyeRight    = [eyeCenter(1)+const*t*cos(theta) eyeCenter(2)+const*t*sin(theta)];
    
    
    % ----------------------------
    %   Now we need the vertical axis of the face!
    
    % theta: The angle of slope of vertical axis (y=mx+b)
    theta = atan((noseTip(2)-mouth(2))/(noseTip(1)-mouth(1)));
    
    % distance bw eye center (eyeCenter) and vertical center (verticalCenter)
    dist = sqrt( (eyeCenter(1)-mouth(1)).^2+(eyeCenter(2)-mouth(2)).^2 );
    t=1*dist;
    
    %plot(vertCenter(1),vertCenter(2), 'o');
    const = 0.70;
    
    % Upper and lower points of the vertical face axis
    midUp     = [eyeCenter(1)-const*t*abs(cos(theta)) eyeCenter(2)-const*t*abs(sin(theta))];
    midBottom = [mouth(1)+const*t*abs(cos(theta)) mouth(2)+const*t*abs(sin(theta))];
    
    % Rotated/Correct face rectangle
    topLeftCorner       = face_corner(eyeLeft, midUp, theta);
    topRightCorner      = face_corner(eyeRight, midUp, theta);
    bottomLeftCorner    = face_corner(eyeLeft, midBottom, theta);
    bottomRightCorner   = face_corner(eyeRight, midBottom, theta);
    
    facePoints = [[topLeftCorner 1]' [topRightCorner 1]' [bottomLeftCorner 1]' [bottomRightCorner 1]'];
    
    % Align face with the correct rotation
    [rIm rPoints] = rotate_im(im, theta, facePoints);
    
    % These two points are sufficient to mark face
    faceTopLeft = round([rPoints(0*3+1) rPoints(0*3+2)]);
    faceBottomRight = round([rPoints(3*3+1) rPoints(3*3+2)]);
%     draw_rect(im, faceTopLeft, faceBottomRight);
    
    x_diff = abs(faceTopLeft(1) - faceBottomRight(1));
    y_diff = abs(faceTopLeft(2) - faceBottomRight(2));
    
    % Face window must be square
    %   If it is rectangular, extend short edges symmetrically
    if x_diff > y_diff
        step = (x_diff-y_diff)/2;
        faceTopLeft(2) = floor(faceTopLeft(2)-step);
        faceBottomRight(2) = floor(faceBottomRight(2)+step);
    else
        step = (y_diff-x_diff)/2;
        faceTopLeft(1) = floor(faceTopLeft(1)-step);
        faceBottomRight(1) = floor(faceBottomRight(1)+step);
    end
    
    % Rectangular to get cropped
    rect = [faceTopLeft(1) faceTopLeft(2) faceBottomRight(1)-faceTopLeft(1) faceBottomRight(2)-faceTopLeft(2)];
    
    [rect frame] = newRects(rect, wrev(size(rIm)));
    positive = imcrop(rIm, frame);
    
    positive = imcrop(positive, rect);
    negatives = 0;
    
    %negatives = get_negatives(im, rect);
    
    % Output must be valid
    if max(rect(:)) == 0
        positive    = 0;
        negatives   = 0;
    end
    
end

% First choose a wider box then crop it. This function widens 
%     the rect appropriately before rotation.
function [rect frame] = newRects( rect, origSize )
    % rect: Rectangle to crop
    % origSize: Original image size
    w1 = rect(1);
    w2 = origSize(1) - (rect(3)+rect(1));
    h1 = rect(2);
    h2 = origSize(2) - (rect(4)+rect(2));
    
    expand = min([w1,  w2, h1, h2]);
    
    y = rect; % Expanded rect coords
    y(1) = y(1) - expand;
    y(2) = y(2) - expand;
    y(3) = y(3) + expand + rect(1) - y(1);
    y(4) = y(4) + expand + rect(2) - y(2);
    
    y = round(y);
    
    frame = y;
    
    % Update the face coords in the cropped image
    rect(1) = expand;
    rect(2) = expand;
        
end

% Produce negatives from a given positive region
function negatives = get_negatives(im, rect)
    width = rect(3);
    
    % leftBox, rightBox, upperBox, lowerBox
    ltBox = [rect(1) rect(2)-width width width];
    rtBox = [rect(1) rect(2)+width width width];
    upBox = [rect(1)-width rect(2) width width];
    lwBox = [rect(1)+width rect(2) width width];
    
    % upperLeftBox, upperRightBox ...
    upLtBox = [rect(1)-width rect(2)-width width width];
    upRtBox = [rect(1)-width rect(2)+width width width];
    lwLtBox = [rect(1)+width rect(2)-width width width];
    lwRtBox = [rect(1)+width rect(2)+width width width];
    
    boxCoords = [ltBox; rtBox; upBox; lwBox; upLtBox; upRtBox; lwLtBox; lwRtBox];
    
    % Max valid negatives number: size(boxCoords,1) = 8 (ltBox ... )
    negatives = zeros(width+1, width+1, size(boxCoords,1)+1);
    
    for i=1:size(boxCoords,1)
        % Box coords must be valid
        if 0 >= min(min(boxCoords(i,:)))
            continue;
        end
        
        imW = size(im,2);
        imH = size(im,1);
        
        curCoords = boxCoords(i,:);
        
        if curCoords(1)+curCoords(3)>imW
            continue;
        end
        
        if curCoords(2)+curCoords(4)>imH
            continue;
        end
        
        negatives(:,:,i) = imcrop(im, boxCoords(i,:));
    end
end

function [ rotatedIm rotatedPoints ] = rotate_im( im, theta, points )
%ROTATE_�M Summary of this function goes here
%   Detailed explanation goes here
    im = double(im);
    
    if theta < 0
        phi = theta+pi/2;
    else
        phi = theta-pi/2;
    end
    
    % Transformation matrix for points
    hp = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    
    % Transformation matrix for image
    h = [cos(-phi) -sin(-phi) 0; sin(-phi) cos(-phi) 0; 0 0 1];
    h=inv(h);
    
    
    [xi yi] = meshgrid(1:size(im,2),1:size(im,1));
    
    % Rotate face points before matrix is inverted
    % NOTE!! -> hp\points = INV(hp)*points
    rotatedPoints = round( hp\points );
    
    xx = (h(1,1)*xi+h(1,2)*yi+h(1,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));
    yy = (h(2,1)*xi+h(2,2)*yi+h(2,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));
    
    rotatedIm = uint8(interp2(im,xx,yy));

end

function [leftEye rightEye noseTip mouth] = get_coords( imName, draw )
%GET_COORDINATES of face partss
%   Detailed explanation goes here
    
    % Read information about current file (Eye coordinates vs...)
    fid = fopen([imName '.gnd']);
    C=fread(fid,'*char');  
    C=C';
    fclose(fid);
    
    % Get image
    im  = imread([imName '.tif']);
    
    % Retrieve each coordinate through REGULAR EXPRESSION
    leftEye     = regexp(C, 'left_eye_coords= (\d+) (\d+)', 'tokens');
    rightEye    = regexp(C, 'right_eye_coords= (\d+) (\d+)', 'tokens');
    noseTip     = regexp(C, 'nose_tip_coords= (\d+) (\d+)', 'tokens');
    mouth       = regexp(C, 'mouth_center_coords= (\d+) (\d+)', 'tokens');
    
    % All data (leftEye, .. ) is needed to draw face region...
    if ( isempty(leftEye) || isempty(rightEye) || isempty(noseTip) || isempty(mouth) )
        display(sprintf('Insufficient DATA for image %s.tif! Image is passed by...', imName));
        leftEye = -1;
        return;
    end
    
    % Convert strings to double and separate X and Y
    leftEye     = [str2double(leftEye{1}{1})    str2double(leftEye{1}{2})];
    rightEye    = [str2double(rightEye{1}{1})   str2double(rightEye{1}{2})];
    noseTip     = [str2double(noseTip{1}{1})    str2double(noseTip{1}{2})];
    mouth       = [str2double(mouth{1}{1})      str2double(mouth{1}{2})];
    
    % Shall we draw points?
    if ( draw )
        figure;
        imshow(im);
        hold on;
        mark_point(leftEye);
        mark_point(rightEye);
        mark_point(noseTip);
        mark_point(mouth);
    end
end

function mark_point( point )
%MARK_PO�NT Summary of this function goes here
%   Detailed explanation goes here
    
    plot(point(1),point(2), 'o');
    
end

function corner = face_corner(p1, p2, beta)
    x1=p1(1);
    x2=p2(1);
    y1=p1(2);
    y2=p2(2);
    
    a = (x2-x1-y1*tan(beta)+y2*tan(beta))/(cos(beta)*(1+(tan(beta).^2)));
    
    corner = [x1+a*cos(beta) y1+a*sin(beta)];
end
