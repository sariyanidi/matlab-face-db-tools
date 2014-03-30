function rotateAndCrop(imagesDir, coordsFile, destDir)
% This is a script creating rotated and cropped faces from face database
% and its face coordinates.
%
% Until now it works separately for multipie and LFIW db.

    %inPlaneRots= [0 30 60];
    inPlaneRots = [330];
    
    profileFacesRotatedCrop(imagesDir, coordsFile, destDir, inPlaneRots);
    %lfiwRotatedCrop(imagesDir, coordsFile, destDir, inPlaneRots);
    %mpRotatedCrop(baseDir, coordsDir, destDir, origImSize, inPlaneRots);
end

function profileFacesRotatedCrop(imagesDir, coordsFile, destDir, inPlaneRots)
% Crop profile faces labelled by us
%
% @imagesDir:   Folder of images
% @coordsFile:  Text file where all image coordinates do exist
% @destDir:   	Destination directory of base
% @inPlaneRots: A vector where in plane rotations are given in degrees
% 
% 
% The steps accomplished using this script
% ------------------------------------------------------------
% 
% (For each image in directory...)
%   1. Get filename, pose and rect coordinates from text file
%   2. Find and read related image
%   3. Find and crop largest rectangle in image -> tempImage
%   4. Extend image, mirror symmetry -> another tempImage
%   5. Rotate image "numel(inPlaneRots)" times
%   6. Save rotated images
% 
% ------------------------------------------------------------
	global FACE_SIZE;
	global RANDOMIZE;

	FACE_SIZE = [24 24];
	RANDOMIZE = false;

	if ~isdir(imagesDir)
		disp('Invalid input image directory!');
	end
	
	totalCropped = zeros(3,1);
	
	clc
	tic
	
	fid = fopen(coordsFile, 'r');
    
    if -1 == fid
        disp('Could not open textfile including face coordinate!');
        return;
    end
    
    line = fgetl(fid);
    
    while line ~= -1
        % read coordinates and pose
        temp = textscan(line,'%s\t%d\t%d\t%d\t%d\t%d',1);
        
        try
            filename = cell2mat(temp{1}); temp(1) = []; % Convert cell to string
            
            % read rectangle
            temp = cell2mat(temp); % Convert cell to integer

            % this is the rectangle to crop
            cropR = [temp(1) temp(2) temp(3) temp(4)]; % original face rectangle
            
            % pose and face size of this sample
            pose = temp(5);
            faceSize = temp(4);
        
            % read image
            im = imread([imagesDir filesep filename]); % is target image there?

            % warn in every 250 crops
            if mod(sum(totalCropped),250) == 1
                disp(num2str(sum(totalCropped)));
            end
            
            padSize = faceSize+1;
            
            % We are enhancing this rectangle because we had cropped it too
            % strict! ATTENTION!! HARDCODED :)
            cropR(1) = cropR(1)-faceSize*0.15;
            cropR(2) = cropR(2)-faceSize*0.15;
            cropR(3) = cropR(3)+faceSize*0.25;
            cropR(4) = cropR(4)+faceSize*0.25;
            
            % pad image to have a rectangle enough to crop
            im = padarray(im, [double(padSize) double(padSize)], 'symmetric', 'both');
 
            % new crop rectangle after the extension of the face box
            largeRect(1) = cropR(1);
            largeRect(2) = cropR(2);
            largeRect(3) = cropR(3)+2*padSize;
            largeRect(4) = cropR(4)+2*padSize;

            % crop the face with a large padding, to be able to rotate it
            paddedFaceIm = imcrop(im, largeRect);
            
            % the face is always at the center of this cropped image
            cropR = ones(1,4);
            cropR(1) = floor(size(paddedFaceIm,1)/3)+1;
            cropR(2) = floor(size(paddedFaceIm,1)/3)+1;
            cropR(3) = floor(size(paddedFaceIm,1)/3)+1;
            cropR(4) = floor(size(paddedFaceIm,1)/3)+1;
            
            % we are saving profile faces one sided
            if (pose > 3)
                paddedFaceIm = flipdim(paddedFaceIm,2);
                pose = 6-pose;
            end
            
            if  pose == 3 
                line = fgetl(fid); % Get next face data
                continue;
            end
            
            % images will be saved as 00001.png ... etc
            newFileName = sprintf('%05d.png',totalCropped(pose));
            
            %saveRotated(paddedFaceIm, cropR, pose, destDir, newFileName, inPlaneRots);
            saveSlided(paddedFaceIm, cropR, pose, destDir, newFileName);

            totalCropped(pose) = totalCropped(pose)+1;
            
            line = fgetl(fid); % Get next face data
        catch e
            disp(['Error: ' e.message ]);
            line = fgetl(fid); % Get next face data
        end
    end
end

function lfiwRotatedCrop(imagesDir, coordsFile, destDir, inPlaneRots)
% Crop Labeled Faces In the Wild face database
%
% @imagesDir:   Folder of iamges
% @coordsFile:  Text file where all image coordinates do exist
% @destDir:   	Destination directory of base
% @inPlaneRots: A vector where in plane rotations are given in degrees
% 
% 
% The steps accomplished using this script
% ------------------------------------------------------------
% 
% (For each image in directory...)
%   1. Get filename, pose and rect coordinates from text file
%   2. Find and read related image
%   3. Find and crop largest rectangle in image -> tempImage
%   4. Extend image, mirror symmetry -> another tempImage
%   5. Rotate image 12 times (0deg, 30deg, 60deg, ... 330deg)
%   6. Save rotated images
% 
% ------------------------------------------------------------

    global IPR_RESOLUTION;
    global FACE_SIZE;
    global RANDOMIZE;
    
    IPR_RESOLUTION = 30;
    FACE_SIZE = [24 24];
    RANDOMIZE = true;
    
    if ~isdir(imagesDir)
        disp('Invalid image directory!');
        return; % 
    end
    
    totalCropped = 0;
    totalSkipped = 0;
    
    clc
    tic
    
    fid = fopen(coordsFile, 'r');
    
    if -1 == fid
        disp('Could not open textfile including face coordinates! Quitting...');
        return;
    end
    
    fid = fopen(coordsFile);
    line = fgetl(fid);
    
    disp('Cropping LFIW faces...');

    while line ~= -1
        % Read filename, coordinates and pose for each line
        temp = textscan(line,'%s\t%d\t%d\t%d\t%d\t%d',1);

        filename = cell2mat(temp{1}); temp(1) = []; % Convert cell to string
        temp = cell2mat(temp); % Convert cell to int
        
        cropR = [temp(1) temp(2) temp(3) temp(4)]; % original face rectangle
        pose = temp(5);
        
        try
            im = imread([imagesDir filesep filename]); % is target image there?
        catch e
            disp(e.message);
            continue;
        end
        
        if length(size(im)) == 3
            im = rgb2gray(im);
        end
        
        imSize = size(im); % Is needed for several things
        [looseR cropR] = getLooseRect( cropR, imSize ); % Get loose rectangle and update face rectangle
        
        % Some face images are very close to borders so rotation becomes
        % impossible, pad image to rotate
        [im cropR] = padImage(im, looseR, cropR);
        
        % Here beats the HEART of this program
        try
            saveRotated( im, cropR, pose, destDir, filename, inPlaneRots );
        catch e
            disp(e.identifier);
            disp(e.message);
            totalSkipped = totalSkipped+1;
            continue; % Don't try to fix it just skip it, there are thousands of faces...
        end
        
        line = fgetl(fid); % Get next face data
        totalCropped = totalCropped+1;
    end
    
    disp([num2str(totalCropped) ' faces are cropped, ' num2str(totalSkipped) ' faces are skipped.']);
end

function mpRotatedCrop(baseDir, coordsDir, destDir, origImSize, inPlaneRots)
% Crop Multi-PIE faces to use for rotation invariant FACE DETECTION
% classifier
%
% @baseDir:     Root folder of multi-pie db. Folders must follow MultiPIE
%               folder hierarchy! (see guide.pdf of MultiPIE db)
% @coordsDir:   Directory where of face coordinates. No folder hierarchy is
%               needed, just text files stacked within coordsDir folder.
% @destDir:   	Destination directory (Subfolders will be created here)
% @origImSize:  Original Multi-PIE image size. Default value is handed as
%               640x480
% @inPlaneRots: A vector where in plane rotations are given in degrees
%
% The faces are saved
%   * in different in plane rotations:
%     N separate folders are created for N separate in-plane rotations
%     (e.g. 30 deg, 60 deg, ... 330 deg rotations-> 12 subfolders). 
%     Rotation resolution is stated below in IPR_RESOLUTION parameter.
%     Output folder will be created next to baseDir
%
% This program mantains:
%   * Small random rotations before cropping
%   * Small random pixel shifts
%   * Small random scalings
%   ... so that the classifier is trained more robustly
% 
% ------------------------------------------------------------------
%
%   Author:     Evangelos Sariyanidi
%   E-mail:     sariyanidi@gmail.com
%   Date:       Mon, 27 Sep 2010 07:38:23 +0000
%   URL:        http://www.sariyanidi.com
%
% ------------------------------------------------------------------

    if 3 == nargin % Default image size in MultiPIE is ...
        origImSize = [480 640];
    end
    
    global IPR_RESOLUTION;
    global FACE_SIZE;
    
    IPR_RESOLUTION = 30; % IPR: In plane rotation
    FACE_SIZE = [24 24]; % Default face size is 24x24
    RANDOMIZE = true;    % Random shit, scale and rotation
    
    % Create destination folder if doesn't exist
    if ~isdir(destDir)
        try % Do we have access?
            mkdir(destDir)
        catch e
            disp(e.identifier);
            disp(e.message);
            return;
        end
    end
    
    root = dir(baseDir); root(2) = []; root(1) = [];
    totalCropped = 0;
    totalSkipped = 0; % Count how many faces are failed to crop
    
    clc
    tic
    for s=1:length(root) % SESSION
        session = root(s).name;

        % currentSessionPath, currentSessionTree
        sessPath = fullfile(baseDir, session, 'multiview');
        sessTree = dir(sessPath); sessTree(2) = []; sessTree(1) = [];
        
        for pe=1:length(sessTree) % PEOPLE
            person = sessTree(pe).name;

            % currentPersonPath, currentPersonTree
            personPath = fullfile(sessPath, person);
            personTree = dir(personPath); personTree(2) = []; personTree(1) = [];

            for e=1:length(personTree) % EXPRESSIONS
                express = personTree(e).name;

                expressPath = fullfile(personPath, express);
                expressTree = dir(expressPath); expressTree(2) = []; expressTree(1) = [];

                expressTree(13) = []; expressTree(6) = [];% these poses won't be used

                for po=1:length(expressTree) % POSES
                    pose = expressTree(po).name;

                    % currentPosePath, currentPoseTree
                    posePath = fullfile(expressPath, pose);
                    poseTree = dir(posePath); poseTree(2) = []; poseTree(1) = [];

                    for l=1:length(poseTree) % LIGHTING CONDITIONS
                        imageName = poseTree(l).name;

                        imagePath = fullfile(posePath, imageName);
                        
                        try % Any file-existence or mistaken coordinate errors will be caught here
                            coordFile = [person '_01_01_' pose(1:2) pose(4) '_00c.txt']; % File to pick face coordinates
                            coordPath = fullfile(coordsDir, coordFile);
                            
                            coords = load(coordPath);
                            image = imread(imagePath);
                            
                            [image coords] = preprocess( image, coords, po );
                            [looseR cropR] = getRects( coords, po, size(image), origImSize );
                        catch e
                            disp(e.identifier);
                            disp(e.message);
                            totalSkipped = totalSkipped+1;
                            continue; % Don't try to fix it just skip it, there are thousands of faces...
                        end
                        
                        % Add random rotations, scalings and shiftings to face rect.
                        if RANDOMIZE
                            image = randRotate(image, looseR); % Random rotation is like preprocessing
                            cropR = randScale(cropR, size(image)); % Scale the box randomly
                            cropR = randShift(cropR, size(image)); % Shift it...
                        end
                        
                        % Here beats the HEART of this program
                        try
                            saveRotated( image, cropR, pose, destDir, imageName, inPlaneRots );
                        catch e
                            disp(e.identifier);
                            disp(e.message);
                            totalSkipped = totalSkipped+1;
                            continue; % Don't try to fix it just skip it, there are thousands of faces...
                        end
                        
                        totalCropped = totalCropped+1;
                    end % LIGHTING
                end % POSES
            end % EXPRESSIONS
        end % PEOPLE
    end % SESSION
    
    % Results...
    toc
    disp(['Cropped faces: ' totalCropped '; Skipped faces: ' totalSkipped]);
    disp(['Check folder "' destDir '" for results.']);
end

function [im cropR] = padImage(im, looseR, cropR)
    padSize = round(size(im)./2);
    
    im = imcrop(im, looseR);
    im = padarray(im, padSize, 'symmetric', 'both');

    cropR(1) = cropR(1)+padSize(2);
    cropR(2) = cropR(2)+padSize(1);
end

function saveRotated(baseIm, rect, pose, destDir, imageName, inPlaneRots)
% Save the face after rotating it through many angles.

    global FACE_SIZE;
    
    poseDir = fullfile(destDir, ['pose_' num2str(pose)]);
    
    if ~isdir(poseDir) % Poses will be separately saved
        mkdir(poseDir)
    end
    
    % Handle every single rotation...
    for curRot=inPlaneRots
        im = baseIm; % duplicate it
        
        % Rotation dir, each rotated image will be saved to sep. folder
        rotDir = fullfile(poseDir, [num2str(curRot) '_DEG_IPR']);
        
        if ~isdir(rotDir)
            mkdir(rotDir)
        end
        
        im = imrotate(im, curRot, 'bicubic', 'crop');
        im = imcrop(im, rect);
        im = imresize(im, FACE_SIZE);
        
        % Create image name
        [tmp name ext]  = fileparts(imageName);
        nImageName      = sprintf('%s_%d_%d%s', name, pose, curRot, '.png'); % new image name
        nImagePath      = fullfile(rotDir, nImageName);
        
        % gray-level images
        if length(size(im)) == 3
            im = rgb2gray(im);
        end
        
        imwrite(im, nImagePath); % save image
    end
end

function rect = randShift(rect, imSize)
% Add random shifts to face rectangle
    
    global FACE_SIZE;
    
    shiftX = round(1.5*(imSize(1)/FACE_SIZE(1))*normrnd(0,0.30)); % parameters are picked experimentally
    shiftY = round(1.5*(imSize(1)/FACE_SIZE(1))*normrnd(0,0.30));
    rect(1) = rect(1) + shiftX;
    rect(2) = rect(2) + shiftY;
    
    % Make sure rectangle does not violate boundaries 
    rect = validateRect(rect, imSize);
end

function rect = randScale(rect, imSize)
% Add random scaling to face rectangle
    
    global FACE_SIZE;
    
    scale = round(1.5*(imSize(1)/FACE_SIZE(1))*normrnd(0,0.30)); % parameters are picked experimentally
    rect(1) = rect(1)-scale;
    rect(2) = rect(2)-scale;
    rect(3) = rect(3)+scale;
    rect(4) = rect(4)+scale;
    
    % Make sure rectangle does not violate boundaries
    rect = validateRect(rect, imSize);
end

function rect = validateRect(rect, imSize)
% Make sure rectangle does not violate boundaries
    if rect(1) < 1
        rect(1) = 1;
    elseif rect(1) > imSize(2)
        rect(1) = imSize(2);
    end
    
    if rect(2) < 1
        rect(2) = 1;
    elseif rect(2) > imSize(1)
        rect(2) = imSize(1);
    end
end

function im = randRotate( im, looseR )
% Take the original image, crop it so that the center of the face becomes
% the center of the image and rotate it randomly.
% 
% @im:              image to rotate
% @looseR:          looseRectangle, the largest rectangle containing the 
%                   centered face;
% @IPR_RESOLUTION:  what is the sensitivity of inplane rotation? The random
%                   rotation should not exceed IPR_RESOLUTION/2
%                   (e.g. if we are dealing with 30 deg rotation at each
%                   portion, random rotation should not exceed +-15 deg.
    
    global IPR_RESOLUTION;
    
    res = IPR_RESOLUTION/2;
    im = imcrop(im, looseR); % Crop to have the face centered before rotation
    randRot = round(res*normrnd(0,0.36)); % Artificial random rotation

    % Random rotation is limited to +-15 degrees
    if (randRot<-res)
        randRot = -res;
    elseif (randRot>res)
        randRot = res;
    end
    
    im = imrotate(im, randRot, 'bicubic', 'crop');
end

% Fonksiyonun ici Birkanin kodundan alinmistir
function [imge, nokta] = preprocess(imge, nokta, p)
    switch p
        case {2, 3, 4, 5, 9, 10}
            aci = atan((nokta(1,2)-nokta(3,2)) / (nokta(3,1)-nokta(1,1)));
            aci_deg = (180 * aci) / pi;
            imge = imrotate(imge, -aci_deg, 'bicubic', 'crop');
            nokta = nokta_dondur(nokta', aci, 640, 480);
            nokta = nokta';
    end
end

% Fonksiyonun ici disi Birkanin kodundan alinmistir
function donuk = nokta_dondur(nokta, aci, iw, ih)
    R = [cos(aci) sin(aci); -sin(aci) cos(aci)];

    temp = repmat([iw/2;ih/2], 1, size(nokta,2));
    nokta = temp - nokta;
    nokta(1,:) = -nokta(1,:);

    donuk = R * nokta;
    donuk(1,:) = -donuk(1,:);
    donuk = temp - donuk;
end

function [looseR cropR] = getLooseRect( cropR, imSize )
    % nokta: Rectangle to crop
    % origSize: Original image size
    % p: Pose identifier, needed for the function basicRect()
    %
    % looseR: Loose rectangle; the largest possible rectangle within image
    %         centered at the center of the face
    % cropR : Rectangle to crop -> will be randomized in later functions.
    
    %cropR(1) = cropR(1)-round((cropR(4)-cropR(3))/2);
    %cropR(3) = cropR(4);
    
    % All images in database have been scaled down (approx to 10% of
    % original) for total size reduction (400gb->~3gb).
    % Face coordinates are scaled down here too.
    % Works fine for original image inputs too!
    
    w1 = cropR(1);
    w2 = imSize(1) - (cropR(3)+cropR(1));
    h1 = cropR(2);
    h2 = imSize(2) - (cropR(4)+cropR(2));
    
    expand = min([w1,  w2, h1, h2]);
    
    y = cropR; % Expanded rect coords, will be needed for rotations...
    y(1) = y(1) - expand;
    y(2) = y(2) - expand;
    y(3) = y(3) + expand + cropR(1) - y(1);
    y(4) = y(4) + expand + cropR(2) - y(2);
    
    looseR = round(y); % The largest rectangle within image centered at face
    
    % Update the face coords in the cropped image
    cropR(1) = expand;
    cropR(2) = expand;
end


function [looseR cropR] = getRects( nokta, p, imSize, origImSize )
    % nokta: Rectangle to crop
    % origSize: Original image size
    % p: Pose identifier, needed for the function basicRect()
    %
    % looseR: Loose rectangle; the largest possible rectangle within image
    %         centered at the center of the face
    % cropR : Rectangle to crop -> will be randomized in later functions.
    
    try
        cropR = basicRect( nokta, p );
    catch e
        rethrow(e);
    end
    
    cropR(1) = cropR(1)-round((cropR(4)-cropR(3))/2);
    cropR(3) = cropR(4);
    
    % All images in database have been scaled down (approx to 10% of
    % original) for total size reduction (400gb->~3gb).
    % Face coordinates are scaled down here too.
    % Works fine for original image inputs too!
    scaleRatio = imSize(2) / origImSize(2);
    cropR = round( cropR*scaleRatio );
    
    w1 = cropR(1);
    w2 = origImSize(1) - (cropR(3)+cropR(1));
    h1 = cropR(2);
    h2 = origImSize(2) - (cropR(4)+cropR(2));
    
    expand = min([w1,  w2, h1, h2]);
    
    y = cropR; % Expanded rect coords, will be needed for rotations...
    y(1) = y(1) - expand;
    y(2) = y(2) - expand;
    y(3) = y(3) + expand + cropR(1) - y(1);
    y(4) = y(4) + expand + cropR(2) - y(2);
    
    looseR = round(y); % The largest rectangle within image centered at face
    
    % Update the face coords in the cropped image
    cropR(1) = expand;
    cropR(2) = expand;
end

function saveSlided(baseIm, rect, pose, destDir, imageName)
% Save the face after rotating it through many angles.
    
    global FACE_SIZE;

    poseDir = fullfile(destDir, ['pose_' num2str(pose)]);
    poseDirSymm = fullfile(destDir,['pose_' num2str(6-pose)] );
    
    if ~isdir(poseDir) % Poses will be separately saved
        mkdir(poseDir)
    end    
    
    if ~isdir(poseDirSymm) % Poses will be separately saved
        mkdir(poseDirSymm)
    end
    
    im = baseIm; % duplicate it
    
    slidedRect = limitSafeSlidedRect(rect,size(baseIm));
    
    % Create image name
    [~, name, ~]  = fileparts(imageName);
    nImageName      = sprintf('%s-%d%s',name, slidedRect(1),'.png'); % new image name
    nImagePath      = fullfile(poseDir, nImageName);
    nImagePathSymm  = fullfile(poseDirSymm, nImageName);
    
    im = imcrop(im, slidedRect);
    im = imresize(im, FACE_SIZE);
    
    if length(size(im)) == 3
        im = rgb2gray(im);
    end
    
    imwrite(im, nImagePath); % save image
    imwrite(flipdim(im,2), nImagePathSymm);
    
end

function y = limitSafeSlidedRect(rect, imSize)
    actualFaceSize = rect(3);
    
    %slideRatio = normrnd(0,0.40);
    % Generate a random number between -1 and 1
    a = -1;
    b = 1;
    slideRatio = a + (b-a).*rand(1,1);
    
    if slideRatio < -1
        slideRatio = -1;
    elseif slideRatio > 1
        slideRatio = 1;
    end
    
    slideValue = round(actualFaceSize*slideRatio/3.4);
    
    leftReach = slideValue+rect(1);
    rightReach = slideValue+rect(1)+rect(3);
    
    if leftReach < 1
        rect(1) = 1;
        y = rect;
        return;
    elseif rightReach>imSize(2)
        rect(1) = imSize(2)-rect(3);
        y = rect;
        return;
    end
    
    % Otherwise
    rect(1) = rect(1)+slideValue;
    y = rect;
end

function basicR = basicRect( nokta, p )
    try % If points are not marked correctly, program shouldn't stop!
        switch p
            case 1
                dikd = [nokta(6,1)+30 (nokta(1,2)-25) (nokta(3,1)-nokta(6,1)-30) (nokta(5,2)-nokta(2,2)+11)];
                dikd2 = [nokta(6,1)-80 (nokta(1,2)-130) (nokta(3,1)-nokta(6,1)+85) (nokta(5,2)-nokta(2,2)+150)];
            case 2
                dikd = [nokta(7,1)+15 (nokta(3,2)-25) (nokta(4,1)-nokta(7,1)+10) (nokta(6,2)-nokta(2,2)+11)];
                dikd2 = [nokta(7,1)-60 (nokta(3,2)-110) (nokta(4,1)-nokta(7,1)+105) (nokta(6,2)-nokta(2,2)+130)];
            case 3
                dikd = [nokta(7,1)+15 (nokta(3,2)-25) (nokta(4,1)-nokta(7,1)+33) (nokta(6,2)-nokta(2,2)+11)];
                dikd2 = [nokta(7,1)-50 (nokta(3,2)-120) (nokta(4,1)-nokta(7,1)+130) (nokta(6,2)-nokta(2,2)+120)];
            case 4
                dikd = [nokta(7,1)+5 (nokta(3,2)-25) (nokta(8,1)-nokta(7,1)-9) (nokta(6,2)-nokta(2,2)+11)];
                dikd2 = [nokta(7,1)-30 (nokta(3,2)-120) (nokta(8,1)-nokta(7,1)+60) (nokta(6,2)-nokta(2,2)+130)];
            case 5
                dikd = [nokta(4,1)-5 (nokta(3,2)-25) (nokta(7,1)-nokta(4,1)-25) (nokta(6,2)-nokta(2,2)+11)];
                dikd2 = [nokta(4,1)-25 (nokta(3,2)-140) (nokta(7,1)-nokta(4,1)+100) (nokta(6,2)-nokta(2,2)+150)];
            case 6
                dikd = [nokta(3,1) (nokta(2,2)-25) (nokta(6,1)-nokta(3,1)-25) (nokta(5,2)-nokta(1,2)+11)];
                dikd2 = [nokta(3,1)-15 (nokta(2,2)-110) (nokta(6,1)-nokta(3,1)+85) (nokta(5,2)-nokta(1,2)+130)];
            case {7, 8}
                dikd = [nokta(3,1)+3 (nokta(2,2)-25) (nokta(6,1)-nokta(3,1)-28) (nokta(5,2)-nokta(1,2)+11)];
                dikd2 = [nokta(3,1)-10 (nokta(2,2)-140) (nokta(6,1)-nokta(3,1)+85) (nokta(5,2)-nokta(1,2)+160)];
            case 9
                dikd = [(nokta(1,1)-17) (nokta(3,2)-25) (nokta(7,1)-nokta(1,1)+1) (nokta(6,2)-nokta(2,2)+11)];
                dikd2 = [(nokta(1,1)-30) (nokta(3,2)-120) (nokta(7,1)-nokta(1,1)+80) (nokta(6,2)-nokta(2,2)+130)];
            case 10
                dikd = [(nokta(1,1)-25) (nokta(3,2)-25) (nokta(7,1)-nokta(1,1)+18) (nokta(6,2)-nokta(2,2)+11)];
                dikd2 = [(nokta(1,1)-50) (nokta(3,2)-135) (nokta(7,1)-nokta(1,1)+90) (nokta(6,2)-nokta(2,2)+150)];
            case 11
                dikd = [nokta(6,1)+20 (nokta(1,2)-25) (nokta(3,1)-nokta(6,1)-10) (nokta(5,2)-nokta(2,2)+11)];
                dikd2 = [nokta(6,1)-60 (nokta(1,2)-120) (nokta(3,1)-nokta(6,1)+90) (nokta(5,2)-nokta(2,2)+140)];
            case {12, 13}
                dikd = [nokta(6,1)+30 (nokta(1,2)-25) (nokta(3,1)-nokta(6,1)-30) (nokta(5,2)-nokta(2,2)+11)];
                dikd2 = [nokta(6,1)-60 (nokta(1,2)-130) (nokta(3,1)-nokta(6,1)+70) (nokta(5,2)-nokta(2,2)+150)];
        end
    catch e
        rethrow(e); % carry error, rethrow it until main function
    end
    
    % dikd: small rectangle, dikd2: large rectangle
    basicR = round(0.4*dikd+0.6*dikd2);
end
