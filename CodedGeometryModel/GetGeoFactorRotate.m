%--------------------------------------------------------------------------
% * Kaichao Liang, 2022.05.15
% * Coded aperture geometry response for each pixel based on solid angle
% intergration.
% * Update 2022.05.15: Coded-aperture thickness model
%--------------------------------------------------------------------------
function geoFactor =  GetGeoFactorRotate(curPosX, curPosY, pixelSize, BeamWidth,discretePixNum, keyPointX, keyPointY, CodePos,CodeThick, DetPos, rotAngle)
    %%-----------------------------paramter--------------------------------
    % curPosX: the X position of pixel center.
    % curPosY: the Y position of pixel center.
    % pixelSize: object discrete pixel size.
    % discretePixNum: the number of sub pixels for intergration.
    % CodePos: the X position of coded aperture.
    % DetPos: the X position of detector.
    % keyPoint X,Y: the location of a code hole and a detector.
    % rotAngle: object rotation angle.
    %%---------------------------------------------------------------------
    discreteU = [-(discretePixNum-1)/2:(discretePixNum-1)/2]*pixelSize/discretePixNum;
    discreteV = [-(discretePixNum-1)/2:(discretePixNum-1)/2]*pixelSize/discretePixNum;
    discreteZ = [-(discretePixNum-1)/2:(discretePixNum-1)/2]*BeamWidth/discretePixNum;
    discreteX = discreteU'*cos(rotAngle) -  discreteV*sin(rotAngle) + curPosX;
    discreteY = discreteU'*sin(rotAngle) +  discreteV*cos(rotAngle) + curPosY;
    
    discreteZ=repmat(reshape(discreteZ,discretePixNum,1,1),1,discretePixNum,discretePixNum);
    discreteX=repmat(reshape(discreteX,1,discretePixNum,discretePixNum),discretePixNum,1,1);
    discreteY=repmat(reshape(discreteY,1,discretePixNum,discretePixNum),discretePixNum,1,1);
    
    %-------------------------------------
    % Factor in X direction
    %-------------------------------------
    upSect1 = (keyPointX(2)-discreteZ)./(CodePos-CodeThick/2-discreteX)*(DetPos-CodePos+CodeThick/2)+keyPointX(2);
    upSect2 = (keyPointX(2)-discreteZ)./(CodePos+CodeThick/2-discreteX)*(DetPos-CodePos-CodeThick/2)+keyPointX(2);
    upSect = min(keyPointX(4),max(keyPointX(3),min(upSect1,upSect2)));
    
    lowSect1 = (keyPointX(1)-discreteZ)./(CodePos-CodeThick/2-discreteX)*(DetPos-CodePos+CodeThick/2)+keyPointX(1);
    lowSect2 = (keyPointX(1)-discreteZ)./(CodePos+CodeThick/2-discreteX)*(DetPos-CodePos-CodeThick/2)+keyPointX(1);
    lowSect = min(keyPointX(4),max(keyPointX(3),max(lowSect1,lowSect2)));
    
    factorX = mean(upSect(:)-lowSect(:));
    
    if(factorX==0)
        geoFactor=0;
        return;
    end
    
    %-------------------------------------
    % Factor in Y direction
    %-------------------------------------
    rightSect1 = (keyPointY(2)-discreteY)./(CodePos-CodeThick/2-discreteX)*(DetPos-CodePos+CodeThick/2)+keyPointY(2);
    rightSect2 = (keyPointY(2)-discreteY)./(CodePos+CodeThick/2-discreteX)*(DetPos-CodePos-CodeThick/2)+keyPointY(2);
    rightSect = min(keyPointY(4),max(keyPointY(3),min(rightSect1,rightSect2)));
    
    leftSect1 = (keyPointY(1)-discreteY)./(CodePos-CodeThick/2-discreteX)*(DetPos-CodePos+CodeThick/2)+keyPointY(1);
    leftSect2 = (keyPointY(1)-discreteY)./(CodePos+CodeThick/2-discreteX)*(DetPos-CodePos-CodeThick/2)+keyPointY(1);
    leftSect = min(keyPointY(4),max(keyPointY(3),max(leftSect1,leftSect2)));
    factorY = mean(rightSect(:)-leftSect(:));
    
    geoFactor = factorX*factorY/(keyPointY(4)-keyPointY(3))^2;
    
end
