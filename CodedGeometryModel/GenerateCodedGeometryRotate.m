%--------------------------------------------------------------------------
% * Kaichao Liang, 2022.05.15
% * Grid coded aperture XRD imaging system matrix, linear
% transformation.Detector intersection weighting intergral method. 
% * Target domain: [DetZ, DetY], source domain: [PixY, PixX]
% * Update 2022.05.15 Coded aperture thickness model
%--------------------------------------------------------------------------
function CodedSysGeoMatrix = GenerateCodedGeometryRotate(PixX,PixY,PixelSize,BeamWidth,CodePos,CodeThick,GridPosSetY,GridPosSetZ,GridNum, ...
    DetPos,DetY,DetZ,offsetY,offsetZ,DetSize,DetMask,CodeSize,RotAngle)
    tic;
    %%-----------------------------paramter--------------------------------
    % PixX: the number of object pixels in X direction.
    % PixY: the number of object pixels in Y direction.
    % PixelSize: object discrete pixel size.
    % CodePos: the distance from iso-center to coded aperture plane.
    % GridPosSetY: the vector of code Y positions.
    % GridPosSetZ: the vector of code Z positions.
    % GridNum: total number of opened codes.
    % DetPos: the distance from iso-center to detector plane.
    % DetY: the number of detector bins in Y direction.
    % DetZ: the number of detector bins in Z direction.
    % offSetY, offSetZ: off-center placement of detecotr.
    % DetSize: the size of detector bin.
    % DetMask: a mask used detector bin, if 0 ignore the detector bin.
    % CodeSize: the aperture open size.
    % RotAngle: the object rotate angle, [0, 2*pi)
    %%---------------------------------------------------------------------
    %%cell store structure for system matrix
    CodedSysGeoMatrix = zeros(DetZ,DetY,PixY,PixX);
    
    %%discrete phantom, detector, and coded aperture
    PosPixX = [-(PixX-1)/2:(PixX-1)/2]*PixelSize;
    PosPixY = [-(PixY-1)/2:(PixY-1)/2]*PixelSize;
    
    PosDetY = [-(DetY-1)/2:(DetY-1)/2]*DetSize+offsetY;
    PosDetZ = [-(DetZ-1)/2:(DetZ-1)/2]*DetSize+offsetZ;
    PosGridY = GridPosSetY;
    PosGridZ = GridPosSetZ;
    
    discretePixNum = 15;
    %%Generate system matrix loop by detector bin
    for dety = 1:DetY
        parfor detz = 1:DetZ    
            %%Unused pixel
            if(~DetMask(detz,dety))
                continue;
            end
            
            %%Each code loop
            for code = 1:GridNum 
                fprintf('calculate pix,: %d, %d, code %d\n',dety, detz, code);
                
                %%Pre-remove useless codes
                if(abs(PosDetZ(detz)-PosGridZ(code))<abs(PosGridZ(code))/4)
                    continue;
                end
                if((PosDetZ(detz)-PosGridZ(code))*PosGridZ(code)<0)
                    continue;
                end
                
                %-------------------------------------
                % Cauculate effective pixel range in X direction
                %-------------------------------------
                posx1 = CodePos-(DetPos-CodePos)/(PosDetZ(detz)-DetSize/2-PosGridZ(code)-CodeSize/2)*(PosGridZ(code)+CodeSize/2);
                d1 = CodePos-(DetPos-CodePos)/(PosDetZ(detz)-DetSize/2-PosGridZ(code)-CodeSize/2)*(PosGridZ(code)+CodeSize/2+BeamWidth/2);
                d2 = CodePos-(DetPos-CodePos)/(PosDetZ(detz)-DetSize/2-PosGridZ(code)-CodeSize/2)*(PosGridZ(code)+CodeSize/2-BeamWidth/2);
                posx2 = CodePos-(DetPos-CodePos)/(PosDetZ(detz)+DetSize/2-PosGridZ(code)+CodeSize/2)*(PosGridZ(code)-CodeSize/2);
                d3 = CodePos-(DetPos-CodePos)/(PosDetZ(detz)+DetSize/2-PosGridZ(code)+CodeSize/2)*(PosGridZ(code)-CodeSize/2+BeamWidth/2);
                d4 = CodePos-(DetPos-CodePos)/(PosDetZ(detz)+DetSize/2-PosGridZ(code)+CodeSize/2)*(PosGridZ(code)-CodeSize/2-BeamWidth/2);
                
                if(posx1>posx2)
                    tmp = d1;
                    d1=d3;
                    d3=tmp;
                    tmp = d2;
                    d2=d4;
                    d4=tmp;
                end
                posx1 = min(d1,d2);
                posx2 = max(d3,d4);
                
                if(posx2<PosPixX(1)-PixelSize/2 || posx1>PosPixX(end)+PixelSize/2)
                    continue;
                end
                
                
                for x = 1:PixX
                    curPosU=PosPixX(x);
                    
                    for y = 1:PixY
                        curPosV=PosPixY(y);
                        
                        curPosX = curPosU*cos(RotAngle)-curPosV*sin(RotAngle);
                        curPosY = curPosU*sin(RotAngle)+curPosV*cos(RotAngle);
                        
                        if(curPosX<posx1-PixelSize*sqrt(2)/2 || curPosX>posx2+PixelSize*sqrt(2)/2)
                            continue;
                        end
                    
                        %-------------------------------------
                        % Cauculate effective pixel range in Y direction
                        %-------------------------------------
                    
                        posy1 = (PosGridY(code)+CodeSize/2)-(PosDetY(dety)-DetSize/2-PosGridY(code)-CodeSize/2)/(DetPos-CodePos)*(CodePos-curPosX);
                        posy2 = (PosGridY(code)-CodeSize/2)-(PosDetY(dety)+DetSize/2-PosGridY(code)+CodeSize/2)/(DetPos-CodePos)*(CodePos-curPosX);
                
                        if(posy1>posy2)
                            tmp=posy1;
                            posy1=posy2;
                            posy2=tmp;
                        end
                
                        if(posy2<PosPixY(1)-PixelSize/2 || posy1>PosPixY(end)+PixelSize/2)
                            continue;
                        end
                    
                        if(curPosY<posy1-PixelSize*sqrt(2)/2 || curPosY>posy2+PixelSize*sqrt(2)/2)
                            continue;
                        end
                        
                        %-------------------------------------
                        % Discrete intergral geometry factor
                        %-------------------------------------
                        keyPointX = [PosGridZ(code)-CodeSize/2,PosGridZ(code)+CodeSize/2,PosDetZ(detz)-DetSize/2,PosDetZ(detz)+DetSize/2];
                        keyPointY = [PosGridY(code)-CodeSize/2,PosGridY(code)+CodeSize/2,PosDetY(dety)-DetSize/2,PosDetY(dety)+DetSize/2];
                        geoFactor =  GetGeoFactorRotate(curPosX, curPosY, PixelSize, BeamWidth, discretePixNum, keyPointX, keyPointY, CodePos, CodeThick,DetPos, RotAngle);  
 
                        CodedSysGeoMatrix(detz,dety,y,x) = CodedSysGeoMatrix(detz,dety,y,x)+geoFactor;
                    end
                end                                   
            end 
        end
    end
    CodedSysGeoMatrix = reshape(CodedSysGeoMatrix, DetZ*DetY,PixY*PixX);
    toc;
end

              
                
                
                    
                   
                
                
    
    
    
    
    
    