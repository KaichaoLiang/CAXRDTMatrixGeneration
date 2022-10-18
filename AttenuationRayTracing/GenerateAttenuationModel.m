%--------------------------------------------------------------------------
% * Kaichao Liang, 2022.01.10
% * Attenuation line-intergral model for XRD based on ray-tracing model.
% * Given source domain attenuation distribution [PixY, PixX],calculate the
% scatter signal attenuation weighting.
% * Target domain: [DetY，PixY，PixX], source domain: [PixY, PixX].
%--------------------------------------------------------------------------

function Hsys = GenerateAttenuationModel(SourcePos,PixX,PixY,PixelSize,DetPos,DetY,offsetY,DetSize,RotAngle)
    tic;
    %%-----------------------------paramter--------------------------------
    % SourcePos: the distance from source to iso-center.
    % PixX: the number of object pixels in X direction.
    % PixY: the number of object pixels in Y direction.
    % PixelSize: object discrete pixel size;
    % DetPos: the distance from iso-center to detector plane.
    % DetY: the number of detector bins in Y direction;
    % offSetY: off-center placement of detecotr.
    % DetSize: the size of detector bin;
    % RotAngle: the object rotate angle, [0, 2*pi)
    %%---------------------------------------------------------------------
    
    %%Normalize all lengths by pixelsize.
    SourcePos = SourcePos/PixelSize;
    DetPos = DetPos/PixelSize;
    offsetY = offsetY/PixelSize;
    DetSize = DetSize/PixelSize;
    
    %%Generate grid.
    grid.x=linspace(-PixX/2,PixX/2,PixX+1); % Nx+1表示，这些都是点数，Nx+1个点
    grid.y=linspace(-PixY/2,PixY/2,PixY+1);
    
    PosPixX = [-(PixX-1)/2:(PixX-1)/2];
    PosPixY = [-(PixY-1)/2:(PixY-1)/2];
    
    %%X-ray source and detector position.
    srcPosition0=[-SourcePos;0];
    detPosition0=[DetPos*ones(1,DetY);[-(DetY-1)/2:(DetY-1)/2]*DetSize+offsetY];
    
    %%Source and detector clockwise rotation.
    szb.x=srcPosition0(1)*cos(RotAngle)+srcPosition0(2)*sin(RotAngle);%源坐标
    szb.y=-srcPosition0(1)*sin(RotAngle)+srcPosition0(2)*cos(RotAngle);
            
    Hsys = sparse(0);
    
    for x = 1:PixX %Scatter point at pixx
        HPixX =sparse(0);
        midX = PosPixX(x);
        for y =1:PixY %Scatter point at pixx
            HPixY =sparse(0);
            midY = PosPixY(y);
            mzb.x = midX;
            mzb.y = midY;
            sectsIn = get_sect(szb,mzb,grid);%the sects along incident path
            
            %%clean the sects out of szb,mzb range.
            cosvalue = (mzb.x-sectsIn(1,:))*(mzb.x-szb.x) + (mzb.y-sectsIn(2,:))*(mzb.y-szb.y);
            abandonIndex = find(cosvalue<0);
            [maxCos,maxIndex] = max(cosvalue(abandonIndex));
            abandonIndex(maxIndex)=[]; %%Keep the attenuation of current pix;
            sectsIn(:,abandonIndex)=[];
            
            if size(sectsIn,2)>=2
               SMIn = get_li_xyx(sectsIn,grid);
            else
               SMIn = [];
            end
            
            for dety=1:DetY % Detector loop
                fprintf('calculate pix,: %d, %d, code %d\n',x, y, dety);
                HdetTmp = sparse(zeros(1, PixY*PixX));
                xx0 = detPosition0(1,dety);
                yy0 = detPosition0(2,dety);
                dzb.x= xx0*cos(RotAngle)+yy0*sin(RotAngle);%探测器坐标
                dzb.y= yy0*cos(RotAngle)-xx0*sin(RotAngle);
                sectsOut = get_sect(mzb,dzb,grid);%the sects along scatter path
                
                %%clean the sects out of szb,mzb range.
                cosvalue = (mzb.x-sectsOut(1,:))*(mzb.x-dzb.x) + (mzb.y-sectsOut(2,:))*(mzb.y-dzb.y);
                sectsOut(:,find(cosvalue<0))=[];
                
                if size(sectsOut,2)>=2
                    SMOut = get_li_xyx(sectsOut,grid);
                else
                    SMOut = [];
                end
                SM =cat(2,SMIn,SMOut);
                if(size(SM,2)>0)
                    HdetTmp(SM(1,:))=SM(2,:)*PixelSize;
                end
                if(size(HPixY,2)==1)
                    HPixY = HdetTmp;
                else
                    HPixY = [HPixY;HdetTmp];
                end
            end
            if(size(HPixX,2)==1)
                HPixX = HPixY;
            else
                HPixX = [HPixX;HPixY];
            end    
        end
        if(size(Hsys,2)==1)
            Hsys = HPixX;
        else
            Hsys = [Hsys;HPixX];
        end 
        
    end
    toc;
end
    