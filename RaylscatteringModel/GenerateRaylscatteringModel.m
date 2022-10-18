%--------------------------------------------------------------------------
% * Kaichao Liang, 2022.01.26
% * Generate 2D phantom X-ray diffraction matrix PanelData= DiffSys*Object
% online. Only calculate object without rotation considering the matrix
% scale.
% * Target domain: [Energy，DetZ，DetY], source domain: [PixY, PixX].
%--------------------------------------------------------------------------
function DiffSys = GenerateRaylscatteringModel(SourcePos,PixX,PixY,PixelSize, ...
    DetPos,DetY,DetZ,offsetY,offsetZ,DetSize,DetMask,SpeSet,MTSet)
    %%-----------------------------paramter--------------------------------
    % objectXRD: the Raylscattering distribution, the size [numel(MTSet),PixY, PixX]
    % SourcePos: the distance from source to iso-center.
    % PixX: the number of object pixels in X direction.
    % PixY: the number of object pixels in Y direction.
    % PixelSize: object discrete pixel size.
    % DetPos: the distance from iso-center to detector plane.
    % DetY: the number of detector bins in Y direction.
    % DetZ: the number of detector bins in Z direction.
    % offSetY, offSetZ: off-center placement of detecotr.
    % DetSize: the size of detector bin.
    % DetMask: a mask used detector bin, if 0 ignore the detector bin.
    % SpeSet: the discrete vector for spectrum.
    % MTSet: the discrete vector for momentum transfer.
    %%---------------------------------------------------------------------
    tic;
    %%discrete phantom, detector, and coded aperture
    PosPixX = [-(PixX-1)/2:(PixX-1)/2]*PixelSize;
    PosPixY = [-(PixY-1)/2:(PixY-1)/2]*PixelSize;
    
    PosDetY = [-(DetY-1)/2:(DetY-1)/2]*DetSize+offsetY;
    PosDetZ = [-(DetZ-1)/2:(DetZ-1)/2]*DetSize+offsetZ;
    
    discreteNumDet = 5;
    discreteNumE = 5;
    
    DiffSys = cell(PixY,PixX);
    
    for x = 1:PixX %%Loop for pixel
        
        posx = PosPixX(x);
        
        for y = 1:PixY
            
            DiffSysY = sparse([]);
            fprintf('Calculate pix X %d, Y %d\n',x,y);
            posy = PosPixY(y);
            
            for dety = 1:DetY %%Loop for detector
                
                DiffSysDetY = sparse([]);
                posdety = PosDetY(dety);
                discreteDetY = posdety + [-discreteNumDet/2+0.5:discreteNumDet/2-0.5]*DetSize/discreteNumDet;
                discreteDetY = repmat(reshape(discreteDetY,1,1,discreteNumDet),discreteNumE,discreteNumDet,1);
                discreteDetY = reshape(discreteDetY,discreteNumE*discreteNumDet*discreteNumDet,1);
                
                for detz = 1:DetZ/2
                    
                    DiffSysDetZ = sparse(zeros(numel(SpeSet),numel(MTSet)));
                    if(~DetMask(detz,dety))
                        DiffSysDetY = [DiffSysDetY;DiffSysDetZ];
                        continue;
                    end
                    
                    posdetz = PosDetZ(detz);
                    discreteDetZ = posdetz + [-discreteNumDet/2+0.5:discreteNumDet/2-0.5]*DetSize/discreteNumDet;
                    discreteDetZ = repmat(reshape(discreteDetZ,1,discreteNumDet,1),discreteNumE,1,discreteNumDet);
                    discreteDetZ = reshape(discreteDetZ,discreteNumE*discreteNumDet*discreteNumDet,1);
                    
                    inScalar = [posx+SourcePos;posy;0];
                    outScalar = [DetPos-posx;posdety-posy;posdetz];
                    cosAngle=sum(inScalar.*outScalar)/norm(inScalar)/norm(outScalar);
                    for e  = 1:numel(SpeSet) %%Loop for each energy, the detector size and detector energy determine momenturm transfer
                        
                        discreteDetE = SpeSet(e) + [-discreteNumE/2+0.5:discreteNumE/2-0.5]*(SpeSet(2)-SpeSet(1))/discreteNumDet;
                        discreteDetE = repmat(reshape(discreteDetE,discreteNumE,1,1),1,discreteNumDet,discreteNumDet);
                        discreteDetE = reshape(discreteDetE,discreteNumE*discreteNumDet*discreteNumDet,1);
                        
                        %q=E*sin(theta/2)/hv
                        inScalar = [posx+SourcePos;posy;0];
                        outScalar = [DetPos*ones(discreteNumE*discreteNumDet*discreteNumDet,1)-posx,discreteDetY-posy,discreteDetZ];
                        scatterAngle = acos(outScalar*inScalar/norm(inScalar)./sqrt(outScalar(:,1).^2+outScalar(:,2).^2+outScalar(:,3).^2));
                        momentum = discreteDetE.*sin(scatterAngle/2)/1.24;
                        
                        %form the q-E sparse matrix
                        momentumIndex = floor((momentum-MTSet(1))./(MTSet(2)-MTSet(1))+1.5);
                        momentumIndex(find(momentumIndex<=0))=1;
                        momentumIndex(find(momentumIndex>=numel(MTSet)))=numel(MTSet);
                        
                        Counts = tabulate(momentumIndex);
                        weights = (1+cosAngle^2)/2*(DetPos-posx)^2/((DetPos-posx)^2+(posdety-posy)^2+posdetz^2)*(DetPos-posx)/sqrt((DetPos-posx)^2+(posdety-posy)^2)...
                            *(DetPos-posx)/sqrt((DetPos-posx)^2+posdetz^2);
                        DiffSysDetZ(e,Counts(:,1))=Counts(:,2)/discreteNumE/discreteNumDet/discreteNumDet*(DetPos/(DetPos-posx))^2*(SourcePos/(SourcePos+posx))*weights;
                    end
                    DiffSysDetY = [DiffSysDetY;DiffSysDetZ]; 
                end
                DiffSysY = [DiffSysY;DiffSysDetY];
            end
            DiffSys{y,x} = DiffSysY;  
        end
        %save DiffSys DiffSys
    end
    toc;
end

    

    
