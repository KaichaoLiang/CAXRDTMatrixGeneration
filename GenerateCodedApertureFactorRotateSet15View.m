%--------------------------------------------------------------------------
% * Kaichao Liang, 2021.1.11
% * Generate coded aperture XRD imaging system matrix sub factors for different angle
% * Geometry system matrix: Target domain: [DetZ, DetY], source domain: [PixY, PixX]
% * Attenuation system matrix: Target domain: [DetY,PixY,PixX], source domain: [PixY, PixX]
%--------------------------------------------------------------------------
dbstop if error;
addpath('CodedGeometryModel');
addpath('AttenuationRayTracing');
addpath('RaylscatteringModel');
%p=parpool(54);

numAngle = 15;
dAngle = 2*pi/numAngle;
AngleSet = [1:numAngle]*dAngle;

SourcePos = 400;

%物体参数
PixX = 54;
PixY = 54;
PixelSize =1;
BeamWidth = 0.5;

%探测器参数
DetPos = 297.5;
DetY = 64;
DetZ = 64;
offsetY = 0;
offsetZ = 0;
DetSize = 1.6;

SpeSet = 21:85;
MTSet= 0.01:0.02:4;

%编码孔离散模型参数
CodePos = 150;
GridIntervalY = 1.5;
FullGridNumY = 60;
GridIntervalZ = 1.5;
FullGridNumZ = 60;
CodeSize =1;
CodeThick = 1.5;

FullGridPosY=repmat([-(FullGridNumY-1)/2:(FullGridNumY-1)/2]'*GridIntervalY,1,FullGridNumZ);
FullGridPosZ=repmat([-(FullGridNumZ-1)/2:(FullGridNumZ-1)/2]*GridIntervalZ,FullGridNumY,1);


%load CodedOpenFull;
%CodedOpen=CodedOpenFull;
DetMask = ones(DetZ,DetY);
DetMask(31:34,:)=0;
%GridPosY = FullGridPosY(find(CodedOpen==1));
%GridPosZ = FullGridPosZ(find(CodedOpen==1));
load GridPosYUn;
load GridPosZUn;
GridNum=numel(GridPosY);
GridPosY = GridPosY+0.5;

%mkdir('RaylMatrix');
%DiffSys = GenerateRaylscatteringModel(SourcePos,PixX,PixY,PixelSize,DetPos,DetY,DetZ,offsetY,offsetZ,DetSize,DetMask,SpeSet,MTSet);
%save('RaylMatrix/DiffSys4W.mat','DiffSys');


mkdir('DecayMatrix15ViewOffY');
mkdir('CodeGeoMatrix15ViewOffY');

for angle = 1:numAngle %loop for rotation angle
    DecaySys = GenerateAttenuationModel(SourcePos,PixX,PixY,PixelSize, ...
    DetPos,DetY,offsetY,DetSize,AngleSet(angle));
         
    GeoSys = GenerateCodedGeometryRotate(PixX,PixY,PixelSize,BeamWidth,CodePos,CodeThick,GridPosY,GridPosZ,GridNum, ...
    DetPos,DetY,DetZ,offsetY,offsetZ,DetSize,DetMask,CodeSize,AngleSet(angle));
    
    saveDecay = sprintf("DecaySysAngle%d = DecaySys;",angle);eval(saveDecay);
    saveDecayName = sprintf("DecayMatrix15ViewOffY/DecaySysAngle%d.mat",angle);
    decayVarName = sprintf("DecaySysAngle%d",angle);
    save(saveDecayName, decayVarName);

    saveGeo = sprintf("GeoSysAngle%d = GeoSys;", angle);eval(saveGeo);
    saveGeoName = sprintf("CodeGeoMatrix15ViewOffY/GeoSysAngle%d.mat",angle);
    geoVarName = sprintf("GeoSysAngle%d",angle);
    save(saveGeoName,geoVarName);
    
end

