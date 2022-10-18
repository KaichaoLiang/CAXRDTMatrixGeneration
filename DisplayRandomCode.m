%load GridPosY;
%load GridPosZ;

load CodedOpenFull;

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

CodedOpen=CodedOpenFull;

GridPosY = FullGridPosY(find(CodedOpen==1));
%GridPosY=(rand(1800,1)-0.5)*45;
%GridPosZ=(rand(1800,1)-0.5)*45;
GridPosZ = FullGridPosZ(find(CodedOpen==1));
load GridPosYUn; load GridPosZUn;
randomDisplay = zeros(900,900);
numCode = numel(GridPosY);
for i =1:numCode
    indexY = round(GridPosY(i)*10+450);
    indexZ = round(GridPosZ(i)*10+450);
    randomDisplay(max(1,indexZ-4):min(900,indexZ+4),max(1,indexY-4):min(900,indexY+4))=1;
end
imshow(randomDisplay);