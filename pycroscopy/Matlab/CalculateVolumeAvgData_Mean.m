TiffFolder='/network/lustre/iss01/wyart/rawdata/mathilde.lapoix/SCAPE/tiff_stacks/220210/F2_run5subset_for_volumeavg/';
SaveFolder='/network/lustre/iss01/wyart/analyses/mathilde.lapoix/SCAPE/Illustrations/Volume_Avg/';
cd(TiffFolder)

%Find all tif files in the initial directory.
list=dir('*/pla*.tiff');
NPlanes=size(list,1);

%Ouverture des images

cd(list(1).folder)
TifName=list(1).name;
Infos=imfinfo(TifName);

NLines=Infos.Height;
NCol=Infos.Width;
NIma=size(Infos,1);

Avg=zeros(NLines,NCol,NPlanes);
Med=zeros(NLines,NCol,NPlanes);

 
for nn=1:NPlanes %For each plane...
nn

 cd(list(nn).folder)
 TifName=list(nn).name;
 Infos=imfinfo(TifName);

  %...Read timelapse   
  Ima=zeros(NLines,NCol,NIma);
    for tt=1:NIma              
        Ima(:,:,tt)=imread(TifName,tt,'info',Infos);
    end
  
    Avg(:,:,nn)=mean(Ima,3);
    Med(:,:,nn)=median(Ima,3);
end

Avg=uint16(rescale(Avg,0,2^16-1));
Med=uint16(rescale(Med,0,2^16-1));

for nn=1:NPlanes
    if nn==1
        imwrite(Avg(:,:,nn),strcat(SaveFolder,'Volume_Avg.tif'));
        % imwrite(uint16(SumDFF),strcat(SaveFolder,'Volume_DFF.tif'));
        imwrite(Med(:,:,nn),strcat(SaveFolder,'Volume_Med.tif'));
        
        
    else
        imwrite(Avg(:,:,nn),strcat(SaveFolder,'Volume_Avg.tif'),'writeMode','append');
        %imwrite(uint16(SumDFF),strcat(SaveFolder,'Volume_DFF.tif'),'writeMode','append');
         imwrite(uint16(Med(:,:,nn)),strcat(SaveFolder,'Volume_Med.tif'),'writeMode','append');
    end

end
