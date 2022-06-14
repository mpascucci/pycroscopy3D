function skew_correct_one(TiffFolder, SaveFolder, RunName, InfoFile)
    %Open a dataset with volume Avg
    

    TifName=strcat(TiffFolder,'/', RunName);
    
    
    %On ouvre les fichiers .tiff
    Infos=imfinfo(TifName);
    NLines=Infos.Height;
    NCol=Infos.Width;
    NPlanes=size(Infos,1);
    
     
    Ima=zeros(NLines,NCol, NPlanes);
    
    for zz=1: NPlanes             
        Ima(:,:,zz)=imread(TifName,zz,'info',Infos);
    end
    %On 
    
    % load(strcat(TiffFolder,'/', InfoFile));
    load(InfoFile);
       
    ylat = info.GUIcalFactors.y_umPerPix;
    zdep = info.GUIcalFactors.z_umPerPix;
    xfact = info.GUIcalFactors.xK_umPerVolt;
    xwid = xfact*info.daq.scanAngle/info.daq.pixelsPerLine;
    
    conversionFactors = [ylat; zdep; xwid];
    
    SkewAng = -47;
    temp = Ima;
    delta = SkewAng;
    affineMatrix = [1 0 0 0;
                    0 1 0 0;
                    0 cotd(delta) 1 0;
                    0 0 0 1];
    tform = affine3d(affineMatrix);
    temp = flip(flip(temp,2),3);
    temp = permute(temp, [3 1 2]);
    R = imref3d(size(temp), conversionFactors(1), conversionFactors(3), conversionFactors(2));
    [temp, ~] = imwarp(temp, R, tform);
    
    SaveName=[SaveFolder];
    mkdir(SaveName);
    % timeStepName = strsplit(RunName, '_');
    % timeName = char(strcat(timeStepName(2), '_', timeStepName(3), '_', timeStepName(4)));
    % imgToSave2=[SaveName, '/', timeName, '.tif'];
    [RunNamePath,RunNameFileName,RunNameext] = fileparts(RunName)
    imgToSave2=[SaveName, '/', strcat(RunNameFileName,'_skewCorrected'), '.tif'];
    imwrite(uint16(temp(:, :, 1)), imgToSave2, 'tif', 'Compression', 'none');                
    for j = 2:size(temp, 3)
        imwrite(uint16(temp(:, :, j)), imgToSave2, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
    end
            
end
