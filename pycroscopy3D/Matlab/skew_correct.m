function skew_correct(TiffFolder, SaveFolder, InfoFile)
    listing = dir(strcat(TiffFolder,'/*.tif'));
    parfor j = 1:length(listing)
        %disp(listing(j).name);
        skew_correct_one(TiffFolder, SaveFolder, listing(j).name, InfoFile);
    end
end