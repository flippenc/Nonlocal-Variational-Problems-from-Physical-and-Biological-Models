workingDir = "D:";
imageNames = dir(fullfile(workingDir,"videoFiles","*.png"));
imageNames = {imageNames.name}';
outputVideo = VideoWriter(fullfile(workingDir,"lambda1.3.avi"));
outputVideo.FrameRate = 4;
open(outputVideo)
for i = 1:length(imageNames)
   img = imread(fullfile(workingDir,"videoFiles",imageNames{i}));
   writeVideo(outputVideo,img)
end
close(outputVideo)