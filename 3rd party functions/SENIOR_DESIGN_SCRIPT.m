% function record_array_to_wav()
% in this example code we demonstrate how to capture a set of data from the
% array and save that data as wav files. This demonstrates how to
% create a Microphone_Array object and stream data into the matlab
% workspace.

function SENIOR_DESIGN_SCRIPT(measurement_angle, nBlocksToGrab,filepath)

block_size = 8192; %this will define the number of samples returned per channel
% nBlocksToGrab = 25; %this will define how many consecutive blocks we wish to record. 


% **** first we create a microphone array object

MA = Microphone_Array();

% **** the next step is to set any parameters, in this case we will set the 
%      block transfer size to 8192 

MA.block_size = block_size;

% **** set the gain setting for the internal analog circuitry, 
%  1 =  -15dB
%  2 =  0dB
%  3 =  15dB
MA.gainSetting = 2;

% **** No we can initialize the array, this will attempt to open the driver
% and connect to the microphone array.  If the array is not connected to
% the computer then this will fail. 

MA.init_array();

% Now we must start the internal acquisition of the microphone array data.
% The driver will begin acquiring data into the background circular buffer
% immediately after this call. 

MA.start();

%At this point we can begin grabbing data from the array.  In order to grab
%data from the array we must first creat a buffer in Matlab to store all
%audio data.  This array should be of size n_channels x block_size

data = zeros(MA.n_channels*MA.block_size,1);

% And we will create a large matrix to store all of the audio data acquired

totalData = zeros(MA.block_size*nBlocksToGrab,MA.n_channels);

% cleanupObj = onCleanup(@() closeVDAM(measurement_angle, totalData,filepath,MA));

% Now we will grab <nBlocksToGrab> consecutive frames 

for i = 1:nBlocksToGrab
    fprintf('Grabbing frame %d\n',i);
    data = getDataBlock(MA,data);
    totalData((i-1)*MA.block_size+1:i*MA.block_size,:) = reshape(data,[MA.block_size,MA.n_channels]);
end

% 
% % now we can write these to n_channel individual wav files
% plot(totalData(:,1:8))
% for i = 1:MA.n_channels
%     fprintf('Writing file number %d\n',i);
%     filename = sprintf('%s_%0.2d.wav',wavPrefix,i);
%     audiowrite(filename, totalData(:,i),MA.sample_rate);
% %     wavwrite(totalData(:,i),MA.sample_rate,filename);
% end

%now we need to clean up the workspace by stopping the acquisition and then
%closing the Microphone_Array object


closeVDAM(measurement_angle,totalData,filepath,MA);


end