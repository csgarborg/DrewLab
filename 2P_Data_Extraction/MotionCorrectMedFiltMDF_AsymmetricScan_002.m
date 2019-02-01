function [proctime,transfertime]=MotionCorrectMedFiltMDF_AsymmetricScan_002(fname,Imaging_channels,Analog_Channel_Num,temporal_filter_size,ChunkSize,PixelOffset,varargin)
if exist('fname','var')==0
    fname=uigetfile('*.MDF','MultiSelect','off'); %make list of MDF files to be filtered and image registered
end
if exist('Imaging_channels','var')==0
    Imaging_channels=input('How many channels were imaged?');
end
if exist('temporal_filter_size','var')==0
    temporal_filter_size=input('How many frames for median filter?');
end
if exist('ChunkSize','var')==0
    ChunkSize=5000; %number of frames to be exported to memory at same time
end
if exist('Analog_Channel_Num','var')==0
    Analog_Channel_Num=1;
end

if exist('PixelOffset','var')==0
     PixelOffset=0;
end
thepath=cd;
TheFile=[thepath '\' fname];
filebase=fname(1:(end-4));
%ImagingData=[];
save([filebase '_MocoMedfilt.mat'],'-v7.3');
TheImagingData=matfile([filebase '_MocoMedfilt.mat'],'Writable',true);%File that will have imaging data added to
f=figure('Position',[300,300,500,500]); %Arbitrary dimentions used for ActiveX controller
mfile=actxcontrol('MCSX.Data',[0,0,500,500],f);%initialize controller to read in frames
openResult=invoke(mfile,'OpenMCSFile',TheFile);
if openResult ~=0
    fprintf('MDF file open fail!\n');
    pause;
end
%Get imaging parameters
fprintf('Getting imaging parameters\n')
mv_mpP(1).Header.Filename=fname;
mv_mpP(1).Header.Frame_Width = str2double(invoke(mfile,'ReadParameter','Frame Width'));
mv_mpP(1).Header.Frame_Height = str2double(invoke(mfile,'ReadParameter','Frame Height'));
mv_mpP(1).Header.Frame_Duration=invoke(mfile,'ReadParameter','Frame Duration (s)');
mv_mpP(1).Header.Frame_Rate=1/str2double(mv_mpP(1).Header.Frame_Duration(1:(end-1)));
mv_mpP(1).Header.Magnification= invoke(mfile,'ReadParameter','Magnification');
mv_mpP(1).Header.Rotation=invoke(mfile,'ReadParameter','Rotation');
mv_mpP(1).Header.num_frames= str2double(invoke(mfile,'ReadParameter','Frame Count'));
mv_mpP(1).Header.Pixel_Clock=invoke(mfile,'ReadParameter','Pixel Clock');
mv_mpP(1).Header.Microns_Per_Pixel=invoke(mfile,'ReadParameter','Microns per Pixel');
mv_mpP(1).xsize = (mv_mpP(1).Header.Frame_Width);
mv_mpP(1).ysize = (mv_mpP(1).Header.Frame_Height);
mv_mpP(1).rowStart=53;
mv_mpP(1).rowEnd=460;
mv_mpP(1).startframe=1;
mv_mpP(1).endframe=mv_mpP(1).Header.num_frames;
mv_mpP(1).Trial_Duration=mv_mpP(1).Header.num_frames/mv_mpP(1).Header.Frame_Rate;
mv_mpP(1).MedFiltSize=temporal_filter_size;
mv_mpP(1).Analog_Freq=invoke(mfile,'ReadParameter','Analog Acquisition Frequency (Hz)');
% mv_mpP(1).Analog_samples=round(mv_mpP(1).Trial_Duration*str2double(mv_mpP(1).Analog_Freq(1:(end-2))),0);
mv_mpP(1).Analog_samples=str2double(invoke(mfile,'ReadParameter','Analog Sample Count')); %New MCSX has this function built in use instead
mv_mpP(1).BehaviorCamera.xsize=480;
mv_mpP(1).BehaviorCamera.ysize=640;
mv_mpP(1).BehaviorCamera.Frame_num=str2double(invoke(mfile,'ReadParameter','Video Image Count'));

%Get analog signals and parameters
fprintf('Extracting analog input signals\n')
begin=tic;
for chan=1:Analog_Channel_Num
    call=['Analog Ch ' num2str((chan-1)) ' Name'];
    mv_mpP(1).Analog_Channel_Name{chan}=invoke(mfile,'ReadParameter',call);
    thesignal=strrep([mv_mpP(1).Analog_Channel_Name{chan} '_Analog_Signal'],' ', '_');
    mv_mpP(1).AnalogSignals.(thesignal)=double(invoke(mfile,'ReadAnalog',(chan),mv_mpP(1).Analog_samples,0));
%     if isnan(mv_mpP(1).AnalogSignals.(thesignal))
%        samples=round((mv_mpP.Trial_Duration-0.01)*str2double(mv_mpP.Analog_Freq(1:(end-2))),0);
%        mv_mpP(1).AnalogSignals.(thesignal)=double(invoke(mfile,'ReadAnalog',(chan),samples,0));
%        mv_mpP(1).AnalogSignals.(thesignal)((samples+1):mv_mpP(1).Analog_samples)=0;
%     end
end
AnalogTime=toc(begin);
fprintf(['Analog signal upload time ' num2str(AnalogTime) ' seconds\n']);
% Emtpy Data structures for images
HoldFrames((mv_mpP(1).xsize),(mv_mpP(1).ysize),ChunkSize)=int16(0);
proctime(1:(mv_mpP(1).endframe))=0;
transfertime(1:ceil(mv_mpP(1).endframe/ChunkSize))=0;
%generate an n-frame FIFO buffer to do the temporal median filtering on
buffered_frames=NaN*ones(mv_mpP(1).xsize ,mv_mpP(1).ysize , temporal_filter_size);
the_midpoint=ceil(temporal_filter_size/2);
n_frames=length(mv_mpP(1).startframe:mv_mpP(1).endframe);
mv_mpP(1).pixel_shift=zeros(4,n_frames);
FirstImage((1:(mv_mpP.xsize-PixelOffset)),mv_mpP.ysize)=NaN;
NewImage((1:(mv_mpP.xsize-PixelOffset)),mv_mpP.ysize)=NaN;
TempImage((1:mv_mpP.xsize),(1:mv_mpP.ysize))=NaN;
for channel=1:Imaging_channels
    thetext=['Motion correcting and median filtering channel number ' num2str(channel) ' of ' num2str(Imaging_channels)];
    fprintf([thetext '\n'])
    if Imaging_channels==1
%         if strcmpi(thechannel,'GFP')
            chan_name=['Scanning Ch ' num2str(channel) ' Name'];
%         else
%             chan_name=['Scanning Ch ' num2str(channel-1) ' Name'];
%         end
    else
        chan_name=['Scanning Ch ' num2str(channel-1) ' Name'];
    end
    %this will get which channel (RFP,GFP) you are processing
    mv_mpP(1).channel_name{channel}=invoke(mfile,'ReadParameter',chan_name);
    thechannel=mv_mpP(1).channel_name{channel};
    
    %Preallocates memory for current channel
    TheImagingData.(thechannel)((mv_mpP(1).xsize), (mv_mpP(1).ysize) ,(mv_mpP(1).Header.num_frames))=int16(0);
    
    %2-D fft of first frame used for image registration
    if channel==1
        if strcmpi(thechannel,'GFP')
            firstframe=max(0,double(invoke(mfile,'ReadFrame',(channel+1),1)));%Some pixels in background may be slightly less than zero, this forces them to 0.
            EvenColumns=firstframe(:,(2:2:end));
            OddColumns=firstframe(:,(1:2:(end-1)));
            FirstImage(:,(1:2:(end-1)))=OddColumns(1:(end-PixelOffset),:);
            FirstImage(:,(2:2:end))=EvenColumns(((PixelOffset+1):end),:);
            FirstImage((1:(mv_mpP.rowStart)),:)=NaN;
            FirstImage(((mv_mpP.rowEnd):mv_mpP.ysize),:)=NaN;
            First_Frame_FFT=fft2(FirstImage);
        else
            firstframe=max(0,double(invoke(mfile,'ReadFrame',channel,1)));%Some pixels in background may be slightly less than zero, this forces them to 0.
            EvenColumns=firstframe(:,(2:2:end));
            OddColumns=firstframe(:,(1:2:(end-1)));
            FirstImage(:,(1:2:(end-1)))=OddColumns((1:(end-PixelOffset)),:);
            FirstImage(:,(2:2:end))=EvenColumns(((PixelOffset+1):end),:);
            FirstImage((1:(mv_mpP.rowStart)),:)=NaN;
            FirstImage(((mv_mpP.rowEnd):mv_mpP.ysize),:)=NaN;
            First_Frame_FFT=fft2(FirstImage);
        end
    end
    
    
    count_frames=0;
    chunk_count=0;
    for frame_num=1:mv_mpP.Header.num_frames
        start=tic;
        count_frames=count_frames+1;
        %load frame
        if channel==1
            if strcmpi(thechannel,'GFP')
                Theframe=max(0,double(invoke(mfile,'ReadFrame',(channel+1),frame_num)));%Some pixels in background may be slightly less than zero, this forces them to 0.
            else
                Theframe=max(0,double(invoke(mfile,'ReadFrame',channel,frame_num)));%Some pixels in background may be slightly less than zero, this forces them to 0.
            end
        else
            Theframe=max(0,double(invoke(mfile,'ReadFrame',channel,frame_num)));%Some pixels in background may be slightly less than zero, this forces them to 0.
        end
        % Correct Pixel Shift
        EvenColumns=Theframe(:,(2:2:end));
        OddColumns=Theframe(:,(1:2:(end-1)));
        NewImage(:,(1:2:(end-1)))=OddColumns((1:(end-PixelOffset)),:);
        NewImage(:,(2:2:end))=EvenColumns(((PixelOffset+1):end),:);
        TempImage((mv_mpP.rowStart:(mv_mpP.rowEnd-PixelOffset)),:)=NewImage((mv_mpP.rowStart:(mv_mpP.rowEnd-PixelOffset)),:);
        Theframe=TempImage;
        if channel==1
            %2-D fft of frame
            fft_frame=fft2(Theframe);
            %Motion correction
            [mv_mpP(1).pixel_shift(:,frame_num), Greg]=dftregistration(First_Frame_FFT,fft_frame,1);
        end
        moco_frame(max(1,1+mv_mpP(1).pixel_shift(3,frame_num)):min(mv_mpP(1).xsize,mv_mpP(1).xsize+mv_mpP(1).pixel_shift(3,frame_num)),...
            max(1,1+mv_mpP(1).pixel_shift(4,frame_num)):min(mv_mpP(1).ysize,mv_mpP(1).ysize+mv_mpP(1).pixel_shift(4,frame_num)))=...
            Theframe(max(1,1-mv_mpP(1).pixel_shift(3,frame_num)):min(mv_mpP(1).xsize,mv_mpP(1).xsize-mv_mpP(1).pixel_shift(3,frame_num)),...
            max(1,1-mv_mpP(1).pixel_shift(4,frame_num)):min(mv_mpP(1).ysize,mv_mpP(1).ysize-mv_mpP(1).pixel_shift(4,frame_num)));
        framedim=size(moco_frame);
        % Correct for motion correction by appending NaN to shifted columns
        if framedim(1)<mv_mpP.xsize
            fprintf('Something went wrong with image registration\n')
%             if mv_mpP(1).pixel_shift(3,frame_num)<0
%                 moco_frame(((framedim(1)+1):mv_mpP.xsize),:)=NaN;
%             else
%                 moco_frame(((mv_mpP(1).xsize-framedim(1)):end),:)=mocoframe;
%                 moco_frame((1:(mv_mpP(1).xsize-framedim(1))),:)=NaN;
%             end
        end
         if framedim(2)<mv_mpP.ysize
             fprintf('Something went wrong with image registration\n')
%             if mv_mpP(1).pixel_shift(4,frame_num)<0
%                 moco_frame(:,((framedim(2)+1):mv_mpP.ysize))=NaN;
%             else
%                 moco_frame(:,((mv_mpP(1).ysize-framedim(2)):end))=mocoframe;
%                 moco_frame(:,(1:(mv_mpP(1).ysize-framedim(2))))=NaN;
%             end
        end
        if channel==1
            %first-in, first-out buffer
            buffered_frames(:,:,1:temporal_filter_size-1)=buffered_frames(:,:,2:temporal_filter_size);
            %put in the motion corrected frame
            buffered_frames(:,:,temporal_filter_size)=moco_frame;
            %perform an n-point temporal median filter
            %medfilted_frame=squeeze(medfilt1(buffered_frames,temporal_filter_size,[],3, 'omitnan'));
            %HoldFrames(:,:,count_frames)=int16(squeeze(medfilted_frame(:,:,the_midpoint)));
            medfilt_frame=squeeze(gather(median(gpuArray(buffered_frames),3,'omitnan')));
            HoldFrames(:,:,count_frames)=medfilt_frame;
        else
            HoldFrames(:,:,count_frames)=moco_frame;
        end
        proctime(channel,frame_num)=toc(start);
        %Write frame chunks of 2p data to memory
        if count_frames==ChunkSize
            go=tic;
            lastframe=chunk_count*ChunkSize;
            %Write registered and filtered frame to memory
            fprintf('Writing images to memory\n')
            TheImagingData.(thechannel)(:,:,((lastframe+1):frame_num))=HoldFrames;
            transfertime(channel,(chunk_count+1))=toc(go);
            %reset counters
            count_frames=0;
            chunk_count=chunk_count+1;
        end
    end
    fprintf(['Time to process all frames ' num2str((sum(proctime(channel,:))/60)) ' minutes\n'])
    fprintf(['Average per frame process time ' num2str(mean(proctime(channel,:))) ' seconds\n'])
    fprintf(['Total time to transfer data to memory ' num2str((sum(transfertime(channel,:))/60)) ' minutes\n'])
end
fprintf('Saving mv_mpP\n')
save([filebase '_mv_mpP'],'mv_mpP','-v7.3');
behavior_start=tic;
BehaviorCamera=VideoWriter([filebase '_behavior_camera.avi'],'Uncompressed AVI');
BehaviorCamera.FrameRate=mv_mpP(1).Header.Frame_Rate;
open(BehaviorCamera);
for next_frame=1:mv_mpP(1).BehaviorCamera.Frame_num
    try
    movie_frame=invoke(mfile,'ReadVideoFrame',next_frame);
    bmovie_frame = movie_frame < 0;
    nmovie_frame=zeros(size(movie_frame));
    nmovie_frame=double(movie_frame)+256*bmovie_frame;
    writeVideo(BehaviorCamera,uint8(nmovie_frame));
    catch
        fprintf(['Error reading frame ' num2str(next_frame) '\n'])
    end
end
close(BehaviorCamera);
movie_frame=[];
behavior_transfer=toc(behavior_start);
fprintf(['Total time to transfer video to memory ' num2str((behavior_transfer/60)) ' minutes\n'])
close all;
end
