function [mv_mpP]=GetDiametersMDF_preMocoed_BehaviorSeg_002(fname,mv_mpPName,Vessel_Channel,animal,Objective,vessel_ID,varargin)
%this function takes a already motion corrected (and temporally medial
%filtered) image structure (movie) and accompanying .mat file (mv_mpP), and
%extracts the diameters from surface and penetrating vessels with GCaMP
%data

%mv_mpP is the processed output with GCaMP and vessel diameters.

%Patrick Drew 6/2018
%Edited by Kyle Gheres to process MDF files and GCaMP recordings 07/2018

%open file and display first frame
if exist('fname','var')==0
    fprintf('Select median filtered image to analyze\n')
    fname=uigetfile('*MedFilt.mat','MultiSelect','off');
end
if exist('mv_mpPName','var')==0
    %fprintf('Select matching analog data file\n')
    %mv_mpPName=uigetfile('*mv_mpP.mat','MultiSelect','off');
    mv_mpPName=[fname(1:11) 'mv_mpP.mat'];
end
if exist('Vessel_Channel','var')==0
    Vessel_Channel=input('What imaging channel are the blood vessels on?[GFP/RFP]','s');
end
load(mv_mpPName);
Channel_num=numel(mv_mpP.channel_name);
ImageData=matfile(fname,'Writable',true);
theimage=median(double(ImageData.(mv_mpP.channel_name{1})(:,:,(1:120))),3);
plotimage=max(0,(theimage-0.5*median(theimage)));

figure(2)
imagesc(plotimage)
colormap hot
caxis([0 500]);
axis image
axis off
% try
%     walking_analog=mv_mpP.AnalogSignals.Rotary_Encoder_Analog_Signal;
%     Contra_Sol_analog=mv_mpP.AnalogSignals.Contralateral_Solenoid_Analog_Signal;
%     Control_Sol_analog=mv_mpP.AnalogSignals.Control_Solenoid_Analog_Signal;
% catch
%     walking_analog=[];
%     Contra_Sol_analog=[];
%     Control_Sol_analog=[];
% end
if exist('animal','var')==0
    animal=input('animal?','s');
end
nvessels=input('# of vessels on this slice:');
gcamp=input('Did you image GCaMP during trial [y/n]\n','s');
if exist('Objective','var')==0
%     microns_per_pixel=input('microns per pixel @ 1x');
%     mv_mpP(1).Microns_Per_Pixel=microns_per_pixel;
    Objective=input('objective: [1]Nikon 16x, 0.8NA  [2]Olympus 20x,0.95NA  [3]Nikon 4x');
end
switch Objective
    case 1
        mv_mpP(1).Microns_Per_Pixel=0.825;
        mv_mpP(1).Objective='Nikon 16x 0.8NA';
    case 2
        mv_mpP(1).Microns_Per_Pixel=0.64;
        mv_mpP(1).Objective='Olympus 20x 0.95NA';
    case 3
        mv_mpP(1).Microns_Per_Pixel=2.763;
        mv_mpP(1).Objective='Nikon 4x 0.8NA';
end


% the_scanmirrors=2;%input('which mirrors? 1)6210(fast) 2)6215(slow)')
frame_rate=mv_mpP(1).Header.Frame_Rate;%str2num(mv_mpP(1).Header.Frame_Rate);
% Pixel_clock=1/((5/4)*frame_rate*mv_mpP(1).xsize*mv_mpP(1).ysize)
% Pixel_clock=mv_mpP(1).Header.Pixel_Clock;
% time_per_line=(mv_mpP(1).Header.Frame_Width)*(5/4)*Pixel_clock*(.05*1e-6);
mv_mpP(1).Header.time_per_line=1/(frame_rate*(mv_mpP(1).Header.Frame_Height));
Xfactor= mv_mpP(1).Microns_Per_Pixel/((str2double(mv_mpP(1).Header.Magnification(1:(end-1)))));
if (nvessels>0)
    
    ystring=['y'];
    theinput=['n'];
    xsize=size(theimage,2);
    ysize=size(theimage,1);
    Y=repmat([1:ysize]',1,xsize);
    X=repmat([1:xsize],ysize,1);
    
    for vesselnumber=1:nvessels
        mv_mpP(vesselnumber).Xfactor=Xfactor;
        if exist('vessel_ID','var')==0
            vessel_ID=input('vessel ID?' ,'s');
        end
        mv_mpP(vesselnumber).vessel_type=input('What type of Vessel are you measuring? [surface/penetrating]','s');
        mv_mpP(vesselnumber)=mv_mpP(1);
        mv_mpP(vesselnumber).vessel_ID=vessel_ID;
        the_date=date;
        switch mv_mpP(vesselnumber).vessel_type
            case 'surface'
                figure(2)
                area = impoly(gca,[1 1; 1 20;20 20;20 1]);
                while (strcmp(ystring,theinput)~=1)
                    theinput=input('diameter box ok? y/n \n','s');
                end
                
                if    strcmp(ystring,theinput)
                    api = iptgetapi(area);
                    mv_mpP(vesselnumber).Vessel.box_position.xy=api.getPosition();
                    mv_mpP(vesselnumber).Vessel.xsize=xsize;
                    mv_mpP(vesselnumber).Vessel.ysize=ysize;
                    theinput=['n'];
                end
                diam_axis=imline(gca,round(xsize*[.25 .75]),round(ysize*[.25 .75]));
                while (strcmp(ystring,theinput)~=1)
                    theinput=input('diameter axis ok? y/n \n','s');
                end
                if    strcmp(ystring,theinput)
                    api = iptgetapi(diam_axis);
                    mv_mpP(vesselnumber).Vessel.vessel_line.position.xy=api.getPosition();
                    theinput=['n'];
                end
                %mv_mpP(vesselnumber).Xfactor=Xfactor;
                clear vessel_ID
                mv_mpP=GetDiametersFromMovie_MDF_01(mv_mpP, ImageData,Vessel_Channel);
                try
                    [mv_mpP]=FWHMfromMovieProjection2(mv_mpP,[mv_mpP(1).startframe mv_mpP(1).endframe] );
                catch
                    disp('plot fail')
                end
            case 'penetrating'
                mv_mpP(1).first_frame=ImageData.(Vessel_Channel)(:,:,1);
                maxframe=mv_mpP.endframe;
                    area = impoly(gca,[1 1; 1 20;20 20;20 1]);
                while (strcmp(ystring,theinput)~=1)
                    theinput=input('diameter box ok? y/n \n','s');
                end
                
                if    strcmp(ystring,theinput)
                    api = iptgetapi(area);
                    mv_mpP(vesselnumber).Vessel.box_position.xy=api.getPosition();
                    mv_mpP(vesselnumber).Vessel.xsize=xsize;
                    mv_mpP(vesselnumber).Vessel.ysize=ysize;
                    theinput=['n'];
                end
                plot_every=1;
                [mv_mpP]=PenetratingVesselAreaTiff(mv_mpP,ImageData,Vessel_Channel, plot_every,maxframe);
        end
        %     try
        %         mv_mpP(1).Analog.walking=walking_analog;
        %     catch
        %         mv_mpP(1).Analog.walking=[];
        %     end
        
    end
    % ntrials=1;%input
    
end

%Draw ROI for GCaMP analysis
if strcmpi(gcamp,'y')==1
    mv_mpP.GCaMP.GCaMPImage=std(double(ImageData.GFP(:,:,(270:330))),0,3);
    figure(3);imagesc(mv_mpP.GCaMP.GCaMPImage);
    ROI_num=input('How many ROI would you like to draw?');
    for ROI_count=1:ROI_num
        ROI_ok='n';
        while strcmpi(ROI_ok,'n')==1
            GCaMP_ROI=impoly(gca,[1 1; 1 20; 20 20 ; 20 1]);
            ROI_ok=input('ROI OK [y/n]','s');
        end
        mv_mpP.GCaMP.ROI_xy(:,:,ROI_count)=getPosition(GCaMP_ROI);
    end
 %Analyze pixel intensity of GCaMP ROI   
    for p=1:ROI_num
        [mv_mpP]=GCaMPIntensityCalc(mv_mpP,ImageData,p,ROI_num);
    end
end

[mv_mpP]=BehaviorChunkData(mv_mpP,gcamp);

[mv_mpP]=PlotMovieDiameters_walking(mv_mpP,1,gcamp,animal,fname,the_date,plotimage);

try
    save_filename=[animal '_mv_mpP_' fname(1:end-4) '_' the_date];
    save(save_filename,'mv_mpP','-v7.3')%save processed file to current directory
    %    save(['f:\data\crunched\' save_filename],'mv_mpP')%save to crunched directory
catch
    disp('crunched folder save failed!')
end
end


function [mv_mpP]=PlotMovieDiameters_walking(mv_mpP,basefignum,gcamp,animal,fname,the_date,plotimage)
%plots episodic movie diameters
%Binarize behavior
StimTime=(1:length(mv_mpP(1).AnalogSignals.Contralateral_Solenoid_Analog_Signal))/str2double(mv_mpP(1).Analog_Freq(1:(end-2)));
Ball_Vel=mv_mpP(1).AnalogSignals.Rotary_Encoder_Analog_Signal*10^-3;%output is in uV convert to V
[imp_bin]=velocity_binarize(Ball_Vel,str2double(mv_mpP(1).Analog_Freq(1:(end-2))),str2double(mv_mpP(1).Analog_Freq(1:(end-2))),1e-4);
Loco_Stims(1:length(mv_mpP(1).AnalogSignals.Rotary_Encoder_Analog_Signal))=NaN;
Loco_Stims(imp_bin==1)=1;

Contra_Stims(1:length(mv_mpP(1).AnalogSignals.Contralateral_Solenoid_Analog_Signal))=NaN;
Contra_Times=diff(mv_mpP(1).AnalogSignals.Contralateral_Solenoid_Analog_Signal)>100;
Contra_Stims(Contra_Times==1)=1;

Ctrl_Stims(1:length(mv_mpP(1).AnalogSignals.Control_Solenoid_Analog_Signal))=NaN;
Ctrl_Times=diff(mv_mpP(1).AnalogSignals.Control_Solenoid_Analog_Signal)>100;
Ctrl_Stims(Ctrl_Times==1)=1;

%walking data is plotted as well
nvessels=length(mv_mpP);
cm_per_volt=38.4;
% stimcolor={'b' 'r'};
vessel_colors={'b' 'c' 'm' 'y' 'r' 'k' 'b' 'c' 'm' 'y' 'r' 'k' };
GCaMP_colors={'g' 'm' 'c' 'r' 'b' 'g' 'm' 'c' 'r' 'b'};
figure(basefignum)
subplot(3,2,6)
hold off
imagesc(((1:size(mv_mpP(1).first_frame,2))*mv_mpP(1).Microns_Per_Pixel/(str2double(mv_mpP(1).Header.Magnification(1:(end-1))))),((1:size(mv_mpP(1).first_frame,1))*mv_mpP(1).Microns_Per_Pixel/(str2double(mv_mpP(1).Header.Magnification(1:(end-1))))), ...
    plotimage);%medfilt2(double(max(0,(mv_mpP(1).first_frame-.5*median(mv_mpP(1).first_frame(:)))))))
axis image
colormap hot
caxis([0 500]);
hold on
axis off
title('Vessel diameter ROI')
duration=length(mv_mpP(1).Vessel.diameter(2:end))/(mv_mpP(1).Header.Frame_Rate);
the_fs=1/(mv_mpP(1).Header.Frame_Rate);
% the_NW=floor(max(1,the_fs*duration/2));

for mv=1:nvessels
    if strcmpi(mv_mpP(mv).vessel_type,'surface')==1
    mv_mpP(mv).Vessel.diameter=(mv_mpP(1).Microns_Per_Pixel/(str2double(mv_mpP(1).Header.Magnification(1:(end-1)))))*mv_mpP(mv).Vessel.raw_diameter;
    plot(mv_mpP(1).Microns_Per_Pixel/str2double(mv_mpP(1).Header.Magnification(1:(end-1)))*[mv_mpP(mv).Vessel.box_position.xy(:,1)' mv_mpP(mv).Vessel.box_position.xy(1,1)],...
        mv_mpP(1).Microns_Per_Pixel/str2double(mv_mpP(1).Header.Magnification(1:(end-1)))*[mv_mpP(mv).Vessel.box_position.xy(:,2)' mv_mpP(mv).Vessel.box_position.xy(1,2)],vessel_colors{mv},'LineWidth',2);
    else
    plot(mv_mpP(1).Microns_Per_Pixel/str2double(mv_mpP(1).Header.Magnification(1:(end-1)))*[mv_mpP(mv).Vessel.box_position.xy(:,1)' mv_mpP(mv).Vessel.box_position.xy(1,1)],...
        mv_mpP(1).Microns_Per_Pixel/str2double(mv_mpP(1).Header.Magnification(1:(end-1)))*[mv_mpP(mv).Vessel.box_position.xy(:,2)' mv_mpP(mv).Vessel.box_position.xy(1,2)],vessel_colors{mv},'LineWidth',2); 
    end
end
max_plot_time=max(length(mv_mpP(mv).Vessel.diameter(1:end))/(mv_mpP(1).Header.Frame_Rate));
% frame_times=((1:length(mv_mpP(mv).Vessel.diameter(2:end)))/(mv_mpP(1).Header.Frame_Rate));       
lowFilt=1; %Frequency to low pass GCaMP signal for plotting
[b,a]=butter(3,(lowFilt/(0.5*mv_mpP(1).Header.Frame_Rate)),'low');
if strcmpi(gcamp,'y')==1
    ax1=subplot(3,1,1); % Keep for GCaMP signal
    hold on   
    try
        GCaMP_Time=(1:length(mv_mpP.GCaMP.Normalized_ROI))/mv_mpP.Header.Frame_Rate;
 
        for ROI_plot=1:size(mv_mpP.GCaMP.Normalized_ROI,2)
            plot(GCaMP_Time,filtfilt(b,a,mv_mpP.GCaMP.Normalized_ROI(:,ROI_plot)),GCaMP_colors{ROI_plot},'LineWidth',2);
        end
        title('GCaMP intensity');
        xlabel('sec')
        ylabel('deltaF/Fo (%)');
        Loco_GStims=Loco_Stims+(1.1*max(mv_mpP(1).GCaMP.Normalized_ROI(:,ROI_plot)));
        Contra_GStims=Contra_Stims+(1.05*max(mv_mpP(1).GCaMP.Normalized_ROI(:,ROI_plot)));
        Ctrl_GStims=Ctrl_Stims+(1.05*max(mv_mpP(1).GCaMP.Normalized_ROI(:,ROI_plot)));
        scatter(StimTime,Loco_GStims,'filled','k');scatter(StimTime,Contra_GStims,'filled','r');scatter(StimTime,Ctrl_GStims,'filled','c');
        xlim([0 max_plot_time])
        ylim([(0.9*min(mv_mpP.GCaMP.Normalized_ROI(:,ROI_plot))) (1.1*max(mv_mpP(1).GCaMP.Normalized_ROI(:,ROI_plot)))+1])
    catch
        fprintf('Something went wrong plotting GCaMP data\n')
    end
    xlabel('sec')
    ylabel('% change')
    %axis([0 max_plot_time (min(mv_mpP(1).Vessel.diameter)-0.1*min(mv_mpP(1).Vessel.diameter)) 100 ])
end
ax2=subplot(3,1,2);
hold off
params.Fs=(mv_mpP(1).Header.Frame_Rate);
params.tapers=[1 1];%[the_NW 2*the_NW-1];
params.err=[1 .01];
params.fpass=[0.03 3];
for mv=1:nvessels
    plot((1:length(mv_mpP(mv).Vessel.diameter(2:end)))/(mv_mpP(1).Header.Frame_Rate),filtfilt(b,a,mv_mpP(mv).Vessel.diameter(2:end)),[vessel_colors{mv}],'LineWidth',2);
    hold on
    [S{mv},f{mv},Serr{mv}]=mtspectrumc(mv_mpP(mv).Vessel.diameter(2:end)-mean(mv_mpP(mv).Vessel.diameter(2:end)),params);
    [diff_S{mv},diff_f{mv},diff_Serr{mv}]=mtspectrumc(diff(mv_mpP(mv).Vessel.diameter(2:end)),params);
    [volume_S{mv},volume_f{mv},volume_Serr{mv}]=mtspectrumc((mv_mpP(mv).Vessel.diameter(2:end).^2),params);
    mv_mpP(mv).S=S{mv};
    mv_mpP(mv).f=f{mv};
    mv_mpP(mv).diff_S=diff_S{mv};
    mv_mpP(mv).diff_f=diff_f{mv};
    mv_mpP(mv).volume_S=volume_S{mv};
    mv_mpP(mv).volume_f=volume_f{mv};
end
Loco_Stims=Loco_Stims+(1.1*max(mv_mpP(1).Vessel.diameter));
Contra_Stims=Contra_Stims+(1.05*max(mv_mpP(1).Vessel.diameter));
Ctrl_Stims=Ctrl_Stims+(1.05*max(mv_mpP(1).Vessel.diameter));
scatter(StimTime,Loco_Stims,'filled','k');scatter(StimTime,Contra_Stims,'filled','r');scatter(StimTime,Ctrl_Stims,'filled','c');
axis([0 max_plot_time 0.9*min(mv_mpP(mv).Vessel.diameter) (1.1*max(mv_mpP(mv).Vessel.diameter)+1)]);
title(['vessel diameters'])
xlabel('time, seconds')
ylabel('diameter, um')
title(mv_mpP(1).Header.Filename)
if strcmpi(gcamp,'y')==1
linkaxes([ax1,ax2],'x');
subplot(3,2,5)
% hold off
% for mv=1:nvessels
%     loglog(f{mv},S{mv},vessel_colors{mv});
%     hold on
% end
% axis([ 0.03 3  1e-4 5])
% axis square

% subplot(3,2,4)
imagesc(mv_mpP(1).GCaMP.GCaMPImage);
colormap gray
axis image
hold on
axis off

for num_ROI=1:size(mv_mpP(1).GCaMP.ROI_xy,3)
plot([mv_mpP(1).GCaMP.ROI_xy(:,1,num_ROI)' mv_mpP(1).GCaMP.ROI_xy(1,1,num_ROI)],[mv_mpP(1).GCaMP.ROI_xy(:,2,num_ROI)' mv_mpP(1).GCaMP.ROI_xy(1,2,num_ROI)],GCaMP_colors{num_ROI},'LineWidth',2);
end
title('GCaMP ROIs')
end
save_figname=[animal '_mv_mpP_' fname(1:end-4) '_' the_date '_Diameters'];
savefig(save_figname);
end


function mv_mpP=GetDiametersFromMovie_MDF_01(mv_mpP,ImageData,Vessel_Channel)
%this function opens the tiff file and gets the  vessel projections from the defined polygons
%mv_mpP(1).first_frame=imread(fname,'TIFF','Index',1)
substructs=who(ImageData);
dummy=strcmpi(substructs,Vessel_Channel);
The_Channel=substructs(dummy==1);
mv_mpP(1).first_frame=ImageData.(The_Channel{1})(:,:,1);

%mp2mat_getChannelData_narrow_noplot(mpfile,1,mv_mpP(1), [1 1],[1 mv_mpP(1).xsize]);
% fft_first_frame=fft2(double(mv_mpP(1).first_frame));
for mv=1:length(mv_mpP)
    mv_mpP(mv).Vessel.shifted_bounded_raw_frame(mv_mpP(1).xsize,mv_mpP(1).ysize,((mv_mpP(1).startframe):(mv_mpP(1).endframe)))=uint16(0);
    stop((mv_mpP(1).startframe):(mv_mpP(1).endframe))=0;
    Y=repmat([1:mv_mpP(1).ysize]',1,mv_mpP(1).xsize);
    X=repmat([1:mv_mpP(1).xsize],mv_mpP(1).ysize,1);
    mv_mpP(mv).Vessel.projection_angle=atand(diff(mv_mpP(mv).Vessel.vessel_line.position.xy(:,1))/diff(mv_mpP(mv).Vessel.vessel_line.position.xy(:,2)));
    %atand(diff(mv_mpP(mv).Vessel.vessel_line.position.xy(:,1))/diff(mv_mpP(mv).Vessel.vessel_line.position.xy(:,2)))
    
    mv_mpP(mv).Vessel.reshifted_ROI=zeros(mv_mpP(1).xsize,mv_mpP(1).ysize,mv_mpP(1).xsize,length(mv_mpP));
    % mv_mpP(mv).Vessel.reshifted_ROI_medfilt=zeros(mv_mpP(1).xsize,mv_mpP(1).ysize,mv_mpP(1).xsize,length(mv_mpP));
    % register all the frames
    for theframe=(mv_mpP(1).startframe):(mv_mpP(1).endframe)
        start=tic;
        %raw_frame =imread(fname,'Index',theframe);%(mp2mat_getChannelData_narrow_noplot(mpfile,1,mv_mpP(1), [theframe theframe],[1 mv_mpP(1).xsize]));
        raw_frame=uint16(ImageData.(The_Channel{1})(:,:,theframe));
        % NO motion correction
        if theframe==1
            inpoly_frame = inpolygon(X,Y,mv_mpP(mv).Vessel.box_position.xy(:,1),mv_mpP(mv).Vessel.box_position.xy(:,2));
        end
        mv_mpP(mv).Vessel.shifted_bounded_raw_frame(:,:,theframe)=raw_frame.*uint16(inpoly_frame);%need to crop
        mv_mpP(mv).Vessel.projection(theframe,:)=radon(squeeze(mv_mpP(mv).Vessel.shifted_bounded_raw_frame(:,:,theframe)),mv_mpP(mv).Vessel.projection_angle);
        stop(theframe)=toc(start);
    end
    total_radon_time=num2str(sum(stop)/60);
    avg_radon_time=num2str(mean(stop));
    fprintf(['Radon projection time ' total_radon_time ' minutes\n']);
    fprintf(['Average Radon projection time per frame ' avg_radon_time ' seconds\n']);
end

figure(4444)
imagesc( mv_mpP(mv).Vessel.projection)

end


function [mv_mpP]=FWHMfromMovieProjection2(mv_mpP,theframes)
%gets dameters from movie projections
vessel_colors={'b' 'r' 'g' 'k' 'c' 'm' 'y'};
stimcolor={'b' 'r'};
for mv=1:length(mv_mpP)
    for k=min(theframes):max(theframes)
        mv_mpP(mv).Vessel.raw_diameter(k)=calcFWHM(mv_mpP(mv).Vessel.projection(k,:));
    end
    
    mv_mpP(mv).Vessel.diameter= mv_mpP(mv).Vessel.raw_diameter*mv_mpP(1).Xfactor;
end
end


function width = calcFWHM(data,smoothing,threshold)
% function which takes data and calculates the full-width, half max value
% half-max values are found looking in from the sides, i.e., the program will work
% even if the data dips to a lower value in the middle
% 2009-06-17 - no longer subracting the min from the entire data set
% 2009-06-17 - changed convolution to conv2, and added 'valid' parameter

data = double(data(:));     % make sure this is column, and cast to double

% smooth data, if appropriate
if nargin < 2
    % smoothing not passed in, set to default (none)
    smoothing = 1;
end

if smoothing > 1
    data = conv2(data,rectwin(smoothing) ./ smoothing,'valid');
end

% subtract out baseline
%data = data - min(data);   %jd - don't subract out min, in case threshold was set externally

if nargin < 3
    %threshold = max(data)/2;
    offset = min(data);                           % find the baseline
    threshold = max(data - offset) / 2 + offset;  % threshold is half max, taking offset into account
end

aboveI = find(data > threshold);    % all the indices where the data is above half max

if isempty(aboveI)
    % nothing was above threshold!
    width = 0;
    return
end

firstI = aboveI(1);                 % index of the first point above threshold
lastI = aboveI(end);                % index of the last point above threshold

if (firstI-1 < 1) || (lastI+1) > length(data)
    % interpolation would result in error, set width to zero and just return ...
    width = 0;
    return
end

% use linear intepolation to get a more accurate picture of where the max was
% find value difference between the point and the threshold value,
% and scale this by the difference between integer points ...
point1offset = (threshold-data(firstI-1)) / (data(firstI)-data(firstI-1));
point2offset = (threshold-data(lastI)) / (data(lastI+1)-data(lastI));

point1 = firstI-1 + point1offset;
point2 = lastI + point2offset;

width = point2-point1;
end

function [mv_mpP]=GCaMPIntensityCalc(mv_mpP,ImageData,ROI_count,ROI_num)
XVals=repmat([1:mv_mpP(1).xsize],mv_mpP(1).ysize,1);
YVals=repmat([1:mv_mpP(1).ysize]',1,mv_mpP(1).xsize);
%find pixels inside ROI
insidepoly=inpolygon(XVals,YVals,mv_mpP(1).GCaMP.ROI_xy(:,1,ROI_count),mv_mpP(1).GCaMP.ROI_xy(:,2,ROI_count));
PixNum=numel(find(insidepoly==1));
PolyMask=double(insidepoly);
PolyMask(PolyMask==0)=NaN;
if ROI_num==1
mv_mpP(1).GCaMP.ROI_Intensity((1:mv_mpP(1).endframe),(1:ROI_count))=0;
mv_mpP(1).GCaMP.Normalized_ROI((1:mv_mpP(1).endframe),(1:ROI_count))=0;
end
%Get average pixel intensity inside ROI for each frame
start=tic;
for frame_num=1:mv_mpP.endframe
ROI_Pix=double(uint16(ImageData.GFP(:,:,frame_num))).*PolyMask;
mv_mpP(1).GCaMP.ROI_Intensity(frame_num,ROI_count)=sum(sum(ROI_Pix,2,'omitnan'),1,'omitnan')/PixNum;
end
ROI_Intensity_Time=toc(start);
fprintf(['GCaMP signal intensity processing time ' num2str(ROI_Intensity_Time/60) ' minutes\n'])
fprintf(['Average time to per frame process GCaMP signal ' num2str(ROI_Intensity_Time/mv_mpP.endframe) ' seconds\n'])
mv_mpP(1).GCaMP.Normalized_ROI(:,ROI_count)=((mv_mpP(1).GCaMP.ROI_Intensity(:,ROI_count)-mean(mv_mpP(1).GCaMP.ROI_Intensity(:,ROI_count),1))/(mean(mv_mpP(1).GCaMP.ROI_Intensity(:,ROI_count),1)))*100;
end

function [mv_mpP]=BehaviorChunkData(mv_mpP,gcamp)
% Once vessel diameters have been calculated and GCaMP intensities
% measured, chunk them around locomotion and whisker stimulus events


RunTime=[];
WinPad=2;%time before behavior/stimulus to include in stim window
WinFollow=5;%time after behavior/stimulus to include in stim window
T_seg=2; %length of time in seconds to consider running event
T_fuse=1; %length of time in seconds between locomotion events to consider one continuous event
T_beg=5; %length of time in seconds between running events to analyze as distinct locomotion events
PuffLength=0.5; %duration of puff in seconds.
PuffDur=PuffLength*str2double(mv_mpP(1).Analog_Freq(1:(end-2)));%Duation of puff in samples;
NormTime=floor(WinPad*mv_mpP(1).Header.Frame_Rate);

%Binarize locomotion, get locomotion times
[imp_bin]=velocity_binarize(mv_mpP(1).AnalogSignals.Rotary_Encoder_Analog_Signal,str2double(mv_mpP(1).Analog_Freq(1:(end-2))));
[T_run,new_T_run]=motion_cont_3(imp_bin,str2double(mv_mpP(1).Analog_Freq(1:(end-2))),T_seg,T_fuse,T_beg);
mv_mpP(1).BehaviorTimes.imp_bin=imp_bin;
mv_mpP(1).BehaviorTimes.Run_Times=new_T_run;
for v=1:size(T_run,2)
    RunTime=[RunTime (T_run(1,v):T_run(2,v))];
end

%Find contralateral whisker stimulus times
TempSol=diff(mv_mpP(1).AnalogSignals.Contralateral_Solenoid_Analog_Signal,1,2);
Ups=find(TempSol>10);
for q=1:numel(Ups)
    TempSol((Ups(q)+1):(Ups(q)+PuffDur))=0;
end
mv_mpP(1).BehaviorTimes.Contra_Puff_Times=find(TempSol>10);

%find control whisker stimulus times
TempSol=diff(mv_mpP(1).AnalogSignals.Control_Solenoid_Analog_Signal,1,2);
Ups=find(TempSol>10);
for q=1:numel(Ups)
    TempSol((Ups(q)+1):(Ups(q)+PuffDur))=0;
end
mv_mpP(1).BehaviorTimes.Ctrl_Puff_Times=find(TempSol>10);

%Chunk vessel diameter measurements around behaviors
for vessel_num=1:size(mv_mpP,2)
    Run_Contra_Count=0;
    Still_Contra_Count=0;
    Run_Ctrl_Count=0;
    Still_Ctrl_Count=0;
    for w=1:size(new_T_run,2)
        Starttime=floor((((new_T_run(1,w)-(WinPad*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);    
        EndTime=floor((((new_T_run(1,w)+(WinFollow*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        if (Starttime-EndTime)~=217
            EndTime=Starttime+217;
        end
        if Starttime>=1
            if EndTime<=size(mv_mpP(vessel_num).Vessel.diameter,2)
                mv_mpP(vessel_num).BehaviorChunkData.Vessel.Run_Trig_Dia(w,:)=mv_mpP(vessel_num).Vessel.diameter(Starttime:EndTime);
                Norm_Const=mean(mv_mpP(vessel_num).BehaviorChunkData.Vessel.Run_Trig_Dia(w,(1:NormTime)),2);
                mv_mpP(vessel_num).BehaviorChunkData.Vessel.Norm_Run_Trig_Dia(w,:)=((mv_mpP(vessel_num).BehaviorChunkData.Vessel.Run_Trig_Dia(w,:)-Norm_Const)/Norm_Const)*100;
            end
        end
    end
    for z=1:size(mv_mpP(1).BehaviorTimes.Contra_Puff_Times,2)
        Starttime=floor((((mv_mpP(1).BehaviorTimes.Contra_Puff_Times(z)-(WinPad*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        EndTime=floor((((mv_mpP(1).BehaviorTimes.Contra_Puff_Times(z)+(WinFollow*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
       if (Starttime-EndTime)~=217
            EndTime=Starttime+217;
        end
        if Starttime>=1
            if EndTime<=size(mv_mpP(vessel_num).Vessel.diameter,2)
                mv_mpP(vessel_num).BehaviorChunkData.Vessel.Contra_Puff_Dia(z,:)=mv_mpP(vessel_num).Vessel.diameter(Starttime:EndTime);
                puffwin=(Starttime:EndTime);
                LocoSeg=ismember(puffwin,RunTime);
                if max(LocoSeg)==1
                    Run_Contra_Count=Run_Contra_Count+1;
                    mv_mpP(vessel_num).BehaviorChunkData.Vessel.Run_Contra_Puff_Dia(Run_Contra_Count,:)=mv_mpP(vessel_num).Vessel.diameter(Starttime:EndTime);
                    Norm_Const=mean(mv_mpP(vessel_num).BehaviorChunkData.Vessel.Run_Contra_Puff_Dia(Run_Contra_Count,(1:NormTime)),2);
                    mv_mpP(vessel_num).BehaviorChunkData.Vessel.Norm_Run_Contra_Puff_Dia(Run_Contra_Count,:)=((mv_mpP(vessel_num).BehaviorChunkData.Vessel.Run_Contra_Puff_Dia(Run_Contra_Count,:)-Norm_Const)/Norm_Const)*100;
                else
                    Still_Contra_Count=Still_Contra_Count+1;
                    mv_mpP(vessel_num).BehaviorChunkData.Vessel.Still_Contra_Puff_Dia(Still_Contra_Count,:)=mv_mpP(vessel_num).Vessel.diameter(Starttime:EndTime);
                    Norm_Const=mean(mv_mpP(vessel_num).BehaviorChunkData.Vessel.Still_Contra_Puff_Dia(Still_Contra_Count,(1:NormTime)),2);
                    mv_mpP(vessel_num).BehaviorChunkData.Vessel.Norm_Still_Contra_Puff_Dia(Still_Contra_Count,:)=((mv_mpP(vessel_num).BehaviorChunkData.Vessel.Still_Contra_Puff_Dia(Still_Contra_Count,:)-Norm_Const)/Norm_Const)*100;
                end
            end
        end
    end
    for p=1:size(mv_mpP(1).BehaviorTimes.Ctrl_Puff_Times,2)
        Starttime=floor((((mv_mpP(1).BehaviorTimes.Ctrl_Puff_Times(p)-(WinPad*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        EndTime=floor((((mv_mpP(1).BehaviorTimes.Ctrl_Puff_Times(p)+(WinFollow*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        if (Starttime-EndTime)~=217
            EndTime=Starttime+217;
        end
        if Starttime>=1
            if EndTime<=size(mv_mpP(vessel_num).Vessel.diameter,2)
                mv_mpP(vessel_num).BehaviorChunkData.Vessel.Ctrl_Puff_Dia(p,:)=mv_mpP(vessel_num).Vessel.diameter(Starttime:EndTime);
                puffwin=(Starttime:EndTime);
                LocoSeg=ismember(puffwin,RunTime);
                if max(LocoSeg)==1
                    Run_Ctrl_Count=Run_Ctrl_Count+1;
                    mv_mpP(vessel_num).BehaviorChunkData.Vessel.Run_Ctrl_Puff_Dia( Run_Ctrl_Count,:)=mv_mpP(vessel_num).Vessel.diameter(Starttime:EndTime);
                    Norm_Const=mean(mv_mpP(vessel_num).BehaviorChunkData.Vessel.Run_Ctrl_Puff_Dia( Run_Ctrl_Count,(1:NormTime)),2);
                    mv_mpP(vessel_num).BehaviorChunkData.Vessel.Norm_Run_Ctrl_Puff_Dia( Run_Ctrl_Count,:)=((mv_mpP(vessel_num).BehaviorChunkData.Vessel.Run_Ctrl_Puff_Dia( Run_Ctrl_Count,:)-Norm_Const)/Norm_Const)*100;
                else
                    Still_Ctrl_Count=Still_Ctrl_Count+1;
                    mv_mpP(vessel_num).BehaviorChunkData.Vessel.Still_Ctrl_Puff_Dia(Still_Ctrl_Count,:)=mv_mpP(vessel_num).Vessel.diameter(Starttime:EndTime);
                    Norm_Const=mean(mv_mpP(vessel_num).BehaviorChunkData.Vessel.Still_Ctrl_Puff_Dia(Still_Ctrl_Count,(1:NormTime)),2);
                    mv_mpP(vessel_num).BehaviorChunkData.Vessel.Norm_Still_Ctrl_Puff_Dia(Still_Ctrl_Count,:)=((mv_mpP(vessel_num).BehaviorChunkData.Vessel.Still_Ctrl_Puff_Dia(Still_Ctrl_Count,:)-Norm_Const)/Norm_Const)*100;
                end
            end
        end
    end
end
if strcmpi(gcamp,'y')==1
%Chunk GCaMP measurements around behaviors
for GCaMP_num=1:size(mv_mpP(1).GCaMP.ROI_xy,3)
    Run_Contra_GCaMP=0;
    Still_Contra_GCaMP=0;
    Run_Ctrl_GCaMP=0;
    Still_Ctrl_GCaMP=0;
    for w=1:size(new_T_run,2)
        Starttime=floor((((new_T_run(1,w)-(WinPad*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        EndTime=floor((((new_T_run(1,w)+(WinFollow*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        if (Starttime-EndTime)~=217
            EndTime=Starttime+217;
        end
        if Starttime>=1
            if EndTime<=size(mv_mpP(vessel_num).GCaMP.ROI_Intensity,1)
                mv_mpP(1).BehaviorChunkData.GCaMP.Run_Trig_ROI_Intensity(w,:,GCaMP_num)=mv_mpP(vessel_num).GCaMP.ROI_Intensity((Starttime:EndTime),GCaMP_num);
                Norm_Const=mean(mv_mpP(1).BehaviorChunkData.GCaMP.Run_Trig_ROI_Intensity(w,(1:NormTime)),2);
                mv_mpP(1).BehaviorChunkData.GCaMP.Norm_Run_Trig_ROI_Intensity(w,:,GCaMP_num)=((mv_mpP(1).BehaviorChunkData.GCaMP.Run_Trig_ROI_Intensity(w,:,GCaMP_num)-Norm_Const)/Norm_Const)*100;
            end
        end
    end
    for z=1:size(mv_mpP(1).BehaviorTimes.Contra_Puff_Times,2)
        Starttime=floor((((mv_mpP(1).BehaviorTimes.Contra_Puff_Times(z)-(WinPad*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        EndTime=floor((((mv_mpP(1).BehaviorTimes.Contra_Puff_Times(z)+(WinFollow*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        if (Starttime-EndTime)~=217
            EndTime=Starttime+217;
        end
        if Starttime>=1
            if EndTime<=size(mv_mpP(vessel_num).GCaMP.ROI_Intensity,1)
                mv_mpP(1).BehaviorChunkData.GCaMP.Contra_Puff_ROI_Intensity(z,:,GCaMP_num)=mv_mpP(1).GCaMP.ROI_Intensity((Starttime:EndTime),GCaMP_num);
                puffwin=(Starttime:EndTime);
                LocoSeg=ismember(puffwin,RunTime);
                if max(LocoSeg)==1
                    Run_Contra_GCaMP=Run_Contra_GCaMP+1;
                    mv_mpP(1).BehaviorChunkData.GCaMP.Run_Contra_Puff_ROI_Intensity(Run_Contra_GCaMP,:,GCaMP_num)=mv_mpP(1).GCaMP.ROI_Intensity((Starttime:EndTime),GCaMP_num);
                    Norm_Const=mean(mv_mpP(1).BehaviorChunkData.GCaMP.Run_Contra_Puff_ROI_Intensity(Run_Contra_GCaMP,(1:NormTime),GCaMP_num),2);
                    mv_mpP(1).BehaviorChunkData.GCaMP.Norm_Run_Contra_Puff_ROI_Intensity(Run_Contra_GCaMP,:,GCaMP_num)=((mv_mpP(1).BehaviorChunkData.GCaMP.Run_Contra_Puff_ROI_Intensity(Run_Contra_GCaMP,:,GCaMP_num)-Norm_Const)/Norm_Const)*100;
                else
                    Still_Contra_GCaMP=Still_Contra_GCaMP+1;
                    mv_mpP(1).BehaviorChunkData.GCaMP.Still_Contra_Puff_ROI_Intensity(Still_Contra_GCaMP,:,GCaMP_num)=mv_mpP(1).GCaMP.ROI_Intensity((Starttime:EndTime),GCaMP_num);
                    Norm_Const=mean(mv_mpP(1).BehaviorChunkData.GCaMP.Still_Contra_Puff_ROI_Intensity(Still_Contra_GCaMP,(1:NormTime),GCaMP_num),2);
                    mv_mpP(1).BehaviorChunkData.GCaMP.Norm_Still_Contra_Puff_ROI_Intensity(Still_Contra_GCaMP,:,GCaMP_num)=((mv_mpP(1).BehaviorChunkData.GCaMP.Still_Contra_Puff_ROI_Intensity(Still_Contra_GCaMP,:,GCaMP_num)-Norm_Const)/Norm_Const)*100;
                end
            end
        end
    end
    for p=1:size(mv_mpP(1).BehaviorTimes.Ctrl_Puff_Times,2)
        Starttime=floor((((mv_mpP(1).BehaviorTimes.Ctrl_Puff_Times(p)-(WinPad*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        EndTime=floor((((mv_mpP(1).BehaviorTimes.Ctrl_Puff_Times(p)+(WinFollow*str2double(mv_mpP(1).Analog_Freq(1:(end-2))))))/str2double(mv_mpP(1).Analog_Freq(1:(end-2))))*mv_mpP(1).Header.Frame_Rate);
        if (Starttime-EndTime)~=217
            EndTime=Starttime+217;
        end
        if Starttime>=1
            if EndTime<=size(mv_mpP(vessel_num).GCaMP.ROI_Intensity,1)
                mv_mpP(1).BehaviorChunkData.GCaMP.Ctrl_Puff_ROI_Intensity(p,:,GCaMP_num)=mv_mpP(1).GCaMP.ROI_Intensity((Starttime:EndTime),GCaMP_num);
                puffwin=(Starttime:EndTime);
                LocoSeg=ismember(puffwin,RunTime);
                if max(LocoSeg)==1
                    Run_Ctrl_GCaMP=Run_Ctrl_GCaMP+1;
                    mv_mpP(1).BehaviorChunkData.GCaMP.Run_Ctrl_Puff_ROI_Intensity(Run_Ctrl_GCaMP,:,GCaMP_num)=mv_mpP(1).GCaMP.ROI_Intensity((Starttime:EndTime),GCaMP_num);
                    Norm_Const=mean(mv_mpP(1).BehaviorChunkData.GCaMP.Run_Ctrl_Puff_ROI_Intensity(Run_Ctrl_GCaMP,(1:NormTime),GCaMP_num),2);
                    mv_mpP(1).BehaviorChunkData.GCaMP.Norm_Run_Ctrl_Puff_ROI_Intensity(Run_Ctrl_GCaMP,:,GCaMP_num)=((mv_mpP(1).BehaviorChunkData.GCaMP.Run_Ctrl_Puff_ROI_Intensity(Run_Ctrl_GCaMP,:,GCaMP_num)-Norm_Const)/Norm_Const)*100;
                else
                    Still_Ctrl_GCaMP=Still_Ctrl_GCaMP+1;
                    mv_mpP(1).BehaviorChunkData.GCaMP.Still_Ctrl_Puff_ROI_Intensity(Still_Ctrl_GCaMP,:,GCaMP_num)=mv_mpP(1).GCaMP.ROI_Intensity((Starttime:EndTime),GCaMP_num);
                    Norm_Const=mean(mv_mpP(1).BehaviorChunkData.GCaMP.Still_Ctrl_Puff_ROI_Intensity(Still_Ctrl_GCaMP,(1:NormTime),GCaMP_num),2);
                    mv_mpP(1).BehaviorChunkData.GCaMP.Norm_Still_Ctrl_Puff_ROI_Intensity(Still_Ctrl_GCaMP,:,GCaMP_num)=((mv_mpP(1).BehaviorChunkData.GCaMP.Still_Ctrl_Puff_ROI_Intensity(Still_Ctrl_GCaMP,:,GCaMP_num)-Norm_Const)/Norm_Const)*100;
                end
            end
        end
    end
end
end
end

function [mv_mpP]=PenetratingVesselAreaTiff(mv_mpP,ImageData,Vessel_Channel, plot_every,maxframe)
% measures the area using an approximation of FWHM at many angles using
% thresholding in radon space
% thefile: tif file to be read
% plot_every:number of fraems to skip before plotting an example image
% maxfram: last frame to measure
nframes=size(ImageData.(Vessel_Channel),3);
the_angles=1:1:180;
rtd_threhold=.5;%0.5 'half max'
irtd_threhold=.2;%0.2 gets rid of the noise from the inverse radon
area=-ones(min(nframes,maxframe),1);
radoncontours=cell(min(nframes,maxframe),1);

fft_first_frame=fft2(ImageData.(Vessel_Channel)(:,:,1)); % *

for mv=1:length(mv_mpP)
    for f=1:min(nframes,maxframe)
        %raw_hold_image=double(sum(imread(thefile,'tif', 'Index',f),3));
        raw_hold_image=double(sum(ImageData.(Vessel_Channel)(:,:,f),3));
        
        fft_raw_hold_frame=fft2(raw_hold_image); % * 
        [mv_mpP(1).pixel_shift(:,f), Greg]=dftregistration(fft_first_frame,fft_raw_hold_frame,1); % *
        
%         hold_image=raw_hold_image(round(mv_mpP(mv).Vessel.box_position.xy(2):mv_mpP(mv).Vessel.box_position.xy(2)+mv_mpP(mv).Vessel.box_position.xy(4)),...
%             round( mv_mpP(mv).Vessel.box_position.xy(1):mv_mpP(mv).Vessel.box_position.xy(1)+mv_mpP(mv).Vessel.box_position.xy(3)));
        [Ordered, Index]=sort(mv_mpP(mv).Vessel.box_position.xy);
        hold_image=raw_hold_image(round(Ordered(1,1):Ordered(4,1)),...
            round( Ordered(1,2):Ordered(4,2)));
        
        hold_image=hold_image-mean(hold_image(:));
        radon_hold_image=radon(hold_image,the_angles);
        
        for k=1:length(the_angles)
            % normalize so that for each agle, the transfomed image is between
            % 0 and 1
            radon_hold_image(:,k)=radon_hold_image(:,k)-min(radon_hold_image(:,k));
            radon_hold_image(:,k)=radon_hold_image(:,k)/max(radon_hold_image(:,k));
            % find the peak of the projection at this angle
            [maxpoint(k),maxpointlocation(k)]=max(radon_hold_image(:,k));
            % threshold at half maximum
            try
                [~,min_edge(k)]=max(find(radon_hold_image(1:maxpointlocation(k),k)<rtd_threhold));
                [~,max_edge(k)]=max(find(radon_hold_image(maxpointlocation(k)+1:end,k)>rtd_threhold));
            catch
                min_edge(k)=0;
                max_edge(k)=0;
            end
            
            radon_hold_image(1:min_edge(k),k)=0;
            radon_hold_image((maxpointlocation(k)+max_edge(k)):end,k)=0;
        end
        % transform the threhlded image back to an x,y image
        irtd_norm=iradon(double(radon_hold_image>rtd_threhold*max(radon_hold_image(:))),(the_angles),'linear','Hamming');
        %draw boundaries
        [cc,l]=bwboundaries(irtd_norm>irtd_threhold*max(irtd_norm(:)));
        %find the biggest area
        numPixels = cellfun(@length,cc);
        [~,idx] = max(numPixels);
        %figure(44)
        %find the nuber of pixes enclosed by the largest area
        area_filled=regionprops(l,'FilledArea','Image','FilledImage');
        area(f)=length(find(area_filled(idx).FilledImage));
        radoncontours{f}=contourc(double(irtd_norm(1:end,1:end)),double([irtd_threhold irtd_threhold]*max(irtd_norm(:))));
        
        % plot the raw image, the inverese of the radon transfoermed imaage, and the contours
        mv_mpP(mv).Vessel.pixel_area=area;
        mv_mpP(mv).Vessel.radoncontours=radoncontours;
        mv_mpP(mv).Vessel.area=mv_mpP(mv).Vessel.pixel_area*mv_mpP(mv).Xfactor*mv_mpP(mv).Xfactor;
        [hold_hist,d]=hist(mv_mpP(mv).Vessel.area,0:10:10000);
        [~,max_d]=max(hold_hist);
        mv_mpP(mv).Vessel.modal_fixed_area=d(max_d);
        mv_mpP(mv).Vessel.diameter=2*sqrt(mv_mpP(mv).Vessel.area/(pi));
        mv_mpP(mv).Vessel.modal_fixed_diameter=2*sqrt(mv_mpP(mv).Vessel.modal_fixed_area/(pi));
    end
    
    f=mv_mpP(1).Header.num_frames;
    if (mod(f,plot_every)==0)
        okpoints=find(radoncontours{f}(1,:)>1);
        use_points=NaN*zeros(size(radoncontours{f}));
        use_points(1,okpoints)=radoncontours{f}(1,okpoints);
        use_points(2,okpoints)=radoncontours{f}(2,okpoints);
        figure
        subplot(221)
        hold off
        imagesc(hold_image)
        axis image
        axis xy
        hold on
        colormap gray
        subplot(221)
        plot(use_points(1,:),use_points(2,:))
        title(['frame=' num2str(f)])
        subplot(222)
        imagesc(irtd_norm)
        axis image
        axis xy
    end
    subplot(2,1,2)
    plot(1:f,mv_mpP(mv).Vessel.pixel_area(1:f),'ro')
    saveas(gcf,fullfile(['rawDiameter_P_radon_' mv_mpP.Header.Filename(1:(length(mv_mpP.Header.Filename)-4)) '_' num2str(mv)]),'fig');
    close
end
end
        
        
    
    
    
    