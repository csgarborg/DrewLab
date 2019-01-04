function []=Proc2PMDFData(Thefiles,Imaging_channels,Analog_Channel_Num,temporal_filter_size,ChunkSize,PixelOffset,varargin)
if exist('Thefiles','var')==0
Thefiles=uigetfile('*.MDF','MultiSelect','on'); %make list of MDF files to be filtered and image registered
end
if exist('Imaging_channels','var')==0
Imaging_channels=input('How many channels were imaged?');
end
if exist('temporal_filter_size','var')==0
temporal_filter_size=input('How many frames for median filter?');
end
if exist('ChunkSize','var')==0
ChunkSize=9270; %number of frames to be exported to memory at same time %9288
end
if exist('Analog_Channel_Num','var')==0
    Analog_Channel_Num=input('How many analog channels were recorded?');
end
if exist('PixelOffset','var')==0
    PixelOffset=input('How many pixels along y axis does image need shifted?');
end
for q=1:numel(Thefiles)
    close all;
    fname=Thefiles{q};
    MotionCorrectMedFiltMDF_AsymmetricScan_002(fname,Imaging_channels,Analog_Channel_Num,temporal_filter_size,ChunkSize,PixelOffset);
end
end