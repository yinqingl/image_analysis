% process images using imagej within matlab
% using package mij
%
% Yinqing Li
% yinqingl@csail.mit.edu
% 2016

%%
%start mij
% run('C:\\fiji-win64\\Fiji.app\\scripts\\ImageJ.m')
addpath('C:\\fiji-win64\\Fiji.app\\scripts')
ImageJ

%use imagej api
imagej_w = ij.ImageJ;    %image window
imagej_ij = ij.IJ;    %imagej macro

% MIJ.start

%load image in imagej
% MIJ.run('Open...','path=C:\\Users\\yinqing\\Desktop\\20170125_m1-2-5_4-1_1.nd2.tif');
imagej_ij.runMacro('open("C:\\Users\\yinqing\\Desktop\\20170125_m1-2-5_4-1_1.nd2.tif");');

%imagej process
imagej_macro_list = {
    'run("Channels Tool...");'
    'Stack.setDisplayMode("color");'
    'Stack.setChannel(CHANNEL);',
    'run("COLOR");' %Blue, Green, Red, Magenta
    'run("Brightness/Contrast...");'
    'run("Enhance Contrast", "saturated=0.35");'
    'run("Apply LUT");'
    'setMinAndMax(MIN_VAL, MAX_VAL);'
    'saveAs("Tiff", "PATH2FILE");'
    'run("Close All");'
};
for i = 1:size(imagej_macro_list,1)
    disp(sprintf('%d, %s',i,imagej_macro_list{i}))
end

imagej_ij.runMacro(str_macro);

%quit imagej
imagej_ij.runMacro('run("Close All");');
imagej_w.quit();


%to quit mij
%{
MIJ.exit
%}

