function [vol]=imreadBF_nd2(id,zplanes,tframes,channel)
%[vol]=imreadBF_nd2(datname,zplanes,tframes,channel)
%
%imports images nd2 format using the BioFormats package
%you can load multiple z and t slices at once, e.g. zplanes=[1 2 5] loads
%first,second and fifth z-slice in a 3D-Stack 
%
%t is used as series
%
%if loading multiple z slices and tframes, everything is returned in one 3D
%Stack with order ZT. Only one channel can be imported at once
%
%use imreadBFmeta_nd2() to get corresponding metadata of the image file
%
% install bfmatlab
% http://downloads.openmicroscopy.org/bio-formats/
%
% Yinqing Li
% yinqingl@csail.mit.edu

% load the Bio-Formats library into the MATLAB environment
% status = bfCheckJavaPath(autoloadBioFormats);
% assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
%     'to the static Java path or add it to the Matlab path.']);

% initialize logging
% javaMethod('enableLogging', 'loci.common.DebugTools', 'INFO');

% Get the channel filler
r = bfGetReader(id, 0);

width = r.getSizeX();
height = r.getSizeY();

channel=channel-1;
zplane=zplanes-1;
tframe=tframes-1;

vol=zeros(height,width,length(zplane)*length(tframe));

zahler = 0;
for s = 1:length(tframe)
    r.setSeries(tframe(s));
    for z=1:length(zplane)
        p = r.getIndex(zplane(z),channel,0) + 1;
        arr = bfGetPlane(r, p); 
        zahler=zahler+1;
        vol(:,:,zahler)=arr;
    end
end
    
    
end
    
    
    
function [result] = versionCheck(v, maj, min)

tokens = regexp(v, '[^\d]*(\d+)[^\d]+(\d+).*', 'tokens');
majToken = tokens{1}(1);
minToken = tokens{1}(2);
major = str2num(majToken{1});
minor = str2num(minToken{1});
result = major > maj || (major == maj && minor >= min);
end    
    
            
            
            
            
            
  
            
            

        
   
   
   
   
   
            
            