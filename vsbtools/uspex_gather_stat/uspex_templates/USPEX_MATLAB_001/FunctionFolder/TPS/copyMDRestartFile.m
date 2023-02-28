function error = copyMDRestartFile(restartFile, varargin)



global TPS_STRUC
global ORG_STRUC

%-- If we dont specify the resource file, we copy it from result folder
if isempty(restartFile)
    restartFile = [ORG_STRUC.homePath '/' TPS_STRUC.resFolder '/' ORG_STRUC.restartFile];
end

%-- The new name we want to use for the new MD restart file
 
if ~isempty(varargin)
    newName = varargin{1};
else
    newName = ORG_STRUC.restartFile;
end


%
if ~exist(restartFile,'file')
    error = 2;
else
    [error, nothing] = unix(['cp ' restartFile  '  ' newName ]);
end

if error > 0
   disp( ['Failed to copy  the ' restartFile ' file. Please check your input files ....'] );
   %save ([ ORG_STRUC.homePath '/' ORG_STRUC.resFolder '/ERROR_param.txt'],'error')
   quit();
end
