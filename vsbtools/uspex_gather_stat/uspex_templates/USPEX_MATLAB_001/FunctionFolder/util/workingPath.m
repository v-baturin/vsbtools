function [homePath, USPEXPath] = workingPath(homePath) 

if nargin == 0
	homePath = pwd;
end
	
USPEXPath = homePath;
uspexmode = getenv('UsPeXmOdE');
if strcmp(uspexmode,'exe')
    USPEXPath = getenv('USPEXPATH');
end

