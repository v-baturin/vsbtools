function result = callAWK(scriptname, filename, varargin)


global ORG_STRUC

option = '';
if ~isempty(varargin)
   option=varargin{:};
end

[homePath, USPEXPath] = workingPath(ORG_STRUC.homePath);

awkScript=[USPEXPath, '/FunctionFolder/Tool/', scriptname];

%disp(['awk -f ' awkScript '  ' option '  '  filename ]);

[nothing, resultStr] = unix(['awk -f ' awkScript '  ' option '  '  filename ]);
try
   result=str2num(resultStr);
   if isempty(result)          % this condition is necessary for last
      result = resultStr;      % step symmetrization
   end
catch
   result = [];
   USPEXmesage(1010,resultStr,0);
   disp(resultStr);
end
