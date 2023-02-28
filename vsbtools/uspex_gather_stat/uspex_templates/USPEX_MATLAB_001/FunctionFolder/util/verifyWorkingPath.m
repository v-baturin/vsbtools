function verifyWorkingPath(homePath, USPEXPath)


global ORG_STRUC

if ~strcmp(ORG_STRUC.homePath, homePath)
    mesg=[ORG_STRUC.homePath ' --> ' homePath ];
    ORG_STRUC.homePath=homePath;
    USPEXmessage(1011, mesg, 0);
end
ORG_STRUC.USPEXPath=USPEXPath;
