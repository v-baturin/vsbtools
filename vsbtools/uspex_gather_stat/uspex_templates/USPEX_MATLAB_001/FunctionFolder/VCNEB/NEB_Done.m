function NEB_Done()


global ORG_STRUC

cd( ORG_STRUC.resFolder );

unixCmd(['echo -e " " >>  VCNEBReports' ]);
unixCmd(['echo -e "=====================================  VCNEB  DONE ===================================="   >> VCNEBReports']);
unixCmd(['echo -e " " >>  VCNEBReports' ]);

