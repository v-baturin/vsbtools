function SuccessRate()

global ORG_STRUC
global POP_STRUC

if ORG_STRUC.abinitioCode(1) > 0                                                 
    good_and_fit = 0;                                                            
    for fit_loop = 1 : length(POP_STRUC.POPULATION)                              
        if POP_STRUC.POPULATION(fit_loop).Enthalpies(ORG_STRUC.conv_till) < 90000                
            good_and_fit = good_and_fit + 1;                                     
        end                                                                      
    end                                                                          
else                                                                             
    good_and_fit = length(POP_STRUC.POPULATION);                                 
end                                                                              
                                                                                 
if good_and_fit < floor(length(POP_STRUC.POPULATION)/3)                          
    USPEXmessage(561, '', 0, 'Error');
elseif good_and_fit/length(POP_STRUC.POPULATION) < ORG_STRUC.bestFrac            
    USPEXmessage(555, '', 0);                                                      
end                                                                              
                                                                                 
POP_STRUC.good_and_fit = good_and_fit;    
