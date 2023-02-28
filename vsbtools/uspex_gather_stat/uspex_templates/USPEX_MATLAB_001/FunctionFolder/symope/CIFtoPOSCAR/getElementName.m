function [ element_name ] = getElementName( raw_string )
% This function is to read the element name from a string
% ex: raw_string ='Cu12' element_name = 'Cu'

    %disp(raw_string);
    element_name = '';
    if length(raw_string)==1
       element_name = raw_string;
    else
       for i=1:length(raw_string)
          if (isletter(raw_string(i)))
           element_name = strcat(element_name, raw_string(i));
          else 
              break;
          end
       end
end


