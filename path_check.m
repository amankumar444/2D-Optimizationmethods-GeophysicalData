%Function to create path in any system::

function [path] = path_check(str1, str2func)

    str = str1;
    for i = strlength(str):-1:1
        str_new = str(i);
        if str_new == '/' || str_new == '\'
            break;
        end
    end 
    opts=str2func;
    res = strcat(str,str_new,opts{1});

    for i=2:1:length(opts)
        res = strcat(res,str_new,opts{i});
    end
    path = res;
end