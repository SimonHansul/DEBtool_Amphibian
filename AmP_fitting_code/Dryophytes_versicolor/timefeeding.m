% function with the feeding sewitches from Beachy et al. 1999
function [f] = timefeeding(t, f1, f2, switch1)
    if t < 15+4
        f = f1;
    else
       switch switch1
           case 1
              f=f2; % to be improved
           case 0 
              if t>=34+4
                  f=f2; % to be improved
              else 
                  f=f1;
              end
           case -1
               f = f1;
        end
    end
end