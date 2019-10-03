function Distortion = distortion_4(f , y_1 , y_2 , codebook , delta , Pr_z , T)
summation = 0 ;
parfor x_3 = 1 : 2
    for x_4 = 1 : 2
        x = (x_3 - 1) * 2 + x_4 ;
        for y_3 = 1 : 2
            for y_4 = 1 : 2
                y = (y_3 - 1) * 2 + y_4 ;
                u_index = find (T(: , 3) == x) ;
                for u_i = 1 : length(u_index)
                    x_prime = T(u_index(u_i) , 2) ; 
                    if (mod(x_prime - 1 , 2) == 0 )
                        x_2 = 1 ; 
                    else
                        x_2 = 2 ; 
                    end 
                    summation = summation + Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                        * Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 )...
                        * delta * f(u_index(u_i)).* (T(u_index(u_i) , 1) - codebook(y)) .^ 2 ;
                end
            end
        end
    end
end
Distortion = summation ;
end