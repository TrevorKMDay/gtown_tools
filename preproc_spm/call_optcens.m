
function call_optcens = call_optcens(spm_matrix)
    
    % maxNumCompThreads(1) ;

    disp(spm_matrix)

    approach = 1 ; % task
    docalc = 1 ; % censor
    min_removal = 0 ; % min % of data remove - what??
    max_removal = 50 ; % max % of data to remove 
    
    zeropad = 0 ; % zeropad/adapt = do not conform matrices to one another
    adapt = 0 ;
    
    [out_cens, out_int, out_both] = func_optcens(spm_matrix, approach, ...
        docalc, min_removal, max_removal, zeropad, adapt) ;

    disp(out_cens)

end