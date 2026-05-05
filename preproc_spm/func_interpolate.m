function func_interpolate(reconFileRaw, badVols, output)

    % Warning: MATLAB indexes at 1
    % convert to numeric and covert to array (')
    badVolNums = str2double(split(badVols))' ;

    % TO DO: pass badVolNums as space delim string and sep

    if ~isempty(badVolNums)

        disp('Interpolating bad volumes...')
        disp(badVolNums)

        img = niftiread(reconFileRaw);
        hdr = niftiinfo(reconFileRaw);

        keep_vols = setdiff(1:size(img, 4), badVolNums) ;
        for v = badVolNums % exclude

            disp(v) ;

            % the closest previous non-excluded volume
            before_vol = max(keep_vols(keep_vols < v)) ;

            % the closest following non-excluded volume
            after_vol = min(keep_vols(keep_vols > v)) ;

            if isempty(before_vol) && ~isempty(after_vol)
                 % no volumes before, so copy closest following non-excluded
                 % volume
                img(:, :, :, v) = img(:, :, :, after_vol);
            elseif isempty(after_vol) && ~isempty(before_vol)
                % no volumes after, so copy closest previous non-excluded
                % volume
                img(:, :, :, v) = img(:, :, :, before_vol);
            elseif isempty(before_vol) && isempty(after_vol)
                % this should not happen
                error('Error while replacing bad volumes');
            else
                % weighted average of closest before and after non-excluded
                % volumes
                f = (v - before_vol) / (after_vol - before_vol);
                img(:, :, :, v) = (1 - f) .* img(:, :, :, before_vol) + f .* img(:, :, :, after_vol);
            end
        end

        % Overwrite raw with interpolated 4D file...
        niftiwrite(feval(hdr.Datatype, img), output, hdr)

    else

        disp('No bad volumes to interpolate, so skipping bad volume interpolation.')

    end

end