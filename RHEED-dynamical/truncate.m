function [crystalout] = truncate(crystal, h, t, dir)
% truncate
% Simple function to truncate the crystal structure at a specified height
% and thickness.
%
% Inputs:
% crystal       Input crystal structure
% h             Height at which to truncate crystal
% t             Thickness of truncated crystal
% dir           String specifying whether to take a slice of thickness t
%               above (up) or below (down) height h.
%
% Outputs:
% crystalout    Output crystal structure

% Initialize variables.
crystalout = crystal;
atoms0 = crystal.atoms;
r0 = crystal.r;

% Make an array of empty strings and empty coordinates.
atoms = strings(size(atoms0));
x = NaN(size(atoms0)); y = x; z = x;

switch dir
    case 'up'
        % Loop on all atoms and coordinates in the input list.
        for i=1:length(atoms0)
            % Find the list entries with r between h and h+t.
            if (r0(i,3) >= h) && (r0(i,3) <= h+t)
                atoms(i) = atoms0(i);
                x(i) = r0(i,1);
                y(i) = r0(i,2);
                z(i) = r0(i,3);
            end
        end

        % Remove empty values.
        atoms = atoms(~ismember(atoms,''));
        x = x(~isnan(x));
        y = y(~isnan(y));
        z = z(~isnan(z));
        r = [x, y, z];

    case 'down'
        % Loop on all atoms and coordinates in the input list.
        for i=1:length(atoms0)
            % Find the list entries with r between h-t and h.
            if (r0(i,3) >= h-t) && (r0(i,3) <= h)
                atoms(i) = atoms0(i);
                x(i) = r0(i,1);
                y(i) = r0(i,2);
                z(i) = r0(i,3);
            end
        end

        % Remove empty values.
        atoms = atoms(~ismember(atoms,''));
        x = x(~isnan(x));
        y = y(~isnan(y));
        z = z(~isnan(z));
        r = [x, y, z];

    otherwise
        error('Invalid direction specified.');
end

% Re-reference the highest coordinate position as z = 0.
z0 = max(r(:,3));
r(:,3) = r(:,3) - z0;

% Pack up the output crystal structure.
crystalout.atoms = atoms;
crystalout.r = r;

end