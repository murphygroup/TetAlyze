
q = 1;
T = zeros(1,18);

Str = fileread('Summary.txt');
Keys = ["Anterior Area (um^2):"; "# Anterior BBs:"; "Medial Area (um^2):"; "# Medial BBs:"; "Posterior Area (um^2):";...
    "# Posterior BBs:"; "# assigned BBs:"; "# BB rows:"; "Avergae number of BBs per row:"; "Cell height(um):"; "Cell width(um):";...
    "Cell volume(um^3):"; "Average neighbor BB pairwise distance(um) (anterior, medial, posterior):"; "Average BB row pairwise distance(um) (anterior, medial, posterior):";];

for p=1:length(Keys)
    Key = convertStringsToChars(Keys(p));
    Index = strfind(Str, Key);
    if p==13 || p==14
        Value = sscanf(Str(Index(1) + length(Key):end), '%g, %g, %g', 3);
        for j=1:3
            T(q)=Value(j);
            q = q+1;
        end
    else
        Value = sscanf(Str(Index(1) + length(Key):end), '%f', 1);
        T(q)=Value;
        q = q+1;
    end
end
