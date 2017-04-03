function [var1_out,var2_out,stats] = par_vec(var1,var2)
%[var1_out,var2_out] = par_vec(var1,var2)
%Function for vectorizing 2 variables to work in parfor loop

stats = struct;
len_col = length(var2);
len_row = length(var1);
par_iter = len_row*len_col;
var1_out = nan(1,par_iter); 
var2_out = var1_out;
for a=1:len_col
    var2_out((a-1)*len_row + 1:a*len_row) = var2(a)*ones(1,len_row);
    var1_out((a-1)*len_row + 1:a*len_row) = var1;
end

stats.len_col = len_col;
stats.len_row = len_row;
stats.iter = par_iter;