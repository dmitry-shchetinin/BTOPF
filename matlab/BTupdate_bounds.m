function bounds=BTupdate_bounds(bounds, bounds_temp)
%update bounds that have been tightened
bounds.min=max(bounds_temp.min,bounds.min);
bounds.max=min(bounds_temp.max,bounds.max);
end

