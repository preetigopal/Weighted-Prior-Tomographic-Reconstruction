function nm = normc_fcn(m)
nm = sqrt(m.^2 ./ sum(m.^2)) .* sign(m);
end