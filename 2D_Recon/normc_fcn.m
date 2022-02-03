function nm = normc_fcn(m)

% This code is part of the following work which has been submitted to IEEE Transactions on Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade  

nm = sqrt(m.^2 ./ sum(m.^2)) .* sign(m);
end