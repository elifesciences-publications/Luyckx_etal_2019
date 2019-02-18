function [y] = sigmoid(x, minimum,range, inflection, steepness)
%function [y] = sigmoid(x, minimum, range, inflection, steepness)

    y = minimum + range ./ (1 + exp(-(x-inflection)/steepness));

end