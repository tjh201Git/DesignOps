% normalise values
function normalisedMatrix = normaliseMatrix(matrix, minBounds, maxBounds)
    normalisedMatrix = (matrix - minBounds) ./ (maxBounds - minBounds);
end