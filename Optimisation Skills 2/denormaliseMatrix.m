% denormaliseMatrix
function denormalisedMatrix = denormaliseMatrix(matrix, minBounds, maxBounds)
    denormalisedMatrix = minBounds + matrix .* (maxBounds - minBounds);
end