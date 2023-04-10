mutable struct BoundingBox
    mini::Int
    maxi::Int
    minj::Int
    maxj::Int
end

"""Add a position to a BoundingBox by enlarging it."""
function addpos!(bb::BoundingBox, i::Integer, j::Integer)
    addposi!(bb, i)
    addposj!(bb, j)
end

# When updating the hamiltonian we can speed up the update process by skipping unecessary min and max calls 
function addposi!(bb::BoundingBox, i::Integer)
    bb.mini = min(bb.mini, i)
    bb.maxi = max(bb.maxi, i)
end

function addposj!(bb::BoundingBox, j::Integer)
    bb.minj = min(bb.minj, j)
    bb.maxj = max(bb.maxj, j)
end

"""Remove a pos from a BoundingBox by decreasing the corresponding edge.

The BoundingBox might end up being larger than necessary.
"""
function removepos!(bb::BoundingBox, i::Integer, j::Integer)
    removeposi!(bb, i)
    removeposj!(bb, j)
end

function removeposi!(bb::BoundingBox, i::Integer)
    if bb.mini == i
        bb.mini += 1
    elseif bb.maxi == i
        bb.maxi -= 1
    end
end

function removeposj!(bb::BoundingBox, j::Integer)
    if bb.minj == j
        bb.minj += 1
    elseif bb.maxj == j
        bb.maxj -= 1
    end
end

function iteratepos(boundingbox::BoundingBox)
    ((i, j) for i in boundingbox.mini:boundingbox.maxi, j in boundingbox.minj:boundingbox.maxj)
end