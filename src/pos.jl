"A position on the simulation matrix."
struct Pos{T<:Number}
    x::T
    y::T
end
getx(pos::Pos) = pos.x
gety(pos::Pos) = pos.y
Base.round(postype::Type, pos::Pos) = postype(round(getx(pos)), round(gety(pos)))
Base.round(pos::Pos) = round(typeof(pos), pos)
Base.convert(postype::Type{<:Pos}, pos::NTuple{2}) = postype(pos[1], pos[2])

"Returns true if the position lies above a line specified by 'm' and 'n' and false otherwise."
function isabove(pos::Pos, m, n)
    y = m * getx(pos) + n
    gety(pos) - y > 0
end

"Splits the positions in two vectors, those above and below the line defined by 'm' and 'n'."
function split(positions, m, n)
    above, below = MatrixPos[], MatrixPos[]
    for pos in positions
        if isabove(pos, m, n)
            push!(above, pos)
        else
            push!(below, pos)
    end end
    above, below
end

function split(positions)
    center = getcenter(positions)
    m = tan(rand() * 2π)
    n = gety(center) - m * getx(center)
    split(positions, m, n)
end

function getcenter(positions)
    centerx, centery = 0., 0.
    for pos in positions
        centerx += getx(pos)
        centery += gety(pos)
    end
    Pos(centerx / length(positions), centery/length(positions))
end

"Often in CPM simulations positions are represented as matrix indices, so we provide this alias definition."
const MatrixPos = Pos{Int}
Base.getindex(matrix::Matrix, pos::MatrixPos) = matrix[pos.x, pos.y]
Base.getindex(matrix::Matrix, positions::MatrixPos...) = [matrix[pos] for pos in positions]
Base.setindex!(matrix::Matrix, value, pos::MatrixPos) = matrix[pos.x, pos.y] = value
Base.setindex!(matrix::Matrix, value, positions::MatrixPos...) = for pos in positions matrix[pos] = value end
Base.checkbounds(::Type{Bool}, matrix::Matrix, index::MatrixPos) = checkbounds(Bool, matrix, CartesianIndex(getx(index), gety(index)))
Base.convert(::Type{MatrixPos}, pos::Pos) = MatrixPos(getx(pos), gety(pos))
moore_neighbors(pos::MatrixPos) = [MatrixPos(pos.x - 1, pos.y),
                                   MatrixPos(pos.x + 1, pos.y),
                                   MatrixPos(pos.x, pos.y - 1),
                                   MatrixPos(pos.x, pos.y + 1),
                                   MatrixPos(pos.x - 1, pos.y - 1),
                                   MatrixPos(pos.x - 1, pos.y + 1),
                                   MatrixPos(pos.x + 1, pos.y - 1),
                                   MatrixPos(pos.x + 1, pos.y + 1)]

function getrandompos(matrix::Matrix, borderpadding::Integer=0)
    xrange = range(1 + borderpadding, size(matrix)[1] - borderpadding)
    yrange = range(1 + borderpadding, size(matrix)[2] - borderpadding)
    MatrixPos(rand(xrange), rand(yrange))
end

const Edge = Pair{MatrixPos, MatrixPos}

"Removes the edges specified by 'pos1' and 'pos2' from the edge set. Throws error if positions are not in the set."
function removeedges!(edgeset, pos1, pos2)
    pop!(edgeset, Edge(pos1, pos2))
    pop!(edgeset, Edge(pos2, pos1))
end

"Adds the edges specified by 'pos1' and 'pos2' to the edge set.
Returns true if the positions were not already in the set and false otherwise."
function addedges!(edgeset, pos1, pos2)
    edge1 = Edge(pos1, pos2)  # We only check for the first one for efficiency
    if edge1 in edgeset
        false
    else
        push!(edgeset, edge1, Edge(pos2, pos1))
        true
end end
