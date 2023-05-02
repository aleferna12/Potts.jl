@enum Tau veg=1 div

mutable struct AttrSet
    sigma::Int
    tau::Tau
    area::Int
    center::Pos
end

mutable struct Genome{}
    sigma::Int
end

struct Cells
    attrsets::Vector{AttrSet}
    genomes::Vector{Genome}
    matrix::Matrix{Int}
    edgeset::Set{Edge}
end
Cells(fieldsize::Integer) = Cells(AttrSet[], Genome[], createcellmatrix(fieldsize), Set{Edge}())
Base.length(cells::Cells) = length(cells.attrsets)
Base.lastindex(cells::Cells) = lastindex(cells.attrsets)
Base.size(cells::Cells) = size(cells.matrix)
Base.size(cells::Cells, dim::Integer) = size(cells.matrix, dim)
getsigma(cells, pos::MatrixPos) = cells.matrix[pos]
setsigma(cells, pos::MatrixPos, val) = cells.matrix[pos] = val
gettau(cells, sigma::Integer) = cells.attrsets[sigma].tau
getarea(cells, sigma::Integer) = cells.attrsets[sigma].area

function getadhesion(cells, sigma1, sigma2, adhesiontable, adhesionmedium)::Float64
    if sigma1 == sigma2
        return 0.
    end
    if sigma1 == 0 || sigma2 == 0
        return adhesionmedium
    end
    adhesiontable[Int(gettau(cells, sigma1)), Int(gettau(cells, sigma2))]
end

function createcellmatrix(fieldsize)
    m = zeros(Int, fieldsize, fieldsize)
    for i in 1:fieldsize
        m[i, 1] = -1
        m[i, fieldsize] = -1
        m[1, i] = -1
        m[fieldsize, i] = -1
    end
    m
end

function spawncell!(cells::Cells, cell_len::Integer)
    sigma = lastindex(cells) + 1
    spawncenter = getrandompos(cells.matrix, cell_len)
    area, bb = initcellpositions!(cells.matrix, sigma, spawncenter, cell_len)
    attrset = AttrSet(
        sigma,
        veg,
        area,
        getcenter(bb)
    )
    push!(cells.attrsets, attrset)
    union!(cells.edgeset, getedges(cells, bb))
end

"Initialize a cell on a matrix and return its area and bounding box."
function initcellpositions!(matrix::Matrix, sigma::Integer, center::MatrixPos, cell_len::Integer)
    area = 0
    bb = BoundingBox(center, center)
    rangebb = BoundingBox(
        MatrixPos(center.x - floor(cell_len / 2), 
                  center.y - floor(cell_len / 2)),
        MatrixPos(center.x + ceil(cell_len / 2) - 1, 
                  center.y + ceil(cell_len / 2) - 1)
    )
    for pos in iterpositions(rangebb)
        if matrix[pos.x, pos.y] == 0
            matrix[pos.x, pos.y] = sigma
            area += 1
            addpos!(bb, pos)
        end
    end
    area, bb
end

function getedges(cells, bb::BoundingBox)
    edges = Set{Edge}()
    for pos in iterpositions(bb)
        sigmapos = getsigma(cells, pos)
        for neigh in moore_neighbors(pos)
            sigmaneigh = getsigma(cells, neigh)
            if sigmapos > -1 && sigmaneigh > -1 && sigmapos != sigmaneigh
                addedges!(edges, pos, neigh)
    end end end
    edges
end