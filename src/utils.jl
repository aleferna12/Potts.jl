allsame(x) = all(y -> y == first(x), x)

ordered(x1, x2) = x2 < x1 ? (x2, x1) : (x1, x2)

function fullyconcrete(type::Type)
    stack = Type[type]
    while !isempty(stack)
        t = pop!(stack)
        if !isconcretetype(t)
            return false
        end
        paramtypes = filter(pt -> pt isa Type, collect(t.parameters))
        append!(stack, fieldtypes(t), paramtypes)
    end
    true
end