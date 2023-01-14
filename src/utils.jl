
ziptreesyms(iter) = zip(treesyms(iter), iter)
treesyms(x) = treesyms(length(x))
treesyms(i::Integer) = Iterators.flatten((Iterators.repeated('├', i-1), '└'))

function treestyle_string(d::AbstractDict)
    reprf = x -> repr(MIME("text/plain"), x)
    keylength = mapreduce(k->lastindex(reprf(k)), max, keys(d))

    res = mapreduce(*, ziptreesyms(pairs(d))) do (ts, pair)
        ts * " " * lpad(reprf(pair.first), keylength) * " = " * repr(pair.second)*"\n"
    end
end
