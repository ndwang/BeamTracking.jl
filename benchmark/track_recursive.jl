@inline function track_lattice_fused(p::PhaseSpace6D, lattice::Tuple{})
    return p
end

@inline function track_lattice_fused(p::PhaseSpace6D, lattice)
    # Track through the first element
    p_after_first = track!(p, first(lattice))
    
    # Recurse on the tail of the tuple
    return track_lattice_fused(p_after_first, Base.tail(lattice))
end

function unrolled_track(p, lattice)
    p1 = track!(p,  lattice[1])
    p2 = track!(p1, lattice[2])
    p3 = track!(p2, lattice[3])
    p4 = track!(p3, lattice[4])
    return p4
end