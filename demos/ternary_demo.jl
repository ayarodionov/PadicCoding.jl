# Demo: is compression better on a ternary (trit) machine than on a binary one?
#
#   julia --project demos/ternary_demo.jl
#
# Each source is coded three ways, in base 2 and in base 3:
#   * the p-adic arithmetic coder (P = 2 and P = 3),
#   * a D-ary Huffman code (D = 2 and D = 3, one codeword per symbol).
# Output sizes are converted to bits (1 trit = log2(3) = 1.585 bits) so the
# two bases can be compared.  The last column is the "radix economy" hardware
# cost of storing the output: digits * radix (states that must be built).

using PadicCoding, Random, Printf

# ---- D-ary Huffman code lengths -------------------------------------------
function huffman_lengths(p::Vector{Float64}, D::Int)
    n = length(p)
    pad = (D - 1 - (n - 1) % (D - 1)) % (D - 1)     # dummy symbols so the tree is full
    nodes = [(p[i], [i]) for i in 1:n]
    append!(nodes, [(0.0, Int[]) for _ in 1:pad])
    len = zeros(Int, n)
    while length(nodes) > 1
        sort!(nodes; by = first)
        merged = nodes[1:D]
        nodes = nodes[D+1:end]
        leaves = reduce(vcat, last.(merged))
        len[leaves] .+= 1
        push!(nodes, (sum(first.(merged)), leaves))
    end
    return len
end

# ---- one experiment ---------------------------------------------------------
function run(name, freqs, msg)
    n = length(msg)
    p = freqs ./ sum(freqs)
    H = -sum(log2(p[s]) for s in msg) / n          # ideal rate for this model

    ac = Dict{Int,Int}()                           # P -> number of P-digits
    for P in (2, 3)
        m = PadicModel(freqs; P = P, N = P == 2 ? 40 : 26)
        out = padic_encode(m, msg)
        @assert padic_decode(m, out) == msg
        ac[P] = length(out)
    end
    hl = Dict(D => sum(huffman_lengths(p, D)[s] for s in msg) for D in (2, 3))

    b(P, k) = k * log2(P) / n                      # digits -> bits/symbol
    @printf("\n%s   (%d symbols, ideal %.4f bits/sym)\n", name, n, H)
    @printf("  %-22s %10s %10s %12s\n", "coder", "digits", "bits/sym", "cost=d*radix")
    for (lbl, P, k) in (("arithmetic, binary", 2, ac[2]), ("arithmetic, ternary", 3, ac[3]),
                        ("Huffman,    binary", 2, hl[2]), ("Huffman,    ternary", 3, hl[3]))
        @printf("  %-22s %10d %10.4f %12d\n", lbl, k, b(P, k), k * P)
    end
end

Random.seed!(1)
N = 30_000
draw(p, n) = [findfirst(>=(rand()), cumsum(p)) for _ in 1:n]

# 1. Uniform over 3 symbols: the "native" ternary source.
run("A. uniform {1/3,1/3,1/3}", [1000, 1000, 1000], draw([1/3, 1/3, 1/3], N))

# 2. Dyadic source: the "native" binary source.
run("B. dyadic {1/2,1/4,1/4}", [2000, 1000, 1000], draw([1/2, 1/4, 1/4], N))

# 3. Uniform random bytes (incompressible).
run("C. uniform bytes", fill(1000, 256), rand(1:256, N))

# 4. English-like text, static byte model.
let text = repeat("For more than forty years I've been speaking in prose " *
                  "without even knowing it! ", 400)[1:N]
    msg = Int.(Vector{UInt8}(text)) .+ 1
    freqs = ones(Int, 256)
    for s in msg; freqs[s] += 1; end
    run("D. text bytes", freqs, msg)
end

println("\nRadix economy b/ln(b) (hardware cost per unit of information):")
for b in (2, ℯ, 3, 4)
    @printf("  b = %-6.3f  %.4f\n", b, b / log(b))
end
@printf("  ternary vs binary: %.2f%% cheaper\n", 100 * (1 - (3 / log(3)) / (2 / log(2))))
