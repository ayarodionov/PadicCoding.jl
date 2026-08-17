# Demo: compress a text with the p-adic arithmetic coder (P = 2 and P = 3).
#
#   julia --project demo.jl [file]
#
# Without an argument a built-in sample text is used.  The model is static:
# byte frequencies are counted in a first pass (two-pass coding).

using PadicCoding

text = if isempty(ARGS)
    repeat("For more than forty years I've been speaking in prose " *
           "without even knowing it! ", 40)
else
    read(ARGS[1], String)
end

data = Vector{UInt8}(text)
msg = Int.(data) .+ 1                 # bytes 0..255 -> symbols 1..256

freqs = ones(Int, 256)                # +1 so every symbol has a weight
for b in data
    freqs[b + 1] += 1
end

entropy = let p = filter(>(0), (freqs .- 1) ./ length(data))
    -sum(x -> x * log2(x), p)
end

for P in (2, 3)
    m = PadicModel(freqs; P = P, N = P == 2 ? 40 : 26)
    bits = padic_encode(m, msg)
    @assert padic_decode(m, bits) == msg
    outbits = length(bits) * log2(P)  # one P-bit carries log2(P) bits
    println("P = $P: $(length(data)) bytes -> $(length(bits)) $P-bits ",
            "($(round(outbits / 8, digits = 1)) bytes, ",
            "$(round(outbits / length(data), digits = 3)) bits/byte; ",
            "empirical entropy $(round(entropy, digits = 3)) bits/byte)")
end
