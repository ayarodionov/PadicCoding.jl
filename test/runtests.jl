using Test
using Random

using PadicCoding
using PadicCoding: hat, compath, select_path, I128

@testset "p-adic helpers (paper examples)" begin
    # Picture 6 table, P = 2, level 3: path {0,1,1} = p-adic number 6 <-> index 3
    @test hat(I128(3), 2, 3) == 6
    @test hat(I128(6), 2, 3) == 3
    # ^ is idempotent
    for x in 0:15
        @test hat(hat(I128(x), 2, 4), 2, 4) == x
    end
    # p. 12: com_{2,3}(x = 4, y = 6) = ord_2(2) = 1
    @test compath(I128(4), I128(6), 2, 3) == 1
    # p. 16, Picture 8: minimal path point of [g5(4), g10(4)) is g8(4),
    # its path {1,0,0,0} is the p-adic number 1
    @test select_path(I128(5), I128(10), 2, 4, I128(16)) == 1
end

@testset "roundtrip PR-only coder" begin
    # weight interval of the paper: {a:[0,1/2), b:[1/2,3/4), c:[3/4,7/8), d:[7/8,1)}
    m = PadicModel([4, 2, 1, 1] .* 8; P = 2, N = 40)
    for msg in ([1, 2, 1], Int[], [1], [4], [1, 2, 3, 4, 1, 1, 2])
        bits = padic_encode_pr(m, msg)
        @test padic_decode_pr(m, bits) == msg
    end
    rng = MersenneTwister(1)
    for len in (10, 100, 500)
        msg = rand(rng, 1:4, len)
        @test padic_decode_pr(m, padic_encode_pr(m, msg)) == msg
    end
end

@testset "roundtrip PR+AR coder, P = 2" begin
    m = PadicModel([4, 2, 1, 1]; P = 2, N = 32)
    rng = MersenneTwister(2)
    for msg in (Int[], [1], [4], [1, 2, 1], collect(1:4))
        @test padic_decode(m, padic_encode(m, msg)) == msg
    end
    for len in (1, 7, 100, 1000, 5000)
        msg = rand(rng, 1:4, len)
        @test padic_decode(m, padic_encode(m, msg)) == msg
    end
end

@testset "AR rescaling stress (interval straddling 1/2)" begin
    # Shuffled weights as on p. 22: {b:[0,1/4), a:[1/4,3/4), ...} -- a message
    # of only 'a' straddles the centre point, PR never fires, AR must save us.
    m = PadicModel([1, 2, 1]; P = 2, N = 24)
    msg = fill(2, 2000)
    bits = padic_encode(m, msg)
    @test padic_decode(m, bits) == msg
    # the PR-only coder survives only because integer truncation eventually
    # breaks the exact straddle; it must still roundtrip
    @test padic_decode_pr(m, padic_encode_pr(m, msg)) == msg
end

@testset "roundtrip PR+AR coder, P = 3 and P = 5" begin
    rng = MersenneTwister(3)
    for P in (3, 5)
        m = PadicModel([5, 3, 2, 1, 1]; P = P, N = 20)
        for len in (0, 1, 50, 1000)
            msg = rand(rng, 1:5, len)
            bits = padic_encode(m, msg)
            @test all(b -> 0 <= b < P, bits)
            @test padic_decode(m, bits) == msg
        end
    end
end

@testset "byte alphabet, skewed model" begin
    rng = MersenneTwister(4)
    freqs = [rand(rng, 1:100) for _ in 1:256]
    m = PadicModel(freqs; P = 2, N = 40)
    msg = [rand(rng, 1:256) for _ in 1:3000]
    @test padic_decode(m, padic_encode(m, msg)) == msg
end

@testset "compression approaches entropy" begin
    # freqs {1/2, 1/4, 1/8, 1/8} -> entropy 1.75 bits/symbol; frequencies are
    # scaled up so that the EOM weight of 1 does not distort the model
    m = PadicModel([4, 2, 1, 1] .* 1000; P = 2, N = 32)
    rng = MersenneTwister(5)
    len = 20_000
    msg = [(x = rand(rng); x < 0.5 ? 1 : x < 0.75 ? 2 : x < 0.875 ? 3 : 4)
           for _ in 1:len]
    bits = padic_encode(m, msg)
    rate = length(bits) / len
    @test 1.70 < rate < 1.80
end

@testset "optimized P = 2 coder matches reference" begin
    rng = MersenneTwister(7)
    freqs = [4, 2, 1, 1]
    m = PadicModel(freqs; P = 2, N = 32)
    m2 = PadicModel2(freqs; N = 32)
    for len in (0, 1, 10, 1000, 5000)
        msg = rand(rng, 1:4, len)
        bits2 = padic_encode2(m2, msg)
        @test collect(Int, bits2) == padic_encode(m, msg)  # bit-identical
        @test padic_decode2(m2, bits2) == msg
        @test padic_decode(m, collect(Int, bits2)) == msg  # cross-decodable
    end

    # AR stress: interval straddling 1/2
    ms = PadicModel([1, 2, 1]; P = 2, N = 24)
    m2s = PadicModel2([1, 2, 1]; N = 24)
    msg = fill(2, 2000)
    bits2 = padic_encode2(m2s, msg)
    @test collect(Int, bits2) == padic_encode(ms, msg)
    @test padic_decode2(m2s, bits2) == msg

    # byte alphabet
    bfreqs = [rand(rng, 1:100) for _ in 1:256]
    mb = PadicModel(bfreqs; P = 2, N = 40)
    m2b = PadicModel2(bfreqs; N = 40)
    bmsg = [rand(rng, 1:256) for _ in 1:3000]
    bits2 = padic_encode2(m2b, bmsg)
    @test collect(Int, bits2) == padic_encode(mb, bmsg)
    @test padic_decode2(m2b, bits2) == bmsg
end

println("all tests passed")
