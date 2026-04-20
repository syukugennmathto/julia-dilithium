module Fips202

using SHA

export shake128, shake256, shake128!, shake256!

function shake128(output_len::Int, input::Vector{UInt8})
    return SHA.shake128(input, UInt64(output_len))
end

function shake256(output_len::Int, input::Vector{UInt8})
    return SHA.shake256(input, UInt64(output_len))
end

function shake256!(output::Vector{UInt8}, input::Vector{UInt8})
    result = shake256(length(output), input)
    copyto!(output, result)
    return nothing
end

function shake128!(output::Vector{UInt8}, input::Vector{UInt8})
    result = shake128(length(output), input)
    copyto!(output, result)
    return nothing
end

end
