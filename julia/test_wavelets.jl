using Pkg
Pkg.activate(".")
using Wavelets
using Wavelets: WT

println("Available wavelets in WT module:")
for name in names(WT, all=true)
    if !startswith(string(name), "#")
        println("  $name")
    end
end

# Test specific wavelets
println("\nTesting specific wavelets:")
test_wavelets = ["haar", "db2", "db4", "coif2", "sym4"]
for wf in test_wavelets
    try
        sym = Symbol(lowercase(wf))
        if hasproperty(WT, sym)
            println("  $wf: ✓ Available")
        else
            println("  $wf: ✗ Not available")
        end
    catch e
        println("  $wf: ✗ Error: $e")
    end
end
