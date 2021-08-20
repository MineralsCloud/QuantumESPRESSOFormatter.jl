module PWscf

using AbInitioSoftwareBase.Inputs: asstring
using StructArrays: StructArray
using QuantumESPRESSOBase.Inputs.PWscf

@testset "Test constructing `AtomicSpeciesCard`" begin
    a = ["Al", "As"]
    m = [24590.7655930491, 68285.4024548272]
    pp = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
    card = AtomicSpeciesCard(StructArray{AtomicSpecies}((a, m, pp)))
    @test asstring(card) ==
          "ATOMIC_SPECIES\n      Al 24590.765593049 Al.pbe-n-kjpaw_psl.1.0.0.UPF\n      As 68285.402454827 As.pbe-n-kjpaw_psl.1.0.0.UPF\n      Si 25591.192491355 Si.pbe-n-kjpaw_psl.1.0.0.UPF"
end

@testset "Test constructing `AtomicPositionsCard`" begin
    # Data from https://github.com/QEF/q-e/blob/7be27df/PW/examples/gatefield/run_example#L129-L132.
    a = ["S", "Mo", "S"]
    pos = [
        [0.500000000, 0.288675130, 1.974192764],
        [0.000000000, 0.577350270, 2.462038339],
        [0.000000000, -0.577350270, 2.950837559],
    ]
    card = AtomicPositionsCard(
        StructArray{AtomicPosition}((a, pos, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])),
        "alat",
    )
    @test asstring(card) ==
          "ATOMIC_POSITIONS { alat }\n       S    0.500000000    0.288675130    1.974192764\n      Mo    0.000000000    0.577350270    2.462038339\n       S    0.000000000   -0.577350270    2.950837559"
end

end
