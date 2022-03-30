module PWscf

using StructArrays: StructArray
using Test
using QuantumESPRESSOBase.Inputs.PWscf: AtomicSpecies, AtomicPosition, AtomicSpeciesCard, AtomicPositionsCard
using QuantumESPRESSOFormatter.Inputs.PWscf: format_text

@testset "Test constructing `AtomicSpeciesCard`" begin
    a = ["Al", "As"]
    m = [24590.7655930491, 68285.4024548272]
    pp = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
    card = AtomicSpeciesCard(StructArray{AtomicSpecies}((a, m, pp)))
    @test string(card) ==
          "ATOMIC_SPECIES\n      Al 24590.765593049 Al.pbe-n-kjpaw_psl.1.0.0.UPF\n      As 68285.402454827 As.pbe-n-kjpaw_psl.1.0.0.UPF"
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
    @test string(card) ==
          "ATOMIC_POSITIONS { alat }\n       S    0.500000000    0.288675130    1.974192764\n      Mo    0.000000000    0.577350270    2.462038339\n       S    0.000000000   -0.577350270    2.950837559"
end

@testset "Test `format_text`" begin
    str = raw"""
    &control
    calculation='scf'
    pseudo_dir  = '/home/qe/pseudo'
    prefix='GaN'
    outdir      = './tmp'
    /
    &system
    ibrav = 4
    celldm(1)=5.95484286816
    celldm(3)=1.63011343669
    nat = 4
    ntyp = 2
    ecutwfc = 160
    /
    &electrons
    conv_thr=1.0d-10
    /
    ATOMIC_SPECIES
    Ga  69.723   Ga.pbe-dn-kjpaw_psl.1.0.0.UPF
    N   14.007   N.pbe-n-kjpaw_psl.1.0.0.UPF
    ATOMIC_POSITIONS (crystal)
    Ga       0.666666667   0.333333333  -0.000051966
    N        0.666666667   0.333333333   0.376481188
    Ga       0.333333333   0.666666667   0.499948034
    N        0.333333333   0.666666667   0.876481188
    K_POINTS automatic
    6 6 6 1 1 1
    """
    @test format_text(str) == "&CONTROL\n    prefix = 'GaN'\n    outdir = './tmp'\n    pseudo_dir = '/home/qe/pseudo'\n/\n&SYSTEM\n    nat = 4\n    ecutwfc = 160.0\n    celldm(1) = 5.95484286816\n    celldm(3) = 1.63011343669\n    ecutrho = 640.0\n    ecutfock = 640.0\n    ntyp = 2\n    ibrav = 4\n/\n&ELECTRONS\n    conv_thr = 1.0e-10\n/\n&IONS\n/\n&CELL\n/\nATOMIC_SPECIES\n      Ga   69.723000000 Ga.pbe-dn-kjpaw_psl.1.0.0.UPF\n       N   14.007000000 N.pbe-n-kjpaw_psl.1.0.0.UPF\nATOMIC_POSITIONS { crystal }\n      Ga    0.666666667    0.333333333   -0.000051966\n       N    0.666666667    0.333333333    0.376481188\n      Ga    0.333333333    0.666666667    0.499948034\n       N    0.333333333    0.666666667    0.876481188\nK_POINTS { automatic }\n        6     6     6     1     1     1\n\n\n\n\n\n"
end

end
