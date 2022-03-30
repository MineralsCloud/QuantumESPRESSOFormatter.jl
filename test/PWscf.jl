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

# Data from https://gitlab-hpc.cineca.it/pbonfa01/q-e/-/blob/c9f92947/PW/examples/example11/run_example#L86-117
@testset "Test `format_text`" begin
    str = raw"""
     &control
        calculation = 'scf'
    prefix='Fe',
        pseudo_dir = 'pseudo/',
        outdir='tmp/'
     /
     &system
        ibrav=  3,
        celldm(1) =5.42,
        nat=  1,
        ntyp= 1,
        nr1=27,
        nr2=27,
        nr3=27,
        noncolin=.true.
        lspinorb=.true.
        starting_magnetization(1)=0.5,
        occupations='smearing',
        smearing='mp',
        degauss=0.04,
        ecutwfc =45.0, ecutrho =300.0
     /
     &electrons
        conv_thr =  1.0d-10
     /
    ATOMIC_SPECIES
    Fe  0.0    Fe.rel-pbe-kjpaw.UPF
    ATOMIC_POSITIONS   alat
    Fe  0.0000000   0.00000000   0.0
    K_POINTS AUTOMATIC
    8 8 8 1 1 1
    """
    @test format_text(str) == "&CONTROL\n    prefix = 'Fe'\n    outdir = 'tmp/'\n    pseudo_dir = 'pseudo/'\n/\n&SYSTEM\n    starting_magnetization(1) = 0.5\n    ecutwfc = 45.0\n    celldm(1) = 5.42\n    ecutfock = 300.0\n    occupations = 'smearing'\n    ntyp = 1\n    smearing = 'mp'\n    ibrav = 3\n    degauss = 0.04\n    nr2 = 27\n    nr3 = 27\n    lspinorb = .true.\n    ecutrho = 300.0\n    noncolin = .true.\n    nr1 = 27\n    nat = 1\n/\n&ELECTRONS\n    conv_thr = 1.0e-10\n/\n&IONS\n/\n&CELL\n/\nATOMIC_SPECIES\n      Fe    0.000000000 Fe.rel-pbe-kjpaw.UPF\nATOMIC_POSITIONS { alat }\n      Fe    0.000000000    0.000000000    0.000000000\nK_POINTS { automatic }\n        8     8     8     1     1     1\n\n\n\n\n\n"
end

end
