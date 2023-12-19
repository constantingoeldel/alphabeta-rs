CREATE TABLE default.methylome (
    id UInt64 DEFAULT rand64(),
    chromosome UInt8,
    location UInt32,
    strand Enum8('+' = 1, '-' = 2),
    context Enum8('CG' = 1, 'CHG' = 2, 'CHH' = 3),
    count_methylated UInt8,
    count_total UInt8,
    posteriormax Float32,
    status Enum8('U' = 1, 'M' = 2, 'I' = 3),
    meth_lvl Float32,
    trinucelotide_context String,
    pedigree String,
    generation UInt8,
    line UInt8
) ENGINE = MergeTree SAMPLE BY id PARTITION BY pedigree
ORDER BY (
        pedigree,
        generation,
        chromosome,
        strand,
        location,
        id
    )