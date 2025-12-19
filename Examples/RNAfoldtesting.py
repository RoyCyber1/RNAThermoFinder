import RNA

seq = "AUAUAGCGAACUGCUAUAGAAAUAAUUACACAAUACGGUUUGUUACUGGAAUCAAUCGUGAGCAAGCUUGAGUGAGCCAUUAUG"

# Test different settings
configs = [
    ("Default", {}),
    ("No lonely pairs", {"noLP": 1}),
    ("Dangles=0", {"dangles": 0}),
    ("Dangles=1", {"dangles": 1}),
    ("Dangles=3", {"dangles": 3}),
    ("No GU", {"noGU": 1}),
]

for name, params in configs:
    md = RNA.md()
    for key, val in params.items():
        setattr(md, key, val)

    fc = RNA.fold_compound(seq, md)
    structure, mfe = fc.mfe()
    print(f"{name:20s}: {structure} ({mfe:.2f})")