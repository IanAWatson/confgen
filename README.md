# Confgen

This is a minor adaptation of a nice conformer generation script published by iwatobipen (name looks familiar, but not me).

[wordpress](https://iwatobipen.wordpress.com/2021/01/31/generate-conformers-script-with-rdkit-rdkit-chemoinformatics/)

Uses RDKit EmbedMultipleConfs to generate conformers.

This version reads a file of smiles and generates conformers for each molecule read.

External dependencies are
+ [RDKit](http://rdkit.org/)
# Usage

```
python confgen.py --help

Usage: confgen.py [OPTIONS]

  Performs 2D-3D conversion and conformer generation on smiles input

Options:
  -i, --input TEXT            inputfile smiles file  [required]
  -o, --output TEXT           output sdf file path
  -t, --prunermsthresh FLOAT  RMS threshold for conformer retention
  --numconf INTEGER           Number of conformations to produce (def 50)
  -r, --add_ref BOOLEAN       Also write the starting molecule
  --rmchiral / --normchiral   Remove chirality before processing
  -v, --verbose               Verbosity
  --help                      Show this message and exit.
```

New here is the --rmchiral option, which removes chirality from all input molecules
before processing. This is useful if the input smiles contains invalid chirality,
which then induces failures in AllChem.MMFFOptimizeMoleculeConfs. Unfortunately
this is applied to all input molecules. Better would be to detect the failures,
remove chirality from those, and retry them without chirality.
