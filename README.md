
# msToVCF v1.0

Converts ms style output to VCF format.

## Building msToVCF

The resulting executable will be in the **msToVCF/bin** directory.

```
git clone https://github.com/TQ-Smith/msToVCF.git 
cd msToVCF
make
```

## Options

```
Usage: msToVCF [options] <inFile.ms.gz>
Options:
   --help                  Prints help menu and exits.
   --version               Prints version number and exits.
   --length INT            Sets length of segment in number of base pairs. Default 1,000,000.
   --unphased              If set, the phase is removed from genotypes.
   --missing DOUBLE        Genotypes are missing with supplied probability. Default 0.
   --compress              If set, the resulting files are compressed.
```