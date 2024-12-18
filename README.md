
# msToVCF v1.0

## Building msToVCF

The resulting executable will be in the **msToVCF/bin** directory.

```
git clone https://github.com/TQ-Smith/msToVCF.git 
cd msToVCF
make
```

## Options

msToVCF reads from **stdin**.

```
Usage: msToVCF [options]
Options:
   -h/--help                  Prints help menu and exits.
   -v/--version               Prints version number and exits.
   -p/--ploidy INT            Sets ploidy of samples. Default 2.
   -l/--length INT            Sets length of segment in number of base pairs. Default 1,000,000.
   -u/--unphased              If set, the phase is removed from genotypes.
   -m/--missing DOUBLE        Genotypes are missing with probability of 0 <= DOUBLE < 1. Default 0.
   -c/--compress              If set, the resulting files are compressed.
```