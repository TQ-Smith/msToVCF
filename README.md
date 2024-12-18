
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
   -l/--length INT            Sets length of segment in number of base pairs. Default 1,000,000.
   -u/--unphased              If set, the phase is removed from genotypes.
   -m/--missing DOUBLE        Genotypes are missing with supplied probability. Default 0.
   -c/--compress              If set, the resulting files are compressed.
```