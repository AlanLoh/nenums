# **nenums**

[![nenums](hhttps://img.shields.io/pypi/v/nenums.svg)](
    https://pypi.python.org/pypi/nenums)

**Alan Loh** (LESIA, Paris Observatory), **Julien Girard** (AIM CEA-Saclay, Paris Diderot University) and the *NenuFAR* team

`nenums` is a Python package designed to convert [*NenuFAR*](https://nenufar.obs-nancay.fr) *Cross-correlation Statistics* data (**XST**) to Measurement Sets.

## Installation
### Easy install
Installing *nenums* can be done via [`pip`](https://pypi.org/project/pip/):
```
pip install nenums
```
The package is updated with the addition of `--upgrade`.

### Requirements
* [*astropy*](http://www.astropy.org)
* [*pyrap*](https://github.com/casacore/python-casacore)
* [*makems*](https://github.com/ska-sa/makems)

### See also
The set of radio astronomical software packages [`KERN`](http://kernsuite.info).

## Usage
### Command-line
`nenums` can be called directly via the shell to quickly convert a XST observation into a MS:
```
nenums -xst /path/to/DATA_XST.fits -msname /path/to/new/directory/example.ms
```

### Python Shell
```
from nenums import MS
ms = MS(xst='/path/to/DATA_XST.fits', msname='/path/to/new/directory/example.ms')
ms.createMS()
```
Will perform the same operation as the above commande-line example.

However, one could also access some Measurement-Set specific functions:
```
from nenums.utils import ms
ms.splitMS(msname='/path/to/new/directory/example.ms', remove=True)
```
