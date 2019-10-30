# modbase_utils
script(s) to handle guppy basecalls with modified bases


## get_modbases_from_fast5.py
Gives per read positions of '6mA 5mC' modifications with at least 0.5 caller confidence

Usage:
```
get_modbases_from_fast5.py [-h] [-t FLOAT] [-m STRING] DIR

Extract modified bases from guppy fast5 basecalls. This only applies if
--modbase was enabled for basecalling.

positional arguments:
  DIR                   folder containing fast5 files with modbase basecalls

optional arguments:
  -h, --help            show this help message and exit
  -t FLOAT, --threshold_min_mod FLOAT
                        minimum caller confidence of modified bases that will
                        be extracted (default: 0.5)
  -m STRING, --modbases_string STRING
                        modified base(s) to extract need to be comma-separated
                        strings with structure: name,symbol,canonicalbase
                        (default: ['6mA,Y,A', '5mC,Z,C'])
```

Output on stdout:
```
>[readid] [modbase,...] [metadata]...
6mA pos1,pos2,...
5mC pos1,pos2,...
```
