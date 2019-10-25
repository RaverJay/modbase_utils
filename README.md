# modbase_utils
script(s) to handle guppy basecalls with modified bases


## get_modbases_from_fast5.py
Gives per read positions of '6mA 5mC' modifications with at least 0.5 caller confidence

Usage:
```python3 get_modbases_from_fast5.py [fast5_basecalls.fast5]...```

Output on stdout:
```
>[readid] 6mA 5mC [metadata]...
6mA pos1,pos2,...
5mC pos1,pos2,...
```
