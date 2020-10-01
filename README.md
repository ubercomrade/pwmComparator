# pwmComparator is pipeline for comparing pwm models on the same data


## Requirements

PYTHON:
  * cython: `pip3 install cython`

TOOLS:
  * bedtools: https://bedtools.readthedocs.io/en/latest/  version >= 2.26.0

## Installation

```  
git clone https://github.com/ubercomrade/pipeline.git  
cd pipeline/  
pip3 install -e .  
```

## Usage
The command `pmwComparator.py -h` return:

```
usage: pwmComparator.py [-h] [-nnames N [N ...]] [-f FPR] [-C CPU_COUNT]
                        bed N genome output N [N ...]

positional arguments:
  bed                   path to BED file
  N                     promoters of organism (hg38, mm10)
  genome                path to genome fasta file
  output                output dir
  N                     list of path to PFMs to use (./model1.pfm,
                        ./model2.pfm ...)

optional arguments:
  -h, --help            show this help message and exit
  -names N [N ...]     list of names for models (name1, name2 ...)
  -f FPR, --FPR FPR     FPR, def=1.9*10^(-4)
```

Example run:

```
pwmComparator.py \
path/to/peaks.bed \
hg38 \
/path/to/genomes/hg38.fa \
./output_dir \
./model1.pfm \
./model2.pfm \
./model3.pfm \
```

## License

Copyright (c) 2020 Anton Tsukanov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.