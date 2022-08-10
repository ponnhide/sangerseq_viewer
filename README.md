# sangerseq_viewer Installation and User Manual



sangerseq_viewer is a python package to automatically visualize Sanger sequencing results and the corresponding annotated sequence map. 

Despite the fact that Sanger sequencing is one of the indispensable tasks in constructing and editing DNA sequences, there is no open-source software that provides an intuitive and understandable graphical representation of sequencing results.   Some commercial GUI software (Snapgene, Geneious prime) fulfills this requirement. However, now that it is common to design and construct dozens or hundreds of DNA sequences, manually processing large numbers of Sanger sequencing results is a tedious and time-consuming task.

Here, I provide sangerseq_viewer that allow us to visulize sanger sequncing result with a simple command.   

*sangerseq_viewer is depend on m packages [patchworklib](https://github.com/ponnhide/patchworklib) and [QUEEN](https://github.com/yachielab/QUEEN). Both packaeges provide APIs for handling matplotlib subplots and GenBank files, respectivily. If you are interested in them, please see their documents. *

## Software dependency
- python 3.8.0 or later
	

## Installation
Plese execute the following commands.

```
pip install patchworklib
pip install python-queen
pip install sangerseq-viewer --prefix=PREFIX
```
`PREFIX` is the executable path of sangerseq_viewer. 
By specifying `--prefix`, you can run `sangerseq_viewer`  from your command line.

## Example code
**Example command 1**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/Spec-2xU6gRNA-1.ab1 -o output/example1.png --dpi 300`

**Output figure 1** 
<img src="output/example1.png" width="1200x1200">

---

**Example command 2**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/Spec-2xU6gRNA-1.ab1 -o output/example2.png -l 200 --dpi 300`

**Output figure 2**

<img src="output/example2.png" width="1200x1200">

---

**Example command 3**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/Spec-2xU6gRNA-1.ab1 -o output/example3.png-l 200 -rs 1700 -re 2100 --dpi 300`

**Output figure 3**

<img src="output/example3.png" width="1200x1200">

---

**Example command 4**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/ -o output/example4.png --dpi 300`

**Output figure 4**

<img src="output/example4.png" width="1200x1200">

---

**Example command 5**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/ -o output/example5.png -l 200 --dpi 300`

**Output figure 5**

<img src="output/example5.png" width="1200x1200">

---

**Example command 6**

`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/ -o output/example6.png -l 200 -rs 1700 -re 2100 --dpi 300`

**Output figure 6**

<img src="output/example6.png" width="1200x1200">

## Usage

```
SYNOPSIS
sangerseq_viewer [-h] [-q QUERY] [-s SUBJECT] [-l LINEBREAK] [-o OUTPUT] [-rs START] [-re END] [-wq {True,False}] [-d DPI]

  -h, --help            
     show this help message and exit
  -q QUERY, --query QUERY
     Abi file path.
  -s SUBJECT, --subject SUBJECT
     Genbank file path.
  -l LINEBREAK, --linebreak LINEBREAK
     Sequence length for line break.
  -o OUTPUT, --output OUTPUT
     Output file path. The output image format can be specified by the filename extension.
  -rs START, --start START
     Start position of the subject sequence region to be visualized. 
  -re END, --end END    
     End position of the subject sequence region to be visualized.
  -wq {True,False}, --quality {True,False}
     If True, display bar plot representing Quality value at each nucleotide position.
  -d DPI, --dpi DPI     
     Resolution of the output image. If output format is pdf, the value is ignored.
```



If you want to use `sanger_seqviewer` as python module, please import `sangerseq_viewer.sangerseq_viewer` and use `view_sanger()` fucntion. It takes same parameters with `sangerseq_viewer` command and returns `matplotlib.figure` object. 

`sangerseq_viewer.sangerseq_viewer` provides other useful functions such as `generate_consensusseq()` and `ab1_to_dict()` for  handling ab1 file. I will add the document for them later.



