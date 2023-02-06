# sangerseq_viewer Installation and User Manual
sangerseq_viewer is a Python package for visualizing Sanger sequencing results and the corresponding annotated sequence maps automatically.  
Despite being an essential task in DNA sequence construction and editing, there is a lack of open-source software that provides a user-friendly graphical representation of Sanger sequencing results. While commercial GUI software like Snapgene and Geneious Prime offer this functionality, manually processing large amounts of data from Sanger sequencing results can be a tedious and time-consuming task.  
Sangerseq_viewer provides a solution by allowing you to visualize sequencing results with just a simple command.

Note: Sangerseq_viewer relies on the packages [patchworklib](https://github.com/ponnhide/patchworklib) and [QUEEN](https://github.com/yachielab/QUEEN), which provide APIs for managing matplotlib subplots and GenBank files respectively. For more information, please refer to their respective documentation.

You can test Sangerseq_viewer through [Google colaboratory](https://colab.research.google.com/drive/1HysWJv6Afi-LFbowWL8_zgEAZcoMpP0X?usp=sharing).

## Software dependency
- python 3.8.0 or later
	

## Installation
Plese execute the following commands.

```
pip install patchworklib
pip install python-queen
pip install sangerseq-viewer 
```

If you cannnot use `sangerseq_viewr` command after the installation, please add `--prefix=PREFIX` option to `pip install sangerseq-viewer` and try re-installation. `PREFIX` is the executable path of sangerseq_viewer.

## Example code
**Example command 1**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/Spec-2xU6gRNA-1.ab1 -o output/example1.png --dpi 200`

**Output figure 1** 
<img src="output/example1.png" width="1200x1200">

---

**Example command 2**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/Spec-2xU6gRNA-1.ab1 -o output/example2.png -l 200 --dpi 200`

**Output figure 2**

<img src="output/example2.png" width="1200x1200">

---

**Example command 3**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/Spec-2xU6gRNA-1.ab1 -o output/example3.png-l 200 -rs 1700 -re 2100 --dpi 200`

**Output figure 3**

<img src="output/example3.png" width="1200x1200">

---

**Example command 4**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/ -o output/example4.png --dpi 200`

**Output figure 4**

<img src="output/example4.png" width="1200x1200">

---

**Example command 5**
`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/ -o output/example5.png -l 200 --dpi 200`

**Output figure 5**

<img src="output/example5.png" width="1200x1200">

---

**Example command 6**

`sangerseq_viewer -s example_data/puc19_spec_2xu6grna.gb -q example_data/ab1/ -o output/example6.png -l 200 -rs 1700 -re 2100 --dpi 200`

**Output figure 6**

<img src="output/example6.png" width="1200x1200">

## Usage

```
SYNOPSIS
sangerseq_viewer [-h] [-q QUERY] [-s SUBJECT] [-l LINEBREAK] [-o OUTPUT] [-rs START] [-re END] [-wq {True,False}] [-d DPI]

  -h, --help            
     Show this help message and exit
  -q, --query, str
     Ab1 file path, path to the directory containing ab1 files, or txt file path describing a ab1 file path on each line.
  -s, --subject, str
     Genbank file path.
  -l, --linebreak, int (optional, default: No line break)
     Sequence length for line break.
  -o, --output, str (optional, default: No output) 
     Output file path. The output image format can be specified by the filename extension.
  -rs, --start, int (optional, default: 0) 
     Start position of the subject sequence region to be visualized. 
  -re, --end, int (optional, default: End position of the alinged region.)
     End position of the subject sequence region to be visualized.
  -wq, --quality, {True,False} (optional, default: True)
     If True, display bar plot representing Quality value at each nucleotide position.
  -d, --dpi, float (optional, default: 200)
     Resolution of the output image. If output format is pdf, the value is ignored.
  -c, --output_cromatogram (optional, defualt: None) 
     Output table file path. If the option is given, the values of Sanger sequencing cromatogram will be output as a csv file.
  -f, --output_fasta (optional, default: None)
     Output FASTA file path. If the option is given, aligned sequences will be output as a Fasta file.
```



If you want to use `sanger_seqviewer` as python module, please import `sangerseq_viewer.sangerseq_viewer` and use `view_sanger()` fucntion. It takes same parameters with `sangerseq_viewer` command and returns `matplotlib.figure` object. 

`sangerseq_viewer.sangerseq_viewer` provides other useful functions such as `generate_consensusseq()` and `ab1_to_dict()` for  handling ab1 file. I will add the document for them later.



