# DoRothEA v2


DoRothEA (Discriminant Regulon Expression Analysis) is a framework to estimate single sample TF activities from gene expression data and consensus TF-target binding networks. The approach assumes that the activity of a TF can be estimated from the mRNA levels of its direct target genes.  

This new version of DoRothEA (v2) provides **updated scored TF regulons** derived from a broader collection of resources and strategies. The new TF regulons are **signed** (to account for activation/repression), when possible, and accompanied by a **confidence score**. The manuscript is published in [Genome Research](https://genome.cshlp.org/content/early/2019/07/24/gr.240663.118.abstract). 

An earlier version of the TF regulons and the code to estimate TF activities, as described in [Garcia-Alonso et al 2018](http://cancerres.aacrjournals.org/content/early/2017/12/09/0008-5472.CAN-17-1679), can be found in [DoRothEA v1 version](https://github.com/saezlab/DoRothEA/releases/tag/version1).


## About the TF regulons

DoRothEA makes use of consensus TF regulon (CTFR), composed of TF–target regulatory interactions derived from about 20 sources and strategies covering different lines of evidence, including: 13 manually curated repositories, interactions derived from ChIP-seq binding data (from ReMap), in silico prediction of TF binding on gene promoters (using TF binding motifs from JASPAR and HOCMOCO) and the prediction of transcriptional interactions from GTEx (gene expression across human tissues) via ARACNe. 
Each TF-target interaction has been assigned a confidence score, ranging from A-E, being A the most confident interactions (see table below).

| Confidence score  | #Interactions |
| ----------------- | ----- | 
| A                 |  5,869         |
| B                 |  8,991         |
| C                 |  17,519         | 
| D                 |  281,632        |
| E                 |  763,110       |
| Total             |  1,077,121       |

See [Garcia-Alonso 2018 et al.](https://www.biorxiv.org/content/early/2018/06/03/337915) for more information.

## Accessing the TF regulons

We provide the collection of consensus TF regulons in three different ways:
1. Via [omnipathdb](https://www.nature.com/articles/nmeth.4077?proof=trueIn) pypath or webservice.
2. As csv files ``data/TFregulons/consensus/table``.
3. As R objects ready to be used by the VIPER method ([Alvarez et al. 2016](https://www.nature.com/articles/ng.3593)) ``data/TFregulons/consensus/Robjects_VIPERformat``

Please visit our [GitHub page](https://saezlab.github.io/DoRothEA/) for more information. 

### Loading TF regulons in ``pypath``

[``pypath``](https://github.com/saezlab/pypath) is our Python module for building molecular networks.
By loading TF regulons in ``pypath`` you will be able to manipulate it as an ``igraph`` network object,
combine it with annotations from other data sources and also with other networks.
See [here](https://github.com/saezlab/pypath/blob/master/tfregulons_tutorial.md) how to do it.
Briefly, you can build a network of the `A` and `B` confidence level TF-target relationships like this:

```
import pypath

transc = pypath.data_formats.transcription
transc['tfregulons'].inputArgs['levels'] = {'A', 'B'}

pa = pypath.PyPath()
pa.init_network(transc)
```

### Query TF regulons by webservice

TF regulons data is accessible also in the webservice at http://omnipathdb.org/.
Below we show a few example queries, check
[here](https://github.com/saezlab/pypath/blob/master/README.rst) to see some more.

**Important:** Be aware that the server serves not only TF regulons but other datasets as well.
It queries TF regulons only if you explicitely tell to do so by adding ``datasets=tfregulons``
or ``types=TF`` to your query. If you miss to add either of these the returned interactions
will be of other datasets. If you mistype any argument name or value the server returns a
plain text error message pointing out the error.

**Important:** By default the server returns interactions of confidence levels `A` and `B`. If
you want to retrieve other levels you need to explicitely add the argument
``tfregulons_levels=A,B,C,D``. In the webservice only confidence levels `A-D` are available.
The number of interactions in `E` confidence level is too large hence these are available in
static files [here](http://saezlab.org/tfregulons/).

Get all interactions at confidence levels `A-C`:

http://omnipathdb.org/interactions?datasets=tfregulons&tfregulons_levels=A,B,C&genesymbols=1&fields=sources,tfregulons_level

Interactions at confidence level `A` translated to mouse identifiers by homology using NCBI Homologene:

http://omnipathdb.org/interactions?datasets=tfregulons&tfregulons_levels=A&genesymbols=1&fields=sources,ncbi_tax_id&organisms=10090

Interactions from ChIP-Seq and expression based inference methods:

http://omnipathdb.org/interactions?datasets=tfregulons&tfregulons_methods=chipseq,coexp&genesymbols=1&fields=sources,tfregulons_level

All targets of some of the forkhead box transcription factors:

http://omnipathdb.org/interactions?datasets=tfregulons&sources=FOXA1,FOXA2,FOXA3,FOXB1,FOXB2,FOXC1,FOXH1&genesymbols=1&fields=sources,tfregulons_level

All transcription factors regulating EGFR:

http://omnipathdb.org/interactions?datasets=tfregulons&targets=P00533&genesymbols=1&fields=sources,tfregulons_level

Transcriptional regulation of EGFR with its protein-protein interactions and miRNA regulators:

http://omnipathdb.org/interactions?datasets=tfregulons,omnipath,kinaseextra,mirnatarget&partners=EGFR&genesymbols=1&fields=sources,references,type

The same in JSON format:

http://omnipathdb.org/interactions?datasets=tfregulons,omnipath,kinaseextra,mirnatarget&partners=EGFR&genesymbols=1&fields=sources,references,type&format=json

Include 4 additional columns with ``True`` and ``False`` values according to which of the
4 approaches (literature curation, ChIP-Seq, expression based inference, sequence based
binding site prediction) confirmed the TF-target interactions:

http://omnipathdb.org/interactions?datasets=tfregulons&genesymbols=1&fields=databases,tfregulons_curated,tfregulons_chipseq,tfregulons_coexp,tfregulons_tfbs,tfregulons_level


## Usage: estimating TF activities

See ``src/example.r`` for several examples on how to get TF activities based on DoRothEA regulons. We recommend the [VIPER](https://www.bioconductor.org/packages/release/bioc/html/viper.html) R method ([Alvarez et al. 2016](https://www.nature.com/articles/ng.3593)) to compute protein activities. This approach considers the effect signe of the TF-target interaction. Alternatively, you can use any Gene Set Enrichment Analysis method of your preference such as the ones included in the [GSVA package](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html) ([Hänzelmann et al 2013](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7)) or [AUCell](https://bioconductor.org/packages/release/bioc/html/AUCell.html) ([Aibar et al 2017](https://www.nature.com/articles/nmeth.4463)). Please, cite the method accordingly.


## Citation

### DoRothEA v2
[Garcia-Alonso et al 2019](https://genome.cshlp.org/content/early/2019/07/24/gr.240663.118.abstract)
Benchmark and integration of resources for the estimation of human transcription factor activities.
Genome Research, 29(8), 1363-1375. https://doi.org/10.1101/gr.240663.118

```
@article{garcia2019benchmark,
  doi = {10.1101/gr.240663.118},
  url = {https://genome.cshlp.org/content/early/2019/07/24/gr.240663.118.abstract},
  year  = {2019},
  month = {aug},
  volume = {29},
  number = {8},
  pages = {1363-1375},
  title={Benchmark and integration of resources for the estimation of human transcription factor activities},
  author={Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D. and Saez-Rodriguez, J.}
  journal={Genome Research}
```

### DoRothEA v1
[Garcia-Alonso et al 2018](https://www.ncbi.nlm.nih.gov/pubmed/29229604)
Transcription Factor Activities Enhance Markers of Drug Sensitivity in Cancer.
Cancer Res February 1 2018 (78) (3) 769-780; 
DOI: 10.1158/0008-5472.CAN-17-1679


```
@article{garcia2018transcription,
  doi = {10.1158/0008-5472.can-17-1679},
  url = {https://doi.org/10.1158/0008-5472.can-17-1679},
  year  = {2018},
  month = {feb},
  publisher = {American Association for Cancer Research ({AACR})},
  volume = {78},
  number = {3},
  pages = {769--780},
  author = {Luz Garcia-Alonso and Francesco Iorio and Angela Matchan and Nuno Fonseca and Patricia Jaaks and Gareth Peat and Miguel Pignatelli and Fiammetta Falcone and Cyril H. Benes and Ian Dunham and Graham Bignell and Simon S. McDade and Mathew J. Garnett and Julio Saez-Rodriguez},
  title = {Transcription Factor Activities Enhance Markers of Drug Sensitivity in Cancer},
  journal = {Cancer Research}
}
```


## License

Distributed under the GNU GPLv2 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/DoRothEA/blob/master/LICENSE.txt) or copy at https://www.gnu.org/licenses/gpl-2.0.html.
