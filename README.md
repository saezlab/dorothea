# DoRothEA v2


DoRothEA (Discriminant Regulon Expression Analysis) is a framework to estimate single sample TF activities from gene expression data and consensus TF-target DNA binding networks. The approach assumes that the activity of a TF can be estimated from the mRNA levels of its direct target genes.  


This new version of DoRothEA (v2) provides updated TF regulons derived from a broader collection of resources and strategies. The new TF regulons are signed (to account for activation/repression), when possible, and accompanied by a confidence score. The preprint of the corresponding paper can be found on [bioRxiv](https://www.biorxiv.org/content/early/2018/06/03/337915). 


An earlier version of the TF regulons and the code to estimate TF activities, as described in [Garcia-Alonso et al 2018](http://cancerres.aacrjournals.org/content/early/2017/12/09/0008-5472.CAN-17-1679), can be found in [DoRothEA v1 version](https://github.com/saezlab/DoRothEA/releases/tag/version1).



## Usage

See ``src/example.r`` for several examples.



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


The collection of consensus TF regulons, scored according our A-E criteria, is available at  ``data/TFregulons/Robjects_VIPERformat/normal/``

Please visit our [GitHub page](https://saezlab.github.io/dorothea/) for more information. 

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

## Citation

### DoRothEA v2
[Garcia-Alonso et al 2018](https://www.biorxiv.org/content/early/2018/06/03/337915)
Benchmark and integration of resources for the estimation of human transcription factor activities.
BioRxiv;
DOI: 10.1101/337915

```
@article{garcia2018benchmark,
  doi = {10.1101/337915},
  url = {https://www.biorxiv.org/content/early/2018/06/03/337915},
  year  = {2018},
  month = {jun},
  publisher = {},
  volume = {},
  number = {},
  pages = {},
  title={Benchmark and integration of resources for the estimation of human transcription factor activities},
  author={Garcia-Alonso, Luz and Ibrahim, MM and Turei, D and Saez-Rodriguez, J}
  journal={bioRxiv}
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
