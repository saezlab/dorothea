# DoRothEA v2


DoRothEA (Discriminant Regulon Expression Analysis) is a framework to estimate single sample TF activities from gene expression data and consensus TF-target DNA binding networks. The approach assumes that the activity of a TF can be estimated from the mRNA levels of its direct target genes.  


This new version of DoRothEA (v2) provides updated TF regulons derived from a broader collection of resources and strategies. The new TF regulons are signed (to account for activation/repression), when possible, and accompanied by a confidence score. The preprint of the corresponding paper can be found on [bioRxiv](https://www.biorxiv.org/content/early/2018/06/03/337915) 


An earlier version of the TF regulons and the code to estimate TF activities, as described in [Garcia-Alonso et al 2018](http://cancerres.aacrjournals.org/content/early/2017/12/09/0008-5472.CAN-17-1679), can be found in [DoRothEA v1 version](https://github.com/saezlab/DoRothEA/releases/tag/version1).



## Usage

See ``src/example.r`` for several examples.



## About the TF regulons

DoRothEA makes use of consensus TF regulon (CTFR), composed of TF–target regulatory interactions derived from about 20 sources and strategies covering different lines of evidence, including: 13 manually curated repositories, interactions derived from ChIP-seq binding data (from ReMap), in silico prediction of TF binding on gene promoters (using TF binding motifs from JASPAR and HOCMOCO) and the prediction of transcriptional interactions from GTEx (gene expression across human tissues) via ARACNe. 
Each TF-target interaction has been assigned a confidence score, ranging from A-E, being A the most confident interactions. See (Garcia-Alonso 2018 et al.)[unknown] for mor information.


The collection of consensus TF regulons, scored according our A-E criteria, is available at  ``data/TFregulons/Robjects_VIPERformat/consensus/BEST_viperRegulon.rdata`` and ``data/TFregulons/table/database_20180326.csv.zip``


## Citation

[Garcia-Alonso et al 2018](https://www.ncbi.nlm.nih.gov/pubmed/29229604)
Transcription Factor Activities Enhance Markers of Drug Sensitivity in Cancer.
Cancer Res February 1 2018 (78) (3) 769-780; 
DOI: 10.1158/0008-5472.CAN-17-1679


```
@article{garcia2018transcription,
  title={Transcription factor activities enhance markers of drug sensitivity in cancer},
  author={Garcia-Alonso, Luz and Iorio, Francesco and Matchan, Angela and Fonseca, Nuno and Jaaks, Patricia and Peat, Gareth and Pignatelli, Miguel and Falcone, Fiammetta and Benes, Cyril H and Dunham, Ian and others},
  journal={Cancer research},
  volume={78},
  number={3},
  pages={769--780},
  year={2018},
  publisher={American Association for Cancer Research}
}
```


## License

Distributed under the GNU GPLv2 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/DoRothEA/blob/master/LICENSE.txt) or copy at https://www.gnu.org/licenses/gpl-2.0.html.

