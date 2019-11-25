![ELIXIR_GREECE_and_ESPA_logos](http://genomics-lab.fleming.gr/fleming/reczko/elixir/logos/ELIXIR_GREECE_and_ESPA_logos-1338x218.png)


# haplotypeAnalysisPipeline

This pipeline prepares the population genotypes ped/map files, runs them on SHAPEIT (phasing), optionally uses the results as input for IMPUTE2 (imputation) using the 1000G phase III reference panel, and finally converts them to Chromopainter_v2 formats in order to perform with fineSTRUCTURE fs 2.0.7 two clusterings of donors and recipients respectively, and with Chromopainter_v2 a painting of the resulting recipient clusters out of the donor clusters. The results can be then used as input for Globetrotter to find admixture datings, and for fineSTRUCTURE to find ancestry components of the donor clusters in the recipient clusters.
