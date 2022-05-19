

#Install version 2.4.1 picrust2
wget https://github.com/picrust/picrust2/archive/v2.4.1.tar.gz
tar xvzf  v2.4.1.tar.gz
cd picrust2-2.4.1/
Create and activate the environment (with requirements) and then install PICRUSt2 with pip.

conda env create -f  picrust2-env.yaml
conda activate picrust2
pip install --editable .



#convert the final rep_seqs.qza to fasta format for picrust

qiime tools export \
  --input-path rep_seqs_final.qza \
  --output-path rep_seqs.fasta

###run the pipeline

picrust2_pipeline.py -s dna-sequences.fasta -i twins_otus.tsv -o picrust2_out_pipeline -p 8



###add functional descriptions

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz