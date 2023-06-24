nextflow.enable.dsl = 2

params.path ="/home/cq/visualStudio/fastqc/data" 
params.with_fastqc=false
params.with_summrize = false
files_channel = Channel.fromPath("${params.path}", checkIfExists: true )
params.classifier="${params.path}/gg-13-8-99-515-806-nb-classifier.qza"
params.seqs="${params.path}/qiime-result/rep-seqs.qza"
params.denoisingstats= "${params.path}/qiime-result/denoising-stats.qza"

params.table="${params.path}/qiime-result/table.qza"
params.manifest="manifest.tsv"  // your manifest file name
params.metadata="metadata.tsv" // your metdata file name
metadata_channel=Channel.fromPath("${params.path}/${params.metadata}", checkIfExists: true )

params.trim_left="19"
params.trunc_len="100"



process importing_ch {
 publishDir "${params.path}/qiime-result", mode: 'copy'
  input:
	 file files

  output:
    file "single-end-demux.qza" 

	    
	 """
  qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ${params.path}/${params.manifest}\
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2
     
	 """
	}

process summrize_data {
 publishDir "${params.path}/qiime-result", mode: 'copy'
  input:
	 file "singel-end-demux.qza"

  output:
    file "qualities.qzv" 

	    
	 """
    qiime demux summarize\
    --p-n 100000 \
    --i-data single-end-demux.qza\
    --o-visualization qualities.qzv
	 """
	}

process denoising_step {
 publishDir "${params.path}/qiime-result", mode: 'copy'
  input:
	 file "single-end-demux.qza"

  output:
    file "*.qza" 

	    
	 """
    qiime dada2 denoise-single\
    --i-demultiplexed-seqs single-end-demux.qza\
    --p-trim-left ${params.trim_left}\
    --p-trunc-len ${params.trunc_len}\
    --p-n-threads 2\
    --o-table table.qza\
    --o-representative-sequences rep-seqs.qza\
    --o-denoising-stats denoising-stats.qza

	 """
     
	}

process visulizing_step {
 publishDir "${params.path}/qiime-result", mode: 'copy'
  input:
   file   "*.qza"
   path   seqs
	 path  denoisingstats
   path   table
  

  output:
    file "*.qzv" 

	    
	 """
    qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization table.qzv\

   qiime feature-table tabulate-seqs \
    --i-data rep-seqs.qza \
    --o-visualization rep-seqs.qzv

  qiime metadata tabulate \
   --m-input-file denoising-stats.qza\
   --o-visualization denoising-stats.qzv

	 """
	}

// Taxonomy
//We will use a Bayes classifier trained on the GreenGenes database which can be downloaded from https://docs.qiime2.org/2021.4/data-resources/.
// or use this command wget https://data.qiime2.org/2021.8/common/gg-13-8-99-515-806-nb-classifier.qza
 process taxonomy_step {
 publishDir "${params.path}/qiime-result", mode: 'copy'
  input:
   file "*.qza"
	 path   seqs
   path  classifier

  output:
    file "taxa.qza" 

	    
	 """
    qiime feature-classifier classify-sklearn \
    --i-reads rep-seqs.qza \
    --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
    --p-n-jobs 2 \
    --o-classification taxa.qza
	 """
	}
 
process look_taxonomy {
 publishDir "${params.path}/qiime-result", mode: 'copy'
  input:
	 path table 
   file taxa
   file metadata

  output:
    file "taxa_barplot.qzv"  
	 """
    qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy taxa.qza \
    --m-metadata-file metadata.tsv\
    --o-visualization taxa_barplot.qzv
	 """
	}


workflow {  
if(params.path== null) {
        print "ERROR: Please provide an path to your *.fastq.gz folder(e.g. --path  home/fastQfiles)"
        System.exit(1)
    }
 

 import_channel=importing_ch(files_channel)  
 importing_channel=import_channel
 if(params.with_summrize != false) {
      summrize_data(import_channel) 
    }
 
 denosing_channel=denoising_step(importing_channel)
 dada2_ch=denosing_channel
 visulaizing_channel=visulizing_step(denosing_channel,Channel.from(params.seqs),Channel.from(params.table),Channel.from(params.denoisingstats))
 taxonomy_ch=taxonomy_step(dada2_ch,Channel.from(params.seqs),Channel.from(params.classifier))
 taxonomy_ch=look_taxonomy(Channel.from(params.table),metadata_channel,taxonomy_ch)
 
}


//Wedad Ahmed