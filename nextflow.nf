$HOSTNAME = ""
params.outdir = 'results'  

def pathChecker(input, path, type){
	def recursive = (type == "folder") ? "--recursive" : ""
	def cmd = "mkdir -p check && mv ${input} check/. "
	if (!input || input.empty()){
		input = file(path).getName().toString()
		cmd = "mkdir -p check && cd check && ln -s ${path} ${input} && cd .."
		if (path.indexOf('s3:') > -1 || path.indexOf('S3:') >-1){
			cmd = "mkdir -p check && cd check && aws s3 cp ${recursive} ${path} ${workDir}/${input} && ln -s ${workDir}/${input} . && cd .."
		} else if (path.indexOf('http') > -1){
			slashCount = path.count("/")
			cutDir = slashCount - 2;
			cmd = "mkdir -p check && cd check && wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDir -R 'index.html*' -r --no-parent --directory-prefix=\$PWD/${input} ${path}  && cd .."
		} else if (path.indexOf('gs:') > -1 || path.indexOf('GS:') >-1){
			if (type == "folder"){
				cmd = "mkdir -p check ${workDir}/${input} && cd check && gsutil rsync -r ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			} else {
				cmd = "mkdir -p check && cd check && gsutil cp ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			}
		} else if (path.indexOf('/') == -1){
			cmd = ""
		}
}
	return [cmd,input]
}
if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.custom_additional_genome){params.custom_additional_genome = ""} 
if (!params.Metadata){params.Metadata = ""} 
if (!params.custom_additional_gtf){params.custom_additional_gtf = ""} 
if (!params.mask_gtf){params.mask_gtf = ""} 
if (!params.db_feather){params.db_feather = ""} 
if (!params.motif_db){params.motif_db = ""} 
if (!params.tf_lists){params.tf_lists = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)
ch_empty_file_3 = file("$baseDir/.emptyfiles/NO_FILE_3", hidden:true)
ch_empty_file_4 = file("$baseDir/.emptyfiles/NO_FILE_4", hidden:true)

if (params.reads){
Channel
	.fromFilePairs( params.reads,checkExists:true , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 ) 
	.set{g_1_1_g_51}
  } else {  
	g_1_1_g_51 = Channel.empty()
 }

Channel.value(params.mate).set{g_2_2_g_51}
(g_2_1_g_5) = [g_2_2_g_51]
g_18_2_g17_58 = params.custom_additional_genome && file(params.custom_additional_genome, type: 'any').exists() ? file(params.custom_additional_genome, type: 'any') : ch_empty_file_1
g_41_1_g36_0 = params.Metadata && file(params.Metadata, type: 'any').exists() ? file(params.Metadata, type: 'any') : ch_empty_file_1
g_48_3_g17_58 = params.custom_additional_gtf && file(params.custom_additional_gtf, type: 'any').exists() ? file(params.custom_additional_gtf, type: 'any') : ch_empty_file_2
g_59_2_g57_1 = params.mask_gtf && file(params.mask_gtf, type: 'any').exists() ? file(params.mask_gtf, type: 'any') : ch_empty_file_1
g_576_1_g_572 = file(params.db_feather, type: 'any')
g_577_2_g_572 = file(params.motif_db, type: 'any')
g_578_3_g_572 = file(params.tf_lists, type: 'any')

//* @style @array:{bcl_directory,mkfastq_sampleSheet} @multicolumn:{bcl_directory,mkfastq_sampleSheet}
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 2
}
//* autofill

process mkfastq_prep {


output:
 val mergedList  ,emit:g_53_bcl00_g_51 

when:
(params.run_mkfastq && (params.run_mkfastq == "yes")) || !params.run_mkfastq

exec:
mergedList = []
bcl_directory = params.mkfastq_prep.bcl_directory
mkfastq_sampleSheet = params.mkfastq_prep.mkfastq_sampleSheet
for (i = 0; i <bcl_directory.size(); i++) {
   sampleSheet = mkfastq_sampleSheet && mkfastq_sampleSheet[i] && file(mkfastq_sampleSheet[i], type: 'any').exists() ? file(mkfastq_sampleSheet[i]) : ch_empty_file_2
   mergedList.add(tuple(file(bcl_directory[i]), sampleSheet))
}


}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process mkfastq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${bcl_files}_reports$/) "mkfastq/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${bcl_files}_laneBarcode.html$/) "mkfastq/$filename"}
input:
 tuple file(bcl_files), file(sampleSheet)
 tuple val(name),file(reads)
 val mate

output:
 tuple file("fastq_final/*.{R1,R2}.fastq.gz")  ,emit:g_51_reads00_g_5 
 path "${bcl_files}_reports"  ,emit:g_51_outputDir11 
 path "${bcl_files}_laneBarcode.html" ,optional:true  ,emit:g_51_outputHTML22 

when:
params.run_mkfastq == "yes"

script:
cellranger_mkfastq_parameters = params.mkfastq.cellranger_mkfastq_parameters
SampleSheetCSV =  sampleSheet.name.startsWith('NO_FILE') ? "${bcl_files}/SampleSheet.csv": sampleSheet
"""	
cellranger mkfastq --id=mkfastq --run=${bcl_files} --csv=${SampleSheetCSV} --output-dir=fastq --rc-i2-override=true --delete-undetermined ${cellranger_mkfastq_parameters}
mv fastq ${bcl_files}_fastq
mkdir -p ${bcl_files}_reports fastq_final
# find "./${bcl_files}_fastq" -type f \\( -name "*_R1_001.fastq.gz" -o -name "*_R2_001.fastq.gz" \\) -exec mv {} "./fastq_final/." \\;
# Find all files ending with _R1_001.fastq.gz or _R2_001.fastq.gz
# sort -z needed for cat operation, so that reads belong to same sample(different lanes) can be merged
find "./${bcl_files}_fastq" -type f \\( -name "*_R1_001.fastq.gz" -o -name "*_R2_001.fastq.gz" \\) -print0 | sort -z |
while IFS= read -r -d '' file; do
	echo \$file
	base_name=\$(basename \${file})
    # Extract the common part of the file name
    common_part="\${base_name%_L00*_R*_001.fastq.gz}"

    # Determine if it's R1 or R2
    if [[ \$file == *_R1_001.fastq.gz ]]; then
        new_name="\${common_part}.R1.fastq.gz"
    elif [[ \$file == *_R2_001.fastq.gz ]]; then
        new_name="\${common_part}.R2.fastq.gz"
    fi
	echo \$new_name
    # rename them so groovy baseName() can detect sample name easily

    destination_file="fastq_final/\${new_name}"
    echo "source_file:\$file - destination_file: \$destination_file"
    if [ -f "\$destination_file" ]; then
        # Destination file exists, append the content of source file to it
        cat "\$file" >> "\$destination_file" && rm "\$file"
    else
        # Destination file doesn't exist, perform mv operation
        mv "\$file" "\$destination_file"
    fi
done
ls fastq_final


mv ${bcl_files}_fastq/Reports ${bcl_files}_reports/.
mv ${bcl_files}_fastq/Stats ${bcl_files}_reports/.
cp ${bcl_files}_reports/Reports/html/*/all/all/all/laneBarcode.html ${bcl_files}_laneBarcode.html
"""

}

//* params.gtf =  ""  //* @input
//* params.genome =  ""  //* @input
//* params.commondb =  ""  //* @input
//* params.genome_source =  ""  //* @input
//* params.gtf_source =  ""  //* @input
//* params.commondb_source =  ""  //* @input @optional

def downFile(path, task){
	println workDir
    if (path.take(1).indexOf("/") == 0){
      target=path
      if (task.executor == "awsbatch" || task.executor == "google-batch") {
      	a=file(path)
    	fname = a.getName().toString()
    	target = "${workDir}/${fname}"
    	if (!file(target).exists()){
    		a.copyTo(workDir)
    	}
      }
    } else {
      a=file(path)
      fname = a.getName().toString()
      target = "${workDir}/${fname}"
      if (!file(target).exists()){
    		a.copyTo(workDir)
      } 
    }
    return target
}

def getLastName (str){
	if (str.indexOf("/") > -1){
		return  str.substring(str.lastIndexOf('/')+1,str.length())
	} 
	return ""
}

process Check_and_Build_Module_Check_Genome_GTF {


output:
 path "${newNameFasta}"  ,emit:g17_21_genome00_g17_58 
 path "${newNameGtf}"  ,emit:g17_21_gtfFile10_g17_57 

container 'quay.io/viascientific/pipeline_base_image:1.0'

when:
params.run_Download_Genomic_Sources == "yes"

script:
genomeSource = !file("${params.genome}").exists() ? params.genome_source : params.genome
genomeName = getLastName(genomeSource)

gtfSource = !file("${params.gtf}").exists() ? params.gtf_source : params.gtf
gtfName = getLastName(gtfSource)


newNameGtf = gtfName
newNameFasta = genomeName
if (gtfName.contains('.gz')) { newNameGtf =  newNameGtf - '.gz'  } 
if (genomeName.contains('.gz')) { newNameFasta =  newNameFasta - '.gz'  } 

runGzip = ""
if (gtfName.contains('.gz') || genomeName.contains('.gz')) {
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} 

slashCountGenome = params.genome_source.count("/")
cutDirGenome = slashCountGenome - 3;

slashCountGtf = params.gtf_source.count("/")
cutDirGtf = slashCountGtf - 3;

"""
if [ ! -e "${params.genome_source}" ] ; then
    echo "${params.genome_source} not found"
	if [[ "${params.genome_source}" =~ "s3" ]]; then
		echo "Downloading s3 path from ${params.genome_source}"
		aws s3 cp ${params.genome_source} ${workDir}/${genomeName} && ln -s ${workDir}/${genomeName} ${genomeName}
	elif [[ "${params.genome_source}" =~ "gs" ]]; then
		echo "Downloading gs path from ${params.genome_source}"
		gsutil cp  ${params.genome_source} ${genomeName}
	else
		echo "Downloading genome with wget"
		wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDirGenome -R 'index.html*' -r --no-parent  ${params.genome_source}
	fi

else 
	ln -s ${params.genome_source} ${genomeName}
fi

if [ ! -e "${params.gtf_source}" ] ; then
    echo "${params.gtf_source} not found"
	if [[ "${params.gtf_source}" =~ "s3" ]]; then
		echo "Downloading s3 path from ${params.gtf_source}"
		aws s3 cp  ${params.gtf_source} ${workDir}/${gtfName} && ln -s ${workDir}/${gtfName} ${gtfName}
	elif [[ "${params.gtf_source}" =~ "gs" ]]; then
		echo "Downloading gs path from ${params.gtf_source}"
		gsutil cp  ${params.gtf_source} ${gtfName}
	else
		echo "Downloading gtf with wget"
		wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDirGtf -R 'index.html*' -r --no-parent  ${params.gtf_source}
	fi

else 
	ln -s ${params.gtf_source} ${gtfName}
fi

$runGzip

"""




}


process Check_and_Build_Module_convert_gtf_attributes {

input:
 path gtf

output:
 path "out/${gtf}"  ,emit:g17_57_gtfFile01_g17_58 

when:
params.replace_geneID_with_geneName == "yes"

shell:
'''
#!/usr/bin/env perl 

## Replace gene_id column with gene_name column in the gtf file
## Also check if any transcript_id defined in multiple chromosomes.
system("mkdir out");

open(OUT1, ">out/!{gtf}");
open(OUT2, ">notvalid_!{gtf}");
my %transcipt;
my $file = "!{gtf}";
open IN, $file;
while( my $line = <IN>)  {
    chomp;
    @a=split("\\t",$line);
    @attr=split(";",$a[8]);
    my %h;
    for my $elem (@attr) {
        ($first, $rest) = split ' ', $elem, 2;
        $h{$first} = $rest.";";
    }
    my $geneId = "";
    my $transcript_id = "";
    if (exists $h{"gene_name"}){
        $geneId = $h{"gene_name"};
    } elsif (exists $h{"gene_id"}){
        $geneId = $h{"gene_id"};
    }
    if (exists $h{"transcript_id"}){
        $transcript_id = $h{"transcript_id"};
    } elsif (exists $h{"transcript_name"}){
        $transcript_id = $h{"transcript_name"};
    } elsif (exists $h{"gene_id"}){
        $transcript_id = $h{"gene_id"};
    }
    if ($geneId ne "" && $transcript_id ne ""){
        ## check if any transcript_id defined in multiple chromosomes.
        if (exists $transcipt{$transcript_id}){
             if ($transcipt{$transcript_id} ne $a[0]){
               print OUT2 "$transcript_id: $transcipt{$transcript_id} vs $a[0]\\n";
                next;
                }
        } else {
             $transcipt{$transcript_id} = $a[0];
        }
        $a[8]=join(" ",("gene_id",$geneId,"transcript_id",$transcript_id));
        print OUT1 join("\\t",@a), "\\n";
    }  else {
        print OUT2 "$line";
    }
}
close OUT1;
close OUT2;
close IN;
'''
}


process Check_and_Build_Module_Add_custom_seq_to_genome_gtf {

input:
 path genome
 path gtf
 path custom_fasta
 path custom_gtf

output:
 path "${genomeName}_custom.fa"  ,emit:g17_58_genome00_g17_52 
 path "${gtfName}_custom_sorted.gtf"  ,emit:g17_58_gtfFile10_g17_53 

container 'quay.io/viascientific/custom_sequence_to_genome_gtf:1.0'

when:
params.add_sequences_to_reference == "yes"

script:
genomeName = genome.baseName
gtfName = gtf.baseName
is_custom_genome_exists = custom_fasta.name.startsWith('NO_FILE') ? "False" : "True" 
is_custom_gtf_exists = custom_gtf.name.startsWith('NO_FILE') ? "False" : "True" 
"""
#!/usr/bin/env python 
import requests
import os
import pandas as pd
import re
import urllib
from Bio import SeqIO

def add_to_fasta(seq, sqid, out_name):
	new_line = '>' + sqid + '\\n' + seq + '\\n'
	with open(out_name + '.fa', 'a') as f:
		f.write(new_line)

def createCustomGtfFromFasta(fastaFile, outCustomGtfFile):

    fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
    with open(outCustomGtfFile, "w") as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            last = len(sequence)
            line1 = "{gene}\\tKNOWN\\tgene\\t{first}\\t{last}\\t.\\t+\\t.\\tgene_id \\"{gene}\\"; gene_version \\"1\\"; gene_type \\"protein_coding\\"; gene_source \\"KNOWN\\"; gene_name \\"{gene}\\"; gene_biotype \\"protein_coding\\"; gene_status \\"KNOWN\\"; level 1;".format(gene=name, first="1", last=last)
            line2 = "{gene}\\tKNOWN\\ttranscript\\t{first}\\t{last}\\t.\\t+\\t.\\tgene_id \\"{gene}\\"; gene_version \\"1\\"; transcript_id \\"{gene}_trans\\"; transcript_version \\"1\\"; gene_type \\"protein_coding\\"; gene_source \\"KNOWN\\"; transcript_source \\"KNOWN\\"; gene_status \\"KNOWN\\"; gene_name \\"{gene}\\"; gene_biotype \\"protein_coding\\"; transcript_type \\"protein_coding\\"; transcript_status \\"KNOWN\\"; transcript_name \\"{gene}_1\\"; level 1; tag \\"basic\\"; transcript_biotype \\"protein_coding\\"; transcript_support_level \\"1\\";".format(gene=name, first="1", last=last)
            line3 = "{gene}\\tKNOWN\\texon\\t{first}\\t{last}\\t.\\t+\\t.\\tgene_id \\"{gene}\\"; gene_version \\"1\\"; transcript_id \\"{gene}_trans\\"; transcript_version \\"1\\"; exon_number 1; gene_type \\"protein_coding\\"; gene_source \\"KNOWN\\"; transcript_source \\"KNOWN\\"; gene_status \\"KNOWN\\"; gene_name \\"{gene}\\"; gene_biotype \\"protein_coding\\"; transcript_type \\"protein_coding\\"; transcript_status \\"KNOWN\\"; transcript_biotype \\"protein_coding\\"; transcript_name \\"{gene}_1\\"; exon_number 1; exon_id \\"{gene}.1\\"; level 1; tag \\"basic\\"; transcript_support_level \\"1\\";".format(gene=name, first="1", last=last)
            out_file.write("{}\\n{}\\n{}\\n".format(line1, line2, line3))

	
os.system('cp ${genomeName}.fa ${genomeName}_custom.fa')  
os.system('cp ${gtfName}.gtf ${gtfName}_custom.gtf')  

if ${is_custom_genome_exists}:
	os.system("tr -d '\\r' < ${custom_fasta} > ${custom_fasta}_tmp && rm ${custom_fasta} && mv ${custom_fasta}_tmp ${custom_fasta}")
	os.system('cat ${custom_fasta} >> ${genomeName}_custom.fa')
	if ${is_custom_gtf_exists}:
		os.system("tr -d '\\r' < ${custom_gtf} > ${custom_gtf}_tmp && rm ${custom_gtf} && mv ${custom_gtf}_tmp ${custom_gtf}")
		os.system("mv ${custom_gtf} ${custom_fasta}.gtf")
	else:
		createCustomGtfFromFasta("${custom_fasta}", "${custom_fasta}.gtf")
	os.system('cat ${custom_fasta}.gtf >> ${gtfName}_custom.gtf')

	
os.system('samtools faidx ${genomeName}_custom.fa')
os.system('igvtools sort ${gtfName}_custom.gtf ${gtfName}_custom_sorted.gtf')
os.system('igvtools index ${gtfName}_custom_sorted.gtf')

"""
}

//* params.gtf2bed_path =  ""  //* @input
//* params.bed =  ""  //* @input

process Check_and_Build_Module_Check_BED12 {

input:
 path gtf

output:
 path "${gtfName}.bed"  ,emit:g17_53_bed03_g17_54 

container "${ params.IMAGE_BASE ? "${params.IMAGE_BASE}/rnaseq:4.0" : "quay.io/viascientific/rnaseq:4.0" }"

when:
params.run_Download_Genomic_Sources == "yes"

script:
gtfName  = gtf.baseName
beddir = ""
if (params.bed.indexOf('/') > -1){
	beddir  = params.bed.substring(0, params.bed.lastIndexOf('/')) 
}
"""

if [ ! -e "${params.bed}" ] ; then
    echo "${params.bed} not found"
    perl ${params.gtf2bed_path} $gtf > ${gtfName}.bed
else 
	cp -n ${params.bed} ${gtfName}.bed
fi
if [ "${beddir}" != "" ] ; then
	mkdir -p ${beddir}
	cp -n ${gtfName}.bed ${params.bed} 
fi
"""




}

//* params.gtf2bed_path =  ""  //* @input
//* params.genome_sizes =  ""  //* @input

process Check_and_Build_Module_Check_chrom_sizes_and_index {

input:
 path genome

output:
 path "${genomeName}.chrom.sizes"  ,emit:g17_52_genomeSizes02_g17_54 

when:
params.run_Download_Genomic_Sources == "yes"

script:
genomeName  = genome.baseName
genome_sizes_dir = ""
if (params.genome_sizes.indexOf('/') > -1){
	genome_sizes_dir  = params.genome_sizes.substring(0, params.genome_sizes.lastIndexOf('/')) 
}

"""
if [ ! -e "${params.genome_sizes}" ] ; then
    echo "${params.genome_sizes} not found"
    cat ${genome} | awk '\$0 ~ ">" {print c; c=0;printf substr(\$1,2,100) "\\t"; } \$0 !~ ">" {c+=length(\$0);} END { print c; }' > ${genomeName}.chrom.sizes
    ##clean first empty line
    sed -i '1{/^\$/d}' ${genomeName}.chrom.sizes
    if [ "${genome_sizes_dir}" != "" ] ; then
    	mkdir -p ${genome_sizes_dir}
		cp -n ${genomeName}.chrom.sizes ${params.genome_sizes} 
	fi
else 
	cp ${params.genome_sizes} ${genomeName}.chrom.sizes
fi

"""




}


process Check_and_Build_Module_check_files {

input:
 path gtf
 path genome
 path genomeSizes
 path bed

output:
 path "*/${gtf2}" ,optional:true  ,emit:g17_54_gtfFile01_g_19 
 path "*/${genome2}" ,optional:true  ,emit:g17_54_genome10_g_19 
 path "*/${genomeSizes2}" ,optional:true  ,emit:g17_54_genomeSizes22 
 path "*/${bed2}" ,optional:true  ,emit:g17_54_bed33 

container 'quay.io/viascientific/pipeline_base_image:1.0'
stageInMode 'copy'

script:
(cmd1, gtf2) = pathChecker(gtf, params.gtf, "file")
(cmd2, genome2) = pathChecker(genome, params.genome, "file")
(cmd3, genomeSizes2) = pathChecker(genomeSizes, params.genome_sizes, "file")
(cmd4, bed2) = pathChecker(bed, params.bed, "file")
"""
$cmd1
$cmd2
$cmd3
$cmd4
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 3
    $MEMORY = 32
}
//* platform
//* platform
//* autofill

process cell_ranger_mkref {

input:
 path genome
 path gtf

output:
 path "ref"  ,emit:g_19_reference00_g_20 

when:
(params.run_mkref && (params.run_mkref == "yes")) || !params.run_mkref

script:
optional_mkgtf_filtering_parameters = params.cell_ranger_mkref.optional_mkgtf_filtering_parameters
"""
if [ "${params.gtf_type}" == "gencode" ]; then
	# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.
	BIOTYPE_PATTERN="(protein_coding|lncRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene)"
	GENE_PATTERN="gene_type \\"\${BIOTYPE_PATTERN}\\""
	TX_PATTERN="transcript_type \\"\${BIOTYPE_PATTERN}\\""
	READTHROUGH_PATTERN="tag \\"readthrough_transcript\\""
	PAR_PATTERN="tag \\"PAR\\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
	cat "${gtf}" \
    | awk '\$3 == "transcript"' \
    | grep -E "\$GENE_PATTERN" \
    | grep -E "\$TX_PATTERN" \
    | grep -Ev "\$READTHROUGH_PATTERN" \
    | grep -Ev "\$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\\1/' | sort | uniq > "gene_allowlist"
    
    ## IF gene_type or transcript_type columns are not found, it means that line does not belong to gencode. Keep those gene_ids. It might come from custom gtf files.
    cat ${gtf}  | awk '!/gene_type/ || !/transcript_type/ {print \$0}' |awk '\$3 == "transcript"' | grep -Ev "^#" | sed -E 's/.*(gene_id "[^"]+").*/\\1/' | sort | uniq >> "gene_allowlist"

	# Filter the GTF file based on the gene allowlist
	# Copy header lines beginning with "#"
	grep -E "^#" "${gtf}" > filtered_${gtf}
	# Filter to the gene allowlist
	grep -Ff "gene_allowlist" "${gtf}" >> filtered_${gtf}
	
	cellranger mkgtf filtered_${gtf} mkgtf_${gtf} ${optional_mkgtf_filtering_parameters}
else 
	cellranger mkgtf $gtf mkgtf_${gtf} --attribute=gene_biotype:protein_coding ${optional_mkgtf_filtering_parameters}
fi
export TENX_IGNORE_DEPRECATED_OS=1 && TENX_IGNORE_DEPRECATED_OS=1 cellranger mkref --genome=ref --fasta=${genome} --genes=mkgtf_${gtf}
"""

}

//* params.transcriptome =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process cellranger_ref_checker {

input:
 path transcriptome

output:
 path "*/${transcriptome}"  ,emit:g_20_reference02_g_5 

container 'quay.io/viascientific/pipeline_base_image:1.0'
stageInMode 'copy'

script:
(cmd1, transcriptome) = pathChecker(transcriptome, params.transcriptome, "folder")
"""
$cmd1
"""
}

//* params.cellranger_path =  ""  //* @input
cell_ranger_count_parameters = params.Count.cell_ranger_count_parameters
cell_ranger_count_create_bam = params.Count.cell_ranger_count_create_bam
expected_cells = params.Count.expected_cells
force_cells = params.Count.force_cells
//* params.transcriptome =  ""  //*  @input 
chemistry = params.Count.chemistry
publish_barcode_features_matrix_files_into_reports = params.Count.publish_barcode_features_matrix_files_into_reports

def getLastDirName(row){
   firstSec = row.toString().substring(0,row.toString().lastIndexOf('/'))
   secondSec = firstSec.substring(firstSec.lastIndexOf('/')+1, firstSec.length())
return secondSec
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process Count {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_outs$/) "cellranger_count/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_web_summary.html$/) "count_web_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_filtered_feature_bc_matrix.h5$/) "Filtered_Count_Matrix_h5/$filename"}
input:
 tuple val(name), file(reads)
 val mate
 path ref

output:
 path "${name}_outs"  ,emit:g_5_outputDir00_g57_5 
 path "${name}_web_summary.html"  ,emit:g_5_outputHTML11 
 tuple val(name), file("${name}_filtered_feature_bc_matrix") ,optional:true  ,emit:g_5_outputDir22 
 tuple val(name), file("${name}_raw_feature_bc_matrix")  ,emit:g_5_outputDir33 
 tuple val(name), file("${name}_filtered_feature_bc_matrix.h5")  ,emit:g_5_h5_file40_g_49 
 tuple val(name), file("${name}_raw_feature_bc_matrix.h5")  ,emit:g_5_h5_file51_g_49 

disk { 1000.GB * task.attempt }
when:
params.run_Cell_Ranger_Count == "yes"

script:

sample = name
nameAll = reads.toString()
nameArray = nameAll.split(' ')
if (mate == "pair"){
    read1 = nameArray[0]
    read2 = nameArray[1]
    if (nameAll.contains('.gz')) {
        runGzip = ''
    } else {
    	read1raw = read1
        read2raw = read2
        read1 = read1raw + ".gz"
        read2 = read2raw + ".gz"
        runGzip = "gzip -f ${read1raw} ${read2raw}"
    }
    mvReads = "mv " +read1 +" reads/"+ name + "_S1_L001_R1_001.fastq.gz && mv " +read2 +" reads/"+ name + "_S1_L001_R2_001.fastq.gz"  
} else {
    read1 = nameArray[0]
    if (nameAll.contains('.gz')) {
        runGzip = ''
    } else {
    	read1raw = read1
        read1 = read1raw + ".gz"
        runGzip = "gzip -f ${read1raw}"
    }
    mvReads = "mv " +read1 +" reads/"+ name + "_S1_L001_R1_001.fastq.gz"
}
localmem = task.memory.toGiga()
expected_cells_text = (expected_cells.toString() != "") ? "--expect-cells "+ expected_cells : ""
force_cells_text = (force_cells.toString() != "") ? "--force-cells "+ force_cells : ""

"""
$runGzip
mkdir reads
$mvReads

## cell ranger expect this pattern ${name}_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
${params.cellranger_path} count --localcores=${task.cpus} --localmem=${localmem} --id=$name \
                   --transcriptome=\$PWD/${ref} \
                   --fastqs=reads \
                   --sample=$sample \
                   --create-bam=${cell_ranger_count_create_bam} \
                   $expected_cells_text \
                   $force_cells_text \
                   --chemistry=$chemistry $cell_ranger_count_parameters
                   

mv ${name}/outs ${name}_outs
rm -rf ${name}
cp ${name}_outs/web_summary.html ${name}_web_summary.html
cp ${name}_outs/*.bam ${name}_possorted_genome_bam.bam

cp -r ${name}_outs/raw_feature_bc_matrix ${name}_raw_feature_bc_matrix
cp ${name}_outs/*filtered*.h5 ${name}_filtered_feature_bc_matrix.h5
cp ${name}_outs/*raw*.h5 ${name}_raw_feature_bc_matrix.h5



if [ "${publish_barcode_features_matrix_files_into_reports}" == "yes" ]; then
	mv ${name}_outs/filtered_feature_bc_matrix ${name}_filtered_feature_bc_matrix
fi 

"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 200
}
//* platform
//* platform
//* autofill

process Ambient_RNA_QC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_Ambient_RNA_QC.html$/) "Ambient_RNA_Removal_Report/$filename"}
input:
 tuple val(name), file("${name}_filtered_feature_bc_matrix.h5")
 tuple val(name), file("${name}_raw_feature_bc_matrix.h5")

output:
 tuple val(name), file("${name}_corrected_feature_bc_matrix.h5"), file("*raw_placeholder.h5"), file("*tags_placeholder.tsv")  ,emit:g_49_h5_file00_g36_0 
 path "*_Ambient_RNA_QC.html"  ,emit:g_49_outputFileHTML11 

container 'quay.io/viascientific/ambient_rna_removal:1.0'
label 'Ambient_RNA_QC'

shell:
Ambient_RNA_Removal = (params.Ambient_RNA_Removal == "yes") ? "TRUE" : "FALSE"

'''
#!/usr/bin/env perl

# --- create placeholder files for downstream tuple compatibility ---
# If upstream didn’t provide these, make sure they exist (empty is fine)
system("bash -lc ': > !{name}_raw_placeholder.h5'") == 0 or die "failed to create raw placeholder\n";
system("bash -lc ': > !{name}_tags_placeholder.tsv'") == 0 or die "failed to create tags placeholder\n";
# ------------------------------------------------------------------


my $script = <<'EOF';
---
title: "Ambient RNA QC Report"
author: "Zhaorong Li"
output: 
  html_document:
    toc: true
    toc_float:
      toc_collapsed: true
      toc_depth: 3
    number_sections: true
    fig_caption: yes
    theme: cerulean
    code_folding: hide
editor_options: 
  chunk_output_type: console
params:
  Ambient_RNA_Removal: FALSE
  RawCountMatrix: ""
  FilteredCountMatrix: ""
  Outputname: ""

---


```{r setup, include=FALSE}

suppressPackageStartupMessages({
library(Seurat)
library(DropletUtils)
library(celda)
library(SingleCellExperiment)
})

```

# Read in samples

```{r read in samples}

RawData=Read10X_h5(params$RawCountMatrix)
RawData=CreateSeuratObject(RawData)
RawData = as.SingleCellExperiment(RawData)
RawData <- RawData[, colSums(counts(RawData)) > 0]

FilteredData=Read10X_h5(params$FilteredCountMatrix)
FilteredData=CreateSeuratObject(FilteredData)
FilteredData = as.SingleCellExperiment(FilteredData)
FilteredData <- FilteredData[, colSums(counts(FilteredData)) > 0]

```

# deconX

```{r deconX, fig.width=5,fig.height=5}

FilteredData <- celda::decontX(FilteredData,background=RawData)


```

# Plot out Contamination in the data

``` {r plotDecontXContamination, fig.width=5,fig.height=5}
plotDecontXContamination(FilteredData)
#plot(b$contamination)

if (params$Ambient_RNA_Removal) {
	counts=decontXcounts(FilteredData)
	counts=round(counts)
	write10xCounts(params$Outputname,counts,barcodes = colnames(counts),gene.id = rownames(counts),gene.symbol = rownames(counts),type = 'HDF5')
	
} else {
	counts=Read10X_h5(params$FilteredCountMatrix)
	write10xCounts(params$Outputname,counts,barcodes = colnames(counts),gene.id = rownames(counts),gene.symbol = rownames(counts),type = 'HDF5')
}

```


EOF

open OUT, ">!{name}_Ambient_RNA_QC.rmd";
print OUT $script;
close OUT;

runCommand("Rscript -e 'rmarkdown::render(\\"!{name}_Ambient_RNA_QC.rmd\\",\\"html_document\\", output_file = \\"!{name}_Ambient_RNA_QC.html\\",params = list(Ambient_RNA_Removal=as.logical(\\"!{Ambient_RNA_Removal}\\"),RawCountMatrix=\\"!{name}_raw_feature_bc_matrix.h5\\",FilteredCountMatrix=\\"!{name}_filtered_feature_bc_matrix.h5\\",Outputname=\\"!{name}_corrected_feature_bc_matrix.h5\\"))'");

sub runCommand {
            my ($com) = @_;
            my $error = system($com);
            if   ($error) { die "Command failed: $error $com\\n"; }
            else          { print "Command successful: $com\\n"; }
          }

'''




}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 200
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Quality_Control_and_Filtering {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_filtering_report.html$/) "QC_Reports/$filename"}
input:
 tuple val(name), file(input_file), file(raw), file(tags)
 path metadata

output:
 path "${name}.rds"  ,emit:g36_0_rdsFile00_g36_14 
 tuple val(name),file("${name}_filtering_report.html")  ,emit:g36_0_outputFileHTML11 
 path "${name}_filter_summary.tsv"  ,emit:g36_0_outFileTSV20_g36_34 

container "quay.io/viascientific/scrna_seurat:2.1"

when:
(params.run_scRNA_Analysis && (params.run_scRNA_Analysis == "yes")) || !params.run_scRNA_Analysis || params.run_pySCENIC == "yes"

script:

remove_mitochondiral_genes = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.remove_mitochondiral_genes
remove_ribosomal_genes = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.remove_ribosomal_genes

min_genes = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.min_genes
max_genes = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.max_genes
min_UMIs = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.min_UMIs
max_UMIs = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.max_UMIs
percent_mitochondrial = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.percent_mitochondrial
percent_ribosomal = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.percent_ribosomal

doublet_removal = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.doublet_removal
doublet_percentage = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.doublet_percentage
doublet_removal_tool = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.doublet_removal_tool

normalization_method = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.normalization_method
variable_features = params.scRNA_Analysis_Module_Quality_Control_and_Filtering.variable_features

remove_mitochondiral_genes_arg = remove_mitochondiral_genes == 'true' ? '--remove-mitochondrial-genes' : ''
remove_ribosomal_genes_arg = remove_ribosomal_genes == 'true' ? '--remove_ribosomal_genes' : ''

//* @style @condition:{doublet_removal="true", doublet_percentage, doublet_removal_tool},{doublet_removal="false"} @multicolumn:{remove_mitochondiral_genes, remove_ribosomal_genes}, {min_genes, max_genes}, {min_UMIs, max_UMIs}, {percent_mitochondrial, percent_ribosomal}, {doublet_percentage, doublet_removal_tool}, {normalization_method, variable_features}

doublet_percentage_arg = doublet_percentage ? "--doublet-percentage ${doublet_percentage}" : ""

"""

tags_arg=""
if [ -s ${tags} ] && [ \$(wc -l < ${tags}) -gt 1 ]; then
    echo "Valid tags file detected - including in analysis"
    tags_arg="--tags-file ${tags}"
else
    tags_arg="--tags-file ''"
fi

well_arg=""
if [ -s "${name}_raw_placeholder.h5" ]; then
  echo "Valid well file detected - including in analysis"
  well_arg="--input-well ${name}_raw_placeholder.h5"
fi

echo \$tags_arg

build_QC_report.py --output-prefix ${name} --input-file ${input_file} \$well_arg \$tags_arg \
--metadata-file ${metadata} \
--min-genes ${min_genes} --max-genes ${max_genes} \
--min-UMIs ${min_UMIs} --max-UMIs ${max_UMIs} \
--percent-mitochondrial-cutoff ${percent_mitochondrial} --percent-ribosomal-cutoff ${percent_ribosomal} ${remove_mitochondiral_genes_arg} ${remove_ribosomal_genes_arg} \
--variable-features ${variable_features} --normalization-method ${normalization_method} \
--doublet-removal ${doublet_removal} ${doublet_percentage_arg} --doublet-tool ${doublet_removal_tool}
"""



}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Merge_Seurat_Objects {

input:
 path seurat_obj

output:
 path "merged_filtered_seurat.rds"  ,emit:g36_14_rdsFile00_g36_17 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
"""
#!/usr/bin/env Rscript

library(Seurat)
list_of_samples <- list.files(pattern = "*.rds")

if (length(list_of_samples)==1) {
	list_of_seurat=list()
	obj = readRDS(list_of_samples[1])
	list_of_seurat[[list_of_samples[1]]]=obj
} else {
list_of_seurat <- list()
for(i in 1:length(list_of_samples)){
  # print name
  print(list_of_samples[i])
  list_of_seurat[[i]] <- readRDS(list_of_samples[i])
}

}
saveRDS(list_of_seurat, file="merged_filtered_seurat.rds")

"""


}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 140
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction {

input:
 path seurat_object

output:
 path "Reduced_and_Corrected.rds"  ,emit:g36_17_rdsFile00_g36_19 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
varFeatures = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.varFeatures
selmethod = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.selmethod
Batch_Effect_Correction = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.Batch_Effect_Correction
WNN = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.WNN

//* @style @multicolumn:{varFeatures, selmethod},{Batch_Effect_Correction, WNN}

"""
#!/usr/bin/env Rscript

# libraries
library(Seurat)
library(dplyr)
#install.packages("harmony",repos = "http://cran.us.r-project.org")
library(harmony)

selmethod <- "${selmethod}"
varFeatures <- "${varFeatures}"

Data=readRDS("${seurat_object}")
Multi_sample=0
if (length(Data)==1) {
	Data=Data[[1]]
	if (DefaultAssay(Data)=="SCT"){
		Data=RunPCA(Data,npcs=100)
	} else {
		Data <- FindVariableFeatures(Data,selection.method=selmethod,nfeatures=as.numeric(varFeatures))
		if (all(Data[["percent.mt"]]==0)) {
		Data <- ScaleData(Data)
		} else {
		Data <- ScaleData(Data,vars.to.regress="percent.mt")
		}
		Data=RunPCA(Data,npcs=100)
	}
} else {
Multi_sample=1
	if (DefaultAssay(Data[[1]])=="SCT") {
		variable.features=SelectIntegrationFeatures(object.list = Data, nfeatures = as.numeric(varFeatures))
		Data <- merge(Data[[1]],Data[-1])
		VariableFeatures(Data) <- variable.features

		if (all(Data[["percent.mt"]]==0)) {
			Data <- ScaleData(Data)
		} else {
			Data <- ScaleData(Data,vars.to.regress="percent.mt")
			}
		Data=RunPCA(Data,npcs=100)
		if (as.logical("${Batch_Effect_Correction}")){
		Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
		}
	} else {
		Data <- merge(Data[[1]],Data[-1])
		Data <- FindVariableFeatures(Data,selection.method=selmethod,nfeatures=as.numeric(varFeatures))
		if (all(Data[["percent.mt"]]==0)) {
			Data <- ScaleData(Data)
		} else {
			Data <- ScaleData(Data,vars.to.regress="percent.mt")
		}
		Data=RunPCA(Data,npcs=100)
		if (as.logical("${Batch_Effect_Correction}")){
		Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
		}

	}
}

if ("${WNN}"!="") {
original.assay=DefaultAssay(Data)

DefaultAssay(Data)="${WNN}"

Data=NormalizeData(Data,normalization.method = "CLR",margin=2)

VariableFeatures(Data)=rownames(Data)

Data=ScaleData(Data)

Data=RunPCA(Data,reduction.name = "wpca")

if (Multi_sample==1) {
	Data=RunHarmony(Data,group.by.vars = "sample",assay.use = "${WNN}",reduction = "wpca",reduction.save = "wharmony")

}

DefaultAssay(Data)=original.assay

}

saveRDS(Data,"Reduced_and_Corrected.rds")

"""


}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 140
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Clustering_and_Find_Markers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /final_report.html$/) "Final_Report/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "ClusterMarkers/$filename"}
input:
 path seurat_object

output:
 path "final_report.html"  ,emit:g36_19_outputHTML00 
 path "Final_Analysis.rds"  ,emit:g36_19_rdsFile10_g36_36 
 path "*.tsv"  ,emit:g36_19_outFileTSV22 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
min_resolution = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.min_resolution
max_resolution = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.max_resolution
num_pc = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.num_pc
find_markers_for_all_resolution = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.find_markers_for_all_resolution

find_markers_for_all_resolution_arg = find_markers_for_all_resolution == 'true' ? '--all-resolution-cluster-markers' : ''

//* @style @multicolumn:{min_resolution, max_resolution}, {num_pc, find_markers_for_all_resolution}

"""
build_clustering_and_find_markers.py --sample-path ${seurat_object} --min-resolution ${min_resolution} --max-resolution ${max_resolution} --num-pc ${num_pc} ${find_markers_for_all_resolution_arg}
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 4
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_filter_summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "filter_summary/$filename"}
input:
 path input_files

output:
 path "*"  ,emit:g36_34_outputFileHTML00 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
	
"""
build_filtration_report.py --input-dir .

mkdir output
mv by_criteria_summary.tsv output
mv filtration_summary_report.Rmd output
mv overall_filtration_summary.tsv output
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 20
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_sc_annotation {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.rds$/) "scViewer/$filename"}
input:
 path input_rds

output:
 path "*.rds"  ,emit:g36_36_rdsFile01_g36_40 

container 'quay.io/mustafapir/sc_annotation:1.0.1'

script:

tissue_type = params.scRNA_Analysis_Module_sc_annotation.tissue_type

species_map = [
    "human_hg38_gencode_v32_cellranger_v6"      : "human",
    "human_hg38_cellranger_GRCh38-2024-A"       : "human",
    "mouse_mm10_gencode_vm23_cellranger_v6"     : "mouse",
    "mouse_GRCm39_cellranger_GRCm39-2024-A"      : "mouse",
    "zebrafish_GRCz11plus_ensembl"              : "zebrafish",
    "d_melanogaster_BDGP6_32_ensembl_105_cellranger_v6": "d_melanogaster",
    "d_melanogaster_flybase_r6_45_cellranger_v6" : "d_melanogaster"
]

species = species_map.get(params.genome_build, "human")

"""
if [ ${params.run_annotation} = "yes" ]; then
  run_sctype.R \
    --input ${input_rds} \
    --output Final_Analysis_annotated.rds \
    --organism ${species} \
    --tissue "${tissue_type}" && \
  rm -f ${input_rds}
fi

"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_SCEtoLOOM {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /Data.loom$/) "LOOM_File/$filename"}
input:
 path seurat_object

output:
 path "Data.loom" ,optional:true  ,emit:g36_30_outputFileOut00_g_572 

container "quay.io/mustafapir/scrna_seurat:2.0.2"

script:
Generate_loom_file = params.scRNA_Analysis_Module_SCEtoLOOM.Generate_loom_file

"""
#!/usr/bin/env Rscript

#library
library(Seurat)
library(biomaRt)
 
if (as.logical("${Generate_loom_file}") || as.logical("${params.run_pySCENIC == 'yes'}")) {

Data=readRDS("${seurat_object}")

create_annotation_table <- function(organism = c("human", "mouse", "d_melanogaster")) {
  organism <- match.arg(organism)
  
  # Select ENSEMBL dataset depending on organism
  dataset <- switch(
    organism,
    human = "hsapiens_gene_ensembl",
    mouse = "mmusculus_gene_ensembl",
    d_melanogaster = "dmelanogaster_gene_ensembl"
  )
  
  # Connect to Ensembl
  ensembl <- useEnsembl(biomart = "genes", dataset = dataset)
  
  # Retrieve Ensembl gene ID, gene name, gene type
  genes <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
    mart = ensembl
  )
  
  # Ensure column names match your exact requested format
  colnames(genes) <- c("ensembl_id", "gene_name", "gene_type")
  
  return(genes)
}

# annotation=read.csv("https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/raw/main/example_input_files/gene_info_table.csv",header = T,row.names = 1)
annotation=create_annotation_table("${params.species}")

annotation=annotation[annotation\$gene_name%in%names(table(annotation\$gene_name))[table(annotation\$gene_name)==1],]
rownames(annotation)=annotation\$gene_name
metadata=Data@meta.data
matrix=Data@assays\$RNA@counts
matrix=matrix[rowSums(matrix)>0,]
matrix=matrix[intersect(rownames(matrix),rownames(annotation)),]
annotation=annotation[intersect(rownames(matrix),rownames(annotation)),]
matrix=matrix[rownames(matrix)[order(rownames(matrix),decreasing = F)],]
NewData=CreateSeuratObject(matrix,meta.data = metadata)
NewData[["RNA"]]@meta.features\$ensembl_id=annotation[rownames(NewData),"ensembl_id"]
NewData.loom <- SeuratDisk::as.loom(NewData, filename = "Data.loom", verbose = FALSE)

}
"""


}

//* params.db_feather =  ""  //* @input
//* params.motif_db =  ""  //* @input
//* params.tf_lists =  ""  //* @input


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 20
}
//* platform
//* platform
//* autofill

process pySCENIC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /pyscenic_out.zip$/) "pySCENIC_out/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /scenic_integrated.loom$/) "pySCENIC_loom/$filename"}
input:
 path loom_file
 path db_feather
 path motif_db
 path tf_lists

output:
 path "pyscenic_out.zip"  ,emit:g_572_zipFile00 
 path "scenic_integrated.loom"  ,emit:g_572_loom11 

container 'quay.io/viascientific/pyscenic:1.0.1'

when:
params.run_pySCENIC == "yes"

script:

threads = task.cpus
mask_dropouts = params.pySCENIC.mask_dropouts

mask_dropouts_option = mask_dropouts ? '--mask_dropouts' : ''
auc_threshold = params.pySCENIC.auc_threshold
"""

f_db_path=${db_feather}
f_db_names=\$(echo "\$f_db_path"/*.feather)

f_motif_path=${motif_db}
f_tfs=${tf_lists}

pyscenic grn ${loom_file} ${tf_lists} -o adjacencies.csv --num_workers ${threads}


pyscenic ctx adjacencies.csv \
    \$f_db_names \
    --annotations_fname ${motif_db} \
    --expression_mtx_fname ${loom_file} \
    --output regulons.csv \
    ${mask_dropouts_option} \
    --num_workers ${threads}

pyscenic aucell \
    ${loom_file} \
    regulons.csv \
    --output pyscenic_out.loom \
    --num_workers ${threads} \
    --auc_threshold ${auc_threshold}

integrate_pyscenic_output.py \
    -i ${loom_file} \
    -p pyscenic_out.loom \
    -o scenic_integrated.loom \
    --export_auc_csv aucell_matrix.csv

zip pyscenic_out.zip adjacencies.csv regulons.csv aucell_matrix.csv

"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 60
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Create_h5ad {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.h5ad$/) "H5AD_file/$filename"}
input:
 path seurat_obj

output:
 path "*.h5ad"  ,emit:g36_22_h5ad_file01_g57_12 

container "quay.io/viascientific/scrna_seurat:2.0"

script:
"""
#!/usr/bin/env Rscript

# libraries
library(Seurat)
library(SeuratDisk)

# read data
seurat_obj <- readRDS("${seurat_obj}")

for (var in colnames(seurat_obj@meta.data)) {
  if (is.factor(seurat_obj@meta.data[,var])) {
    seurat_obj@meta.data[,var]=as.character(seurat_obj@meta.data[,var])
  }
}

# save h5ad file
seu_name <- gsub(".rds","","${seurat_obj}")
SaveH5Seurat(seurat_obj, filename = paste0(seu_name,".h5Seurat"))
Convert(paste0(seu_name,".h5Seurat"), dest = "h5ad")
"""


}


process scRNA_Analysis_Module_tcr_analysis {

input:
 path final_rds
 path annotated_final_rds

output:
 path "tcr_analysis.rds"  ,emit:g36_40_rdsFile00 


when:
(params.run_tcr_analysis && (params.run_tcr_analysis == "yes")) && !vdj_contigs.toString().contains("NO_FILE")

script:

"""
if [[ "${annotated_final_rds}" == *NO_FILE* ]]; then
    mv ${final_rds} tcr_analysis.rds
else
    mv ${annotated_final_rds} tcr_analysis.rds
fi

"""
}


process RNA_Velocity_Module_prepare_input_velocyto {

input:
 path outs

output:
 tuple val("${new_name}"), file("output_files/input.bam"), file("output_files/input.bam.bai")  ,emit:g57_5_bamFile00_g57_1 
 tuple val("${new_name}"), file("output_files/input_barcodes.tsv.gz")  ,emit:g57_5_inputFileTsv11_g57_1 

when:
params.run_velocity == "yes"

script:

try {
	myVariable = bam
} catch (MissingPropertyException e) {
	bam = ""
}

try {
	myVariable = bai
} catch (MissingPropertyException e) {
	bai = ""
}

try {
	myVariable = barcodes
} catch (MissingPropertyException e) {
	barcodes = ""
}

new_name = outs.toString().startsWith('NO_FILE')
    ? name
    : outs.toString().replaceFirst(/_outs$/, '')

"""
mkdir -p output_files

if [[ ${outs} == NO_FILE* ]]; then
    mv ${bam} output_files/input.bam
    mv ${bai} output_files/input.bam.bai
    mv ${barcodes} output_files/input_barcodes.tsv.gz

else
    mv ${outs}/possorted_genome_bam.bam output_files/input.bam
    mv ${outs}/possorted_genome_bam.bam.bai output_files/input.bam.bai
    mv ${outs}/filtered_feature_bc_matrix/barcodes.tsv.gz output_files/input_barcodes.tsv.gz
fi
echo ${new_name}

"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 16
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process RNA_Velocity_Module_velocyto {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_output.loom$/) "loom_out/$filename"}
input:
 tuple val(name), file(bam), file(bai)
 tuple val(name), file(barcodes)
 path mask_gtf
 path gtf_velo

output:
 path "${name}_output.loom"  ,emit:g57_1_loom00_g57_12 

container "quay.io/biocontainers/velocyto.py:0.17.17--py310h581d4b6_7"

when:
params.run_velocity == "yes"

script:

mask_gtf_option = mask_gtf.name.startsWith('NO_FILE') ? "" : "-m ${mask_gtf}"

"""
echo ${name}
mkdir -p velocyto_out

cp ${bam} ./local_input.bam
velocyto run -b ${barcodes} ${mask_gtf_option} -o velocyto_out local_input.bam ${gtf_velo}

mv velocyto_out/*.loom ${name}_output.loom

"""
}


process RNA_Velocity_Module_process_anndata {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /processed_adata.h5ad$/) "scVelo_out/$filename"}
input:
 path loom_file
 path h5ad_file

output:
 path "processed_adata.h5ad"  ,emit:g57_12_h5ad00 

container 'quay.io/mustafapir/scvelo_shiny:1.1.0'

when:
params.run_velocity == "yes"

script:

"""

preprocess_anndata.py \
    --h5ad ${h5ad_file} \
    --loom ${loom_file} \
    --output 'processed_adata.h5ad'

"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 8
    $MEMORY = 60 
}
//* platform
//* platform
//* autofill

process CellChat2_create_cellchat_obj {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_cellchat.rds$/) "single_sample_analysis/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /signalling_scripts\/.*.rmd$/) "scripts/$filename"}
input:
 path seurat_obj

output:
 path "*_cellchat.rds"  ,emit:g580_1_rdsFile00 
 path "signalling_scripts/*.rmd"  ,emit:g580_1_rMarkdown11 
 path "cellchat_list.rds"  ,emit:g580_1_rdsFile20_g580_5 

container 'quay.io/viascientific/cellchat2:2.0.0'

when:
params.run_cellchat2 == "yes"

script:

// CellChat2 params
//* params.ident =  "sctype_classification"  //* @input @label:"Cell labels" @description:"Cell label column name"
//* params.db_type =  "all"  //* @dropdown @label:"Ligand-receptor type" @options:"all","except_nonprotein","Secreted Signaling","ECM-Receptor","Cell-Cell Contact","Non-protein Signaling" @description:"Subset CellChatDB by ligand-receptor type"
//* params.smooth =  "no"  //* @dropdown @label:"Project expression data onto PPI network?" @options:"yes","no" @description:"A diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors."
//* params.average_method =  "triMean"  //* @dropdown @label:"Method" @options:"triMean","truncatedMean","thresholdedMean","median" @description:"The method for calculating the average gene expression per cell group."
//* params.trim =  "0.1"  //* @input @label:"Trim" @description:"The fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed."
//* params.min_cells =  "10"  //* @input @label:"Min. number of cells" @description:"The minimum number of cells required in each cell group for cell-cell communication."
//* params.cell_groups =  "no"  //* @dropdown @options:"yes","no" @label:"Subset cell groups?" @show_settings:"params.sources_use","params.targets_use"
//* params.sources_use =  ""  //* @input @label:"Source cell groups" @description:"Comma separated index or the name of source cell groups" @optional
//* params.targets_use =  ""  //* @input @label:"Target cell groups" @description:"Comma separated index or the name of target cell groups" @optional

threads = task.cpus
source_cells = params.sources_use ?: 'FALSE'
target_cells = params.targets_use ?: 'FALSE'

// species_map = [
//     "human_hg38_gencode_v32_cellranger_v6"      : "human",
//     "human_hg38_cellranger_GRCh38-2024-A"       : "human",
//     "mouse_mm10_gencode_vm23_cellranger_v6"     : "mouse",
//     "mouse_GRCm39_cellranger_GRCm39-2024-A"      : "mouse"
// ]
// organism = species_map.get(params.genome_build, "human")

"""
rds_file="\$(which create_cellchat_obj.R)"

sed -i 's/\\r${'$'}//' "\${rds_file}"

create_cellchat_obj.R ${seurat_obj} ${params.ident} ${params.species} ${params.db_type} ${threads} ${params.smooth} \
    ${params.cell_groups} ${source_cells} ${target_cells} ${params.trim} ${params.min_cells}

rmd_file="\$(which signalling.rmd)"
mkdir -p signalling_scripts

cp "\${rmd_file}" signalling_scripts/
"""

}


process CellChat2_cellchat_comparison {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /scripts\/.*.rmd$/) "scripts/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /report_input_files\/.*.rds$/) "comparison_analysis_result/$filename"}
input:
 path rds

output:
 path "scripts/*.rmd"  ,emit:g580_5_rMarkdown00 
 path "report_input_files/*.rds"  ,emit:g580_5_rdsFile11 


script:

"""
ls
files=( cellchat_list.rds )
echo "\${files[@]}"
SAMPLES=()

rmd_file="\$(which compare.rmd)"
mkdir -p report_input_files
mkdir -p scripts

cp "\${rmd_file}" scripts/
cp "\${files[@]}" report_input_files/

"""
}


workflow {


mkfastq_prep()
g_53_bcl00_g_51 = mkfastq_prep.out.g_53_bcl00_g_51

g_1_1_g_51= g_1_1_g_51.ifEmpty(ch_empty_file_1) 


if (!(params.run_mkfastq == "yes")){
g_1_1_g_51.set{g_51_reads00_g_5}
g_51_outputDir11 = Channel.empty()
g_51_outputHTML22 = Channel.empty()
} else {

mkfastq(g_53_bcl00_g_51.flatMap(), ( params.reads ? g_1_1_g_51 : g_1_1_g_51.first() ) ,g_2_2_g_51)
g_51_reads00_g_5 = mkfastq.out.g_51_reads00_g_5.flatten().map { file -> tmpBase = file.baseName; return tuple(tmpBase.substring(0,tmpBase.indexOf('.')), file) }.groupTuple()
g_51_outputDir11 = mkfastq.out.g_51_outputDir11
g_51_outputHTML22 = mkfastq.out.g_51_outputHTML22
}


Check_and_Build_Module_Check_Genome_GTF()
g17_21_genome00_g17_58 = Check_and_Build_Module_Check_Genome_GTF.out.g17_21_genome00_g17_58
g17_21_gtfFile10_g17_57 = Check_and_Build_Module_Check_Genome_GTF.out.g17_21_gtfFile10_g17_57


if (!(params.replace_geneID_with_geneName == "yes")){
g17_21_gtfFile10_g17_57.set{g17_57_gtfFile01_g17_58}
} else {

Check_and_Build_Module_convert_gtf_attributes(g17_21_gtfFile10_g17_57)
g17_57_gtfFile01_g17_58 = Check_and_Build_Module_convert_gtf_attributes.out.g17_57_gtfFile01_g17_58
}



if (!(params.add_sequences_to_reference == "yes")){
g17_21_genome00_g17_58.set{g17_58_genome00_g17_52}
(g17_58_genome01_g17_54) = [g17_58_genome00_g17_52]
g17_57_gtfFile01_g17_58.set{g17_58_gtfFile10_g17_53}
(g17_58_gtfFile10_g17_54) = [g17_58_gtfFile10_g17_53]
} else {

Check_and_Build_Module_Add_custom_seq_to_genome_gtf(g17_21_genome00_g17_58,g17_57_gtfFile01_g17_58,g_18_2_g17_58,g_48_3_g17_58)
g17_58_genome00_g17_52 = Check_and_Build_Module_Add_custom_seq_to_genome_gtf.out.g17_58_genome00_g17_52
(g17_58_genome01_g17_54) = [g17_58_genome00_g17_52]
g17_58_gtfFile10_g17_53 = Check_and_Build_Module_Add_custom_seq_to_genome_gtf.out.g17_58_gtfFile10_g17_53
(g17_58_gtfFile10_g17_54) = [g17_58_gtfFile10_g17_53]
}


Check_and_Build_Module_Check_BED12(g17_58_gtfFile10_g17_53)
g17_53_bed03_g17_54 = Check_and_Build_Module_Check_BED12.out.g17_53_bed03_g17_54


Check_and_Build_Module_Check_chrom_sizes_and_index(g17_58_genome00_g17_52)
g17_52_genomeSizes02_g17_54 = Check_and_Build_Module_Check_chrom_sizes_and_index.out.g17_52_genomeSizes02_g17_54

g17_58_gtfFile10_g17_54= g17_58_gtfFile10_g17_54.ifEmpty(ch_empty_file_1) 
g17_58_genome01_g17_54= g17_58_genome01_g17_54.ifEmpty(ch_empty_file_2) 
g17_52_genomeSizes02_g17_54= g17_52_genomeSizes02_g17_54.ifEmpty(ch_empty_file_3) 
g17_53_bed03_g17_54= g17_53_bed03_g17_54.ifEmpty(ch_empty_file_4) 


Check_and_Build_Module_check_files(g17_58_gtfFile10_g17_54,g17_58_genome01_g17_54,g17_52_genomeSizes02_g17_54,g17_53_bed03_g17_54)
g17_54_gtfFile01_g_19 = Check_and_Build_Module_check_files.out.g17_54_gtfFile01_g_19
(g17_54_gtfFile03_g57_1) = [g17_54_gtfFile01_g_19]
g17_54_genome10_g_19 = Check_and_Build_Module_check_files.out.g17_54_genome10_g_19
g17_54_genomeSizes22 = Check_and_Build_Module_check_files.out.g17_54_genomeSizes22
g17_54_bed33 = Check_and_Build_Module_check_files.out.g17_54_bed33


cell_ranger_mkref(g17_54_genome10_g_19,g17_54_gtfFile01_g_19)
g_19_reference00_g_20 = cell_ranger_mkref.out.g_19_reference00_g_20

g_19_reference00_g_20= g_19_reference00_g_20.ifEmpty(ch_empty_file_1) 


cellranger_ref_checker(g_19_reference00_g_20)
g_20_reference02_g_5 = cellranger_ref_checker.out.g_20_reference02_g_5


Count(g_51_reads00_g_5,g_2_1_g_5,g_20_reference02_g_5)
g_5_outputDir00_g57_5 = Count.out.g_5_outputDir00_g57_5
g_5_outputHTML11 = Count.out.g_5_outputHTML11
g_5_outputDir22 = Count.out.g_5_outputDir22
g_5_outputDir33 = Count.out.g_5_outputDir33
g_5_h5_file40_g_49 = Count.out.g_5_h5_file40_g_49
g_5_h5_file51_g_49 = Count.out.g_5_h5_file51_g_49


Ambient_RNA_QC(g_5_h5_file40_g_49,g_5_h5_file51_g_49)
g_49_h5_file00_g36_0 = Ambient_RNA_QC.out.g_49_h5_file00_g36_0
g_49_outputFileHTML11 = Ambient_RNA_QC.out.g_49_outputFileHTML11



scRNA_Analysis_Module_Quality_Control_and_Filtering(g_49_h5_file00_g36_0,g_41_1_g36_0)
g36_0_rdsFile00_g36_14 = scRNA_Analysis_Module_Quality_Control_and_Filtering.out.g36_0_rdsFile00_g36_14
g36_0_outputFileHTML11 = scRNA_Analysis_Module_Quality_Control_and_Filtering.out.g36_0_outputFileHTML11
g36_0_outFileTSV20_g36_34 = scRNA_Analysis_Module_Quality_Control_and_Filtering.out.g36_0_outFileTSV20_g36_34


scRNA_Analysis_Module_Merge_Seurat_Objects(g36_0_rdsFile00_g36_14.collect())
g36_14_rdsFile00_g36_17 = scRNA_Analysis_Module_Merge_Seurat_Objects.out.g36_14_rdsFile00_g36_17


scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction(g36_14_rdsFile00_g36_17)
g36_17_rdsFile00_g36_19 = scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.out.g36_17_rdsFile00_g36_19


scRNA_Analysis_Module_Clustering_and_Find_Markers(g36_17_rdsFile00_g36_19)
g36_19_outputHTML00 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g36_19_outputHTML00
g36_19_rdsFile10_g36_36 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g36_19_rdsFile10_g36_36
(g36_19_rdsFile10_g36_40) = [g36_19_rdsFile10_g36_36]
g36_19_outFileTSV22 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g36_19_outFileTSV22


scRNA_Analysis_Module_filter_summary(g36_0_outFileTSV20_g36_34.collect())
g36_34_outputFileHTML00 = scRNA_Analysis_Module_filter_summary.out.g36_34_outputFileHTML00


scRNA_Analysis_Module_sc_annotation(g36_19_rdsFile10_g36_36)
g36_36_rdsFile01_g36_40 = scRNA_Analysis_Module_sc_annotation.out.g36_36_rdsFile01_g36_40
(g36_36_rdsFile00_g36_22,g36_36_rdsFile00_g36_30,g36_36_rdsFile00_g580_1) = [g36_36_rdsFile01_g36_40,g36_36_rdsFile01_g36_40,g36_36_rdsFile01_g36_40]


scRNA_Analysis_Module_SCEtoLOOM(g36_36_rdsFile00_g36_30)
g36_30_outputFileOut00_g_572 = scRNA_Analysis_Module_SCEtoLOOM.out.g36_30_outputFileOut00_g_572


pySCENIC(g36_30_outputFileOut00_g_572,g_576_1_g_572,g_577_2_g_572,g_578_3_g_572)
g_572_zipFile00 = pySCENIC.out.g_572_zipFile00
g_572_loom11 = pySCENIC.out.g_572_loom11


scRNA_Analysis_Module_Create_h5ad(g36_36_rdsFile00_g36_22)
g36_22_h5ad_file01_g57_12 = scRNA_Analysis_Module_Create_h5ad.out.g36_22_h5ad_file01_g57_12

g36_36_rdsFile01_g36_40= g36_36_rdsFile01_g36_40.ifEmpty(ch_empty_file_1) 


scRNA_Analysis_Module_tcr_analysis(g36_19_rdsFile10_g36_40,g36_36_rdsFile01_g36_40)
g36_40_rdsFile00 = scRNA_Analysis_Module_tcr_analysis.out.g36_40_rdsFile00

g_5_outputDir00_g57_5= g_5_outputDir00_g57_5.ifEmpty(ch_empty_file_1) 


RNA_Velocity_Module_prepare_input_velocyto(g_5_outputDir00_g57_5)
g57_5_bamFile00_g57_1 = RNA_Velocity_Module_prepare_input_velocyto.out.g57_5_bamFile00_g57_1
g57_5_inputFileTsv11_g57_1 = RNA_Velocity_Module_prepare_input_velocyto.out.g57_5_inputFileTsv11_g57_1



RNA_Velocity_Module_velocyto(g57_5_bamFile00_g57_1,g57_5_inputFileTsv11_g57_1,g_59_2_g57_1,g17_54_gtfFile03_g57_1)
g57_1_loom00_g57_12 = RNA_Velocity_Module_velocyto.out.g57_1_loom00_g57_12


RNA_Velocity_Module_process_anndata(g57_1_loom00_g57_12.collect(),g36_22_h5ad_file01_g57_12)
g57_12_h5ad00 = RNA_Velocity_Module_process_anndata.out.g57_12_h5ad00


CellChat2_create_cellchat_obj(g36_36_rdsFile00_g580_1)
g580_1_rdsFile00 = CellChat2_create_cellchat_obj.out.g580_1_rdsFile00
g580_1_rMarkdown11 = CellChat2_create_cellchat_obj.out.g580_1_rMarkdown11
g580_1_rdsFile20_g580_5 = CellChat2_create_cellchat_obj.out.g580_1_rdsFile20_g580_5


CellChat2_cellchat_comparison(g580_1_rdsFile20_g580_5.collect())
g580_5_rMarkdown00 = CellChat2_cellchat_comparison.out.g580_5_rMarkdown00
g580_5_rdsFile11 = CellChat2_cellchat_comparison.out.g580_5_rdsFile11


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
