$HOSTNAME = ""
params.outdir = 'results'  

def pathChecker(input, path, type){
	recursive = (type == "folder") ? "--recursive" : ""
	cmd = "mkdir -p check && mv ${input} check/. "
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

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'

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

container 'quay.io/ummsbiocore/custom_sequence_to_genome_gtf:1.0'

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

container "${ params.IMAGE_BASE ? "${params.IMAGE_BASE}/rnaseq:4.0" : "quay.io/ummsbiocore/rnaseq:4.0" }"

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

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'
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

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'
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
 path "${name}_outs"  ,emit:g_5_outputDir00 
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
 tuple val(name), file("${name}_corrected_feature_bc_matrix.h5")  ,emit:g_49_h5_file00_g36_0 
 path "*_Ambient_RNA_QC.html"  ,emit:g_49_outputFileHTML11 

container 'quay.io/ummsbiocore/ambient_rna_removal:1.0'
label 'Ambient_RNA_QC'

shell:
Ambient_RNA_Removal = (params.Ambient_RNA_Removal == "yes") ? "TRUE" : "FALSE"

'''
#!/usr/bin/env perl

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

process scRNA_Analysis_Module_Load_Data_h5 {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_filtering_report.html$/) "QC_Reports/$filename"}
input:
 tuple val(name), file(Input)
 path Metadata

output:
 path "${name}.rds"  ,emit:g36_0_rdsFile00_g36_14 
 tuple val(name),file("${name}_filtering_report.html")  ,emit:g36_0_outputFileHTML11 

label 'scrna_seurat'

when:
(params.run_scRNA_Analysis && (params.run_scRNA_Analysis == "yes")) || !params.run_scRNA_Analysis

shell:
minUMI = params.scRNA_Analysis_Module_Load_Data_h5.minUMI
maxUMI = params.scRNA_Analysis_Module_Load_Data_h5.maxUMI
minFeature = params.scRNA_Analysis_Module_Load_Data_h5.minFeature
maxFeature = params.scRNA_Analysis_Module_Load_Data_h5.maxFeature
percent_mt = params.scRNA_Analysis_Module_Load_Data_h5.percent_mt
percent_ribo = params.scRNA_Analysis_Module_Load_Data_h5.percent_ribo
varFeatures = params.scRNA_Analysis_Module_Load_Data_h5.varFeatures
normMethod = params.scRNA_Analysis_Module_Load_Data_h5.normMethod
DoubletRemoval = params.scRNA_Analysis_Module_Load_Data_h5.DoubletRemoval
RemoveMitoGenes = params.scRNA_Analysis_Module_Load_Data_h5.RemoveMitoGenes
RemoveRiboGenes = params.scRNA_Analysis_Module_Load_Data_h5.RemoveRiboGenes

'''
#!/usr/bin/env perl

my $script = <<'EOF';
---
title: "Single Cell RNA-Seq QC and Filtering Report"
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
  SampleName: ''
  SamplePath: ''
  RawInput: FALSE
  nFeature_RNA_lower_threshold: 0.01
  nFeature_RNA_higher_threshold: 0.99
  nCount_RNA_lower_threshold: 0.01
  nCount_RNA_higher_threshold: 0.99
  ribosomal_contents_threshold: 25
  mitochondrial_contents_threshold: 50
  doubletpercentage: 0.01
  filtered_output: ""
  normMethod: "LogNormalize"
  varFeatures: 3000
  DoubletRemoval: TRUE
  Metadata: ''
  RemoveMitoGenes: FALSE
  RemoveRiboGenes: FALSE

---


```{r setup, include=FALSE}

suppressPackageStartupMessages({
library(dplyr)
library(Seurat)
library(ggplot2)
library(DropletUtils)
library(DoubletFinder)
library(clustree)

})

```

```{r helper_functions, include=FALSE}

get_present_pseudogenes <- function(seurat_object, pseudogenes) {

  # Ensure that the Seurat object has rownames (gene names)

  if (is.null(rownames(seurat_object))) {

    stop("The Seurat object does not have gene names as rownames.")

  }

  

  gene_names <- rownames(seurat_object)

  pseudogenes_in_dataset <- intersect(pseudogenes, gene_names)

  

  return(pseudogenes_in_dataset)

}


```

# Read in samples

If the data is a raw Count Matrix, meaning that the empty droplets are not filtered out by cellranger pipeline or Drop-Seq pipeline, the emptyDrops algorithm from DropUtils will be run in order to remove the empty droplets.

```{r read in samples}
Data=Read10X_h5(params$SamplePath)

MultiModal=F
if (is.list(Data)) {
  MultiModal=T
  Raw=Data
  Data=Data[["Gene Expression"]]
  
}
#If any cells with 0 in the dataset, raw input is true
if (any(colSums(Data)==0)) {
	RawInput=TRUE
} else {
	RawInput=FALSE
}


if (params$RemoveRiboGenes) {
	
	
if (any(grepl("^RP[SL]",rownames(Data)))) {
	rb.genes <- rownames(Data)[grep("^RP[SL]",rownames(Data))]
	Data=Data[!rownames(Data)%in%rb.genes,]
}
else if (any(grepl("^Rp[sl]",rownames(Data)))) {
	rb.genes <- rownames(Data)[grep("^Rp[sl]",rownames(Data))]
	GTgenes=c("Gm42418","AY036118")
	rb.genes <-c(rb.genes,GTgenes)
	Data=Data[!rownames(Data)%in%rb.genes,]
}
}

mt.genes=c()
if (params$RemoveMitoGenes) {


if (any(grepl("^MT-",rownames(Data)))) {
	mt.genes <- rownames(Data)[grep("^MT-",rownames(Data))]
	Data=Data[!rownames(Data)%in%mt.genes,]
}
else if (any(grepl("^mt-",rownames(Data)))) {
	mt.genes <- rownames(Data)[grep("^mt-",rownames(Data))]
	Data=Data[!rownames(Data)%in%mt.genes,]
}
}

if (RawInput) {
  empty=emptyDrops(Data[!rownames(Data)%in%mt.genes,])
  empty=data.frame(empty)
  empty=empty[!is.na(empty$FDR),]
  empty$DropIdentity=ifelse(empty$FDR<0.001,yes="Non Empty",
                            no="Empty")
  ggplot(empty,aes(x=DropIdentity,y=Total))+geom_bar(stat = 'identity')+xlab("Droplet classification")+ylab("Number of UMIs per cell")+ggtitle("Empty Droplet classification")
  Data=Data[,rownames(empty)[empty$FDR<0.05]]

}


```

# QC {.tabset}

In this section the number of genes, number of UMIs and the percentage of mitochondrial contents and ribosomal contents will be visualized.

Mitochondrial contents and ribosomal contents will be calculated based on the mitochondrial and ribosomal genes.

Cells with very high mitochondrial and ribosomal contents will bias the downstream clustering and differential expression analysis.

```{r QC, fig.width=5,fig.height=5}

Data=CreateSeuratObject(Data)
Data$sample=params$SampleName
Data$orig.ident=params$SampleName

if (any(grepl("^MT-",rownames(Data)))|any(grepl("^mt-",rownames(Data)))) {
	if (any(grepl("^MT-",rownames(Data)))) {
	Data[["percent.mt"]]=PercentageFeatureSet(Data,pattern="^MT-")
	rb.genes <- rownames(Data)[grep("^RP[SL]",rownames(Data))]
	Data[["percent.ribo"]] <- PercentageFeatureSet(Data, features = rb.genes)	


}
else if (any(grepl("^mt-",rownames(Data)))) {
	Data[["percent.mt"]]=PercentageFeatureSet(Data,pattern="^mt-")
	
	rb.genes <- rownames(Data)[grep("^Rp[sl]",rownames(Data))]
	present_pseudogenes <- get_present_pseudogenes(Data, c("Gm42418","AY036118") )
    rb.genes <- c( rb.genes, present_pseudogenes )
	
	Data[["percent.ribo"]] <- PercentageFeatureSet(Data, features = rb.genes)	

}} else {
	Data[["percent.mt"]]=0
		Data[["percent.ribo"]]=0

}

Unfiltered = Data

```

## Violin plot of number of gene per cell

``` {r vnFeature, fig.width=5,fig.height=5}
VlnPlot(Unfiltered,features = c("nFeature_RNA"),pt.size = 0)
```

## Violin plot of number of UMIs per cell

``` {r vnCount, fig.width=5,fig.height=5}
VlnPlot(Unfiltered,features = c("nCount_RNA"),pt.size = 0)
```

## Violin plot of mitochondrial contents per cell

``` {r vmt, fig.width=5,fig.height=5}
VlnPlot(Unfiltered,features = c("percent.mt"),pt.size = 0)
```

## Violin plot of ribosomal contents per cell

``` {r vrb, fig.width=5,fig.height=5}
VlnPlot(Unfiltered,features = c("percent.ribo"),pt.size = 0)
```

## Scatter plot of number of genes per cell vs number of UMIs per cell

``` {r sfu, fig.width=5,fig.height=5}
FeatureScatter(Unfiltered,
               feature1 = "nFeature_RNA",
               feature2 = "nCount_RNA")+ NoLegend()
```

## Scatter plot of mitochondrial contents

``` {r smfu, fig.width=5,fig.height=5}
if (isFALSE(all(Unfiltered[["percent.mt"]]==0))){
  FeatureScatter(Unfiltered,
               feature1 = "nFeature_RNA",
               feature2 = "percent.mt")+ NoLegend()

  FeatureScatter(Unfiltered,
               feature1 = "nCount_RNA",
               feature2 = "percent.mt")+ NoLegend()
  }

```

## Scatter plot of mitochondrial contents

``` {r srfu, fig.width=5,fig.height=5}
if (isFALSE(all(Unfiltered[["percent.ribo"]]==0))){
  FeatureScatter(Unfiltered,
               feature1 = "nFeature_RNA",
               feature2 = "percent.ribo")+ NoLegend()

  FeatureScatter(Unfiltered,
               feature1 = "nCount_RNA",
               feature2 = "percent.ribo")+ NoLegend()
}

```

# Doublet Classification {.tabset}

Doublets, or sometimes called multiplets, are the droplets which include two or more cells. Including these droplets in the downstream analysis will bias the results because these droplets include gene expression profiles of more than 1 cell.

DoubletFinder is used to classify the doublet. The violin plot in this section will show the number features and UMIs of the doublets vs that of non-doublets.

## Doublet Simulation and detection 
```{r Doublet Removal, error=FALSE, fig.height=5, fig.width=5, message=FALSE, warning=FALSE,results = FALSE}

DoubletRemovalHandle=as.logical(params$DoubletRemoval)

Data$Doublet.Classification="SingleLet"

if (is.na(DoubletRemovalHandle)) {
	if (MultiModal) {
		if (any(grepl("Multiplexing",names(Raw)))) {
		
			DoubletRemovalHandle=FALSE
			
		} else {
		
			DoubletRemovalHandle=TRUE
			
		}
	}
	else {
	DoubletRemovalHandle=TRUE

	}
} 

if (DoubletRemovalHandle) {
DoubletRemoval <- NormalizeData(Data)
DoubletRemoval <- FindVariableFeatures(DoubletRemoval)
DoubletRemoval <- ScaleData(DoubletRemoval)

DoubletRemoval <- RunPCA(DoubletRemoval,npcs=100)

pc.changes=diff(diff(DoubletRemoval@reductions$pca@stdev))
pc.changes=abs(pc.changes)
pc.changes=which(pc.changes>=mean(pc.changes))

DoubletRemoval=FindNeighbors(DoubletRemoval,dims = 1:(max(pc.changes)+2),reduction = "pca")

DoubletRemoval=FindClusters(DoubletRemoval,resolution = seq(2.0,0.1,-0.1))

names=paste0(DefaultAssay(DoubletRemoval),"_snn_res.")
SC3_Stability=clustree(DoubletRemoval,prefix = names)
SC3_Stability.results=SC3_Stability$data
SC3_Stability.results=SC3_Stability.results[,c(names,"sc3_stability")]
colnames(SC3_Stability.results)[1]="resolution"
SC3_Stability.results.mean=aggregate(sc3_stability~resolution,SC3_Stability.results,mean)
colnames(SC3_Stability.results.mean)[2]="sc3_stability_mean"
Idents(DoubletRemoval)=paste0(DefaultAssay(DoubletRemoval),"_snn_res.",max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean==max(SC3_Stability.results.mean$sc3_stability_mean)]))

DoubletRemoval$seurat_clusters=Idents(DoubletRemoval)



params_selection <- paramSweep_v3(DoubletRemoval, PCs = 1:(max(pc.changes)+2), sct = F)
params_selection <- summarizeSweep(params_selection, GT = FALSE)
params_selection <- find.pK(params_selection)



annotations <- DoubletRemoval@meta.data$DoubletRemoval$seurat_clusters

homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(params$doubletpercentage*nrow(DoubletRemoval@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


DoubletRemoval <- doubletFinder_v3(DoubletRemoval, PCs = 1:(max(pc.changes)+2), pN = 0.25, pK =as.numeric(as.character(params_selection$pK))[params_selection$BCmetric==max(params_selection$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
doublet.classification.name=colnames(DoubletRemoval@meta.data)[ncol(DoubletRemoval@meta.data)]
DoubletRemoval$Doublet.Classification=DoubletRemoval@meta.data[,doublet.classification.name]

Data$Doublet.Classification=DoubletRemoval$Doublet.Classification
}

```

## Doublet Visualization

```{r Doublet Visualization, error=FALSE, fig.height=5, fig.width=5, message=FALSE, warning=FALSE,results = FALSE }

if (DoubletRemovalHandle) {
VlnPlot(DoubletRemoval,c("nCount_RNA","nFeature_RNA"),group.by = "Doublet.Classification",pt.size = 0)+NoLegend()
}

```

# Filtering {.tabset}

Based on the quantile information, cells with too low or too high number of features and UMIs are filtered out, which can be observed on the violin plot.

## Violin Plots of number of genes and UMIs of doublets and singlets

```{r Filtering, fig.width=5,fig.height=5,error=FALSE,warning=FALSE,message=FALSE}

if (DoubletRemovalHandle){
Data$Filtering = ifelse(
  Data$Doublet.Classification!='Doublet'&
  Data$nFeature_RNA>quantile(Data$nFeature_RNA,params$nFeature_RNA_lower_threshold)&
  Data$nFeature_RNA<quantile(Data$nFeature_RNA,params$nFeature_RNA_higher_threshold)&
  Data$nCount_RNA>quantile(Data$nCount_RNA,params$nCount_RNA_lower_threshold)&
  Data$nCount_RNA<quantile(Data$nCount_RNA,params$nCount_RNA_higher_threshold)&
  Data$percent.mt<params$mitochondrial_contents_threshold&
  Data$percent.ribo<params$ribosomal_contents_threshold,
  yes="Keep",
  no="Drop"
)

} else {
Data$Filtering = ifelse(
  Data$nFeature_RNA>quantile(Data$nFeature_RNA,params$nFeature_RNA_lower_threshold)&
  Data$nFeature_RNA<quantile(Data$nFeature_RNA,params$nFeature_RNA_higher_threshold)&
  Data$nCount_RNA>quantile(Data$nCount_RNA,params$nCount_RNA_lower_threshold)&
  Data$nCount_RNA<quantile(Data$nCount_RNA,params$nCount_RNA_higher_threshold)&
  Data$percent.mt<params$mitochondrial_contents_threshold&
  Data$percent.ribo<params$ribosomal_contents_threshold,
  yes="Keep",
  no="Drop")
}

VlnPlot(Data,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo",pt.size = 0),group.by = "Filtering",ncol = 2,pt.size = 0)+NoLegend()

Data=subset(Data,subset=Filtering=="Keep")

```

## Barplot of filtering results

```{r Filtering barplot, fig.width=5,fig.height=5,error=FALSE,warning=FALSE,message=FALSE}

Filtering_statistics=data.frame(Category=c("Unfiltered","Filtered"),CellNumber=c(ncol(Unfiltered),ncol(Data)))
Filtering_statistics$Category=factor(Filtering_statistics$Category,levels = c("Unfiltered","Filtered"))

ggplot(Filtering_statistics,aes(x=Category,y=CellNumber,label=CellNumber))+geom_bar(stat="identity")+geom_text(size = 3, position = position_stack(vjust = 0.75))

```


# Normalization

As the violin plots shown in the QC section, the sequencing depth and coverage of each cell in a single cell RNA-Seq dataset vary significantly.

The normalization step normalize the gene expression profile of each cell, which makes them comparable to each other in the downstream analysis.

The SCTransform is recommended as it enhances the biological signature in the data, however it is quite time-consuming and memory-consuming. 

The LogNormalize is very standard practice time-efficient.

```{r Normalization, error=FALSE, fig.height=5, fig.width=5, message=FALSE, warning=FALSE,results = FALSE}

tryCatch({if (file.exists(params$Metadata)) {
	Metadata=read.table(params$Metadata,sep="\t",check.names=F,header=T,row.names=NULL)
	if ("Sample"%in%colnames(Metadata)) {
	AttributeList=colnames(Metadata)[colnames(Metadata)!="Sample"]
	if (params$SampleName%in%Metadata$Sample) {
		for (i in 1:length(AttributeList)) {
			Data[[AttributeList[i]]]=as.character(Metadata[Metadata$Sample==params$SampleName,AttributeList[i]])
		}
	} else {
		for (i in 1:length(AttributeList)) {
			Data[[AttributeList[i]]]=""
		}
	}

}
}
},
error=function(err) {
	print("No metadata information is added to dataset")
})

if (params$normMethod=="SCT") {
	if (all(Data[["percent.mt"]]==0)) {
		Data <- SCTransform(Data,variable.features.n=params$varFeatures)
	} else {
		Data <- SCTransform(Data,variable.features.n=params$varFeatures,vars.to.regress="percent.mt")

	}
	} else {
	Data <- NormalizeData(object = Data,normalization.method=params$normMethod)
}

if (MultiModal) {
  for (name in names(Raw)[names(Raw)!="Gene Expression"]) {
    Data[[gsub(" ","",name)]]=CreateAssayObject((Raw[[name]][,colnames(Data)]))
  }

}

saveRDS(Data,params$filtered_output)


```

EOF

open OUT, ">!{name}_filtering_rmark.rmd";
print OUT $script;
close OUT;

runCommand("Rscript -e 'rmarkdown::render(\\"!{name}_filtering_rmark.rmd\\",\\"html_document\\", output_file = \\"!{name}_filtering_report.html\\",
params = list(SampleName=\\"!{name}\\",
SamplePath=\\"!{Input}\\",
nFeature_RNA_lower_threshold=as.numeric(!{minFeature}),
nFeature_RNA_higher_threshold=as.numeric(!{maxFeature}),
nCount_RNA_lower_threshold=as.numeric(!{minUMI}),
nCount_RNA_higher_threshold=as.numeric(!{maxUMI}),
ribosomal_contents_threshold=as.numeric(!{percent_ribo}),
mitochondrial_contents_threshold=as.numeric(!{percent_mt}),
filtered_output=\\"!{name}.rds\\",
normMethod=\\"!{normMethod}\\",
varFeatures=as.numeric(!{varFeatures}),
DoubletRemoval=as.character(\\"!{DoubletRemoval}\\"),
Metadata=\\"!{Metadata}\\",
RemoveMitoGenes=as.logical(\\"!{RemoveMitoGenes}\\"),
RemoveRiboGenes=as.logical(\\"!{RemoveRiboGenes}\\")
))'");

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

label 'scrna_seurat'

shell:

'''
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

'''


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

label 'scrna_seurat'


shell:

varFeatures = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.varFeatures

selmethod = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.selmethod

Batch_Effect_Correction = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.Batch_Effect_Correction

WNN = params.scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.WNN

'''
#!/usr/bin/env Rscript

# libraries
library(Seurat)
library(dplyr)
#install.packages("harmony",repos = "http://cran.us.r-project.org")
library(harmony)

selmethod <- "!{selmethod}"
varFeatures <- "!{varFeatures}"

Data=readRDS("!{seurat_object}")
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
		if (as.logical("!{Batch_Effect_Correction}")){
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
		if (as.logical("!{Batch_Effect_Correction}")){
		Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
		}

	}
}

if ("!{WNN}"!="") {
original.assay=DefaultAssay(Data)

DefaultAssay(Data)="!{WNN}"

Data=NormalizeData(Data,normalization.method = "CLR",margin=2)

VariableFeatures(Data)=rownames(Data)

Data=ScaleData(Data)

Data=RunPCA(Data,reduction.name = "wpca")

if (Multi_sample==1) {
	Data=RunHarmony(Data,group.by.vars = "sample",assay.use = "!{WNN}",reduction = "wpca",reduction.save = "wharmony")

}

DefaultAssay(Data)=original.assay

}




#if (DefaultAssay(Data)=="SCT"){
#
#if (length(unique(Data$sample))==1) {
#	Data=RunPCA(Data,npcs=100)
#} else {
#	Data <- SplitObject(Data, split.by = "sample")
#	variable.features=SelectIntegrationFeatures(object.list = Data, nfeatures = as.numeric(varFeatures))
#	Data <- merge(Data[[1]],Data[-1])
#	VariableFeatures(Data) <- variable
#	if (all(Data[["percent.mt"]]==0)) {
#		Data <- ScaleData(Data)
#	} else {
#		Data <- ScaleData(Data,vars.to.regress="percent.mt")
#	}
#
#	Data=RunPCA(Data,npcs=100)
#	Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
#}
#} else {
#	Data <- FindVariableFeatures(Data,selection.method=selmethod,nfeatures=as.numeric(varFeatures))
#	if (all(Data[["percent.mt"]]==0)) {
#		Data <- ScaleData(Data)
#	} else {
#		Data <- ScaleData(Data,vars.to.regress="percent.mt")
#	}
#	Data=RunPCA(Data,npcs=100)
#	if (length(unique(Data$sample))>1) {
#		Data=RunHarmony(Data,assay.use = DefaultAssay(Data),group.by.vars = "sample",max.iter.harmony = 10000,max.iter.cluster = 10000)
#	}
#}
saveRDS(Data,"Reduced_and_Corrected.rds")

'''


}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 140
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_Clustering_and_Find_Markers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /Final_Report.html$/) "Final_Report/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /Final_Analysis.rds$/) "scViewer/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "ClusterMarkers/$filename"}
input:
 path seurat_object

output:
 path "Final_Report.html"  ,emit:g36_19_outputHTML00 
 path "Final_Analysis.rds"  ,emit:g36_19_rdsFile10_g36_22 
 path "*.tsv"  ,emit:g36_19_outFileTSV22 

label 'scrna_seurat'

shell:
minRes = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.minRes
maxRes = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.maxRes
npcs = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.npcs
runCellFindR = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.runCellFindR
findClusterforallResolution = params.scRNA_Analysis_Module_Clustering_and_Find_Markers.findClusterforallResolution

'''
#!/usr/bin/env perl

my $script = <<'EOF';

---
title: "Single Cell RNA-Seq Clustering Report"
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
  
  SamplePath: ""
  npcs: 25
  minRes: 0.1
  maxRes: 2.0
  filtered_output: ""
  algorithm: 2
---


```{r setup, include=FALSE}

suppressPackageStartupMessages({
library(dplyr)
library(Seurat)
library(ggplot2)
#if(!require(remotes)) install.packages("remotes",repos = "http://cran.us.r-project.org")
#if(!require(data.table)) install.packages("remotes",repos = "http://cran.us.r-project.org")
#if(!require(DT)) install.packages("DT",repos = "http://cran.us.r-project.org")
#if(!require(clustree)) install.packages("clustree",repos = "http://cran.us.r-project.org")

#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)
library(data.table)
library(DT)
library(clustree)

is_cluster <- function(tenx, thresh_genes = 10, thresh_val = log2(2), pval = 1e-4,Idents=NULL){
  Idents(tenx)=Idents
  val = 0 # groups that does not satisfy threshold genes
  counter = 0 # groups that satisfy threshold genes 
  matrix_output <- data.frame(row.names = row.names(tenx))
  tenx=NormalizeData(tenx,assay="RNA")
  for (j in sort(unique(Idents(tenx)))){
    if (sum(Idents(tenx) == j) < 5){
      return(FALSE)
    }
    markers <- FindMarkers(tenx,ident.1=j,group.by=Idents,logfc.threshold=thresh_val,assay="RNA")
    markers <- markers[markers$p_val_adj < pval,]
    
    #find if the 10th biggest is less than log2, sum 
    print(sort(markers$avg_log2FC, decreasing = TRUE)[thresh_genes])
    # if less than 10 significant genes
    
    if (length((markers$avg_log2FC)) < 10){
      val <- val + 1
    } else if (sort(markers$avg_log2FC, decreasing = TRUE)[thresh_genes] < thresh_val){
      #print(val)
      val <- val + 1
    } else{
      counter = counter + 1
    }
    if (val > 1){
      return(FALSE)
    }
  }
  if (val > 1){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

sub_sampling <- function(object,column="seurat_clusters",subsampling_cells=500) {
	subset=c()
    for (i in unique(object@meta.data[,column])) {
      if(sum(object@meta.data[,column]==i)>=subsampling_cells) {
        cells=sample(rownames(object@meta.data)[object@meta.data[,column]==i],size = subsampling_cells,replace = F)
        subset=c(subset,cells)
      } else {
        cells=rownames(object@meta.data)[object@meta.data[,column]==i]
        subset=c(subset,cells)
        
      }
    }
    return(subset(object,cells=subset))

}



find_res <- function(tenx, initial_res = 0.1, jump = 0.1,max=5.0, thresh_genes = 10, thresh_val = log2(2),subsampling_cells=1000) {
  RES_POST <- initial_res # keeping
  RES_IT <- initial_res # iterative
  Keep_loop=T
  while(Keep_loop){
    # also check if theres only 1 cluster/ then can go up higher es
    # Find number of clusters
    temp <- FindClusters(tenx, resolution = RES_IT,algorithm = 2)
    Idents(temp)="seurat_clusters"
    subset=c()
    for (i in unique(temp$seurat_clusters)) {
      if(sum(temp@meta.data$seurat_clusters==i)>=subsampling_cells) {
        cells=sample(rownames(temp@meta.data)[temp@meta.data$seurat_clusters==i],size = subsampling_cells,replace = F)
        subset=c(subset,cells)
      } else {
        cells=rownames(temp@meta.data)[temp@meta.data$seurat_clusters==i]
        subset=c(subset,cells)
        
      }
    }
    temp=subset(temp,cells=subset)
    length_group <- length(unique(Idents(tenx)))

    if (length_group == 1){
      # still not groups at 0.7 res stop and just step as 1
      if (RES_IT == max){
        Keep_loop=F
      }
    } else{
      if (RES_IT==max) {
        Keep_loop=F
      }
      
      
      testing <- is_cluster(temp,thresh_genes = thresh_genes,Idents='seurat_clusters')
      if (testing == FALSE){ # if not real group
        RES_IT <- RES_IT - jump
        RES_POST <- RES_IT
        print(RES_POST)
        Keep_loop=F
      } else{ # valid groups
        RES_POST <- RES_IT
        RES_IT <- RES_IT + jump
      }
    }
  }
  # if there is only 1 group, return 0,
  return(RES_POST)
}


})
library(cluster)
Data=readRDS(params$SamplePath)
```

# Sample Statistics {.tabset}

In this section violin plots will show the number of genes, UMIs, mitochondrial percentages and ribosomal percentages of cells after filtering.

If there is more than 1 sample in the dataset, the violin plots will group the cells by samples. This is a very good way to compare quality of cells.

## Number of genes per cell

```{r Number of genes per cell by sample, fig.width=10,fig.height =10}

VlnPlot(Data,"nFeature_RNA",group.by="sample",pt.size = 0)+NoLegend()

```

## Number of UMIs per cell

```{r Number of UMIs per cell by sample, fig.width=10,fig.height =10}

VlnPlot(Data,"nCount_RNA",group.by="sample",pt.size = 0)+NoLegend()

```

## mitochondrial percentage per cell

```{r mitochondrial percentage per cell by sample, fig.width=10,fig.height =10}

VlnPlot(Data,"percent.mt",group.by="sample",pt.size = 0)+NoLegend()

```

## ribosomal percentage per cell

```{r ribosomal percentage per cell by sample, fig.width=10,fig.height =10}

VlnPlot(Data,"percent.ribo",group.by="sample",pt.size = 0)+NoLegend()

```

# Visualize PCA results

The results of principle component analysis give more insight than people usually realize. For example, the genes that contribute the most to the top principle components can help people to do sanity check of the data: ideally these genes will match the gene markers of the cell sub-populations in the data. This means that the cell heterogeneity is being captured.

```{r Dimension Reduction heatmap, fig.width=10,fig.height =10}

DimHeatmap(Data,dims = 1:10,nfeatures = 9,balanced = T,cells = 500,ncol = 3)

```

The genes that contribute the most to the top principle components can be visualized using Dimention Reduction heatmap shown above.


Another good way to visualize the PCA results is elbow plot, which plot the standard deviation of cells on each principle components. 

```{r Elbow Plot, fig.width=10,fig.height =10}

ElbowpotData=data.frame(stdev=Data@reductions[[ifelse("harmony"%in%names(Data@reductions),yes="harmony",no="pca")]]@stdev,PCs=seq(1,length(Data@reductions[[ifelse("harmony"%in%names(Data@reductions),yes="harmony",no="pca")]]@stdev)))

pc.changes=(diff((ElbowpotData$stdev)))*(-1)
pc.changes.raw=pc.changes
pc.changes=which(pc.changes>=mean(pc.changes.raw[pc.changes.raw>0]))
ggplot(ElbowpotData,aes(x=PCs,y=stdev,label=PCs))+geom_point()+theme_bw()+geom_vline(xintercept = max(pc.changes)+2,color="darkred")+geom_vline(xintercept = params$npcs,color="green")


```

# Visualize the UMAP and TSNE {.tabset}

Before doing any clustering, let us first use tSNE and UMAP to see if the batch effects between samples are removed.

If you only have one sample in this analysis, then please just enjoy the beautiful tSNE and UMAP. 

## tSNE {.tabset}

```{r tSNE, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

Data=RunTSNE(Data, check_duplicates = FALSE, dims = 1:ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs)
,reduction = ifelse("harmony"%in%names(Data@reductions),
                   yes="harmony",
                   no="pca"))
```

### tSNE colored by sample

```{r tSNE colored by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "tsne")+NoLegend()
} else {
  DimPlot(Data,reduction = "tsne",group.by="sample")+NoLegend()

}
```

### tSNE colored by sample without batch effect correction

```{r tSNE without batch effect correction, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

  Temp=Data
  Temp=RunTSNE(Temp,check_duplicates = FALSE, dims = 1:ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs),reduction = "pca")
  DimPlot(Temp,reduction = "tsne",group.by="sample")+NoLegend()

```

### tSNE split by sample

```{r tSNE split by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "tsne")+NoLegend()
} else {
  DimPlot(Data,reduction = "tsne",split.by="sample",ncol=2)+NoLegend()
}
```

### tSNE split by sample without batch effect correction

```{r tSNE split by sample without batch effect correction, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

  DimPlot(Temp,reduction = "tsne",split.by="sample",ncol=2)+NoLegend()
  rm(Temp)

```

## UMAP {.tabset}

```{r umap, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

Data=RunUMAP(Data,dims = 1:ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs),reduction = ifelse("harmony"%in%names(Data@reductions),
                   yes="harmony",
                   no="pca"))
```

### UMAP colored by sample

```{r UMAP colored by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "umap")+NoLegend()
} else {
  DimPlot(Data,reduction = "umap",group.by="sample")+NoLegend()

}
```

### UMAP colored by sample without batch effect correction

```{r UMAP without batch effect correction, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}
  Temp=Data
  Temp=RunUMAP(Temp,dims = 1:(ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs)),reduction = "pca")
  DimPlot(Temp,reduction = "umap",group.by="sample")+NoLegend()
```

### UMAP split by sample

```{r UMAP split by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "umap")+NoLegend()
} else {
  DimPlot(Data,reduction = "umap",split.by="sample",ncol=2)+NoLegend()
}
```

### UMAP split by sample without batch effect correction

```{r UMAP split by sample without batch effect correction, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

  DimPlot(Temp,reduction = "umap",split.by="sample",ncol=2)+NoLegend()
  rm(Temp)
  
```



## Note on tSNE and UMAP.

Both visualizations are widely used. Feel free to choose the one you like.

One thing to take in mind is that it is faster to generate UMAP reduction than to generate tSNE reduction.

# Clustering

## Building the nearest neighborhood graph

In order to cluster the cells, the shared nearest neighborhood graph of cells are constructed using the top principle components (default is 25).


```{r Build snn graph, error=FALSE, fig.height =10, fig.width=10, message=FALSE, warning=FALSE,results = FALSE}
if ("wpca"%in%names(Data@reductions)) {
	ElbowpotData=data.frame(stdev=Data@reductions[[ifelse("wharmony"%in%names(Data@reductions),yes="wharmony",no="wpca")]]@stdev,PCs=seq(1,length(Data@reductions[[ifelse("wharmony"%in%names(Data@reductions),yes="wharmony",no="wpca")]]@stdev)))
	wpc.changes=(diff((ElbowpotData$stdev)))*(-1)
	wpc.changes.raw=wpc.changes
	wpc.changes=which(wpc.changes>=mean(wpc.changes.raw[wpc.changes.raw>0]))
	Data <- FindMultiModalNeighbors(Data,
	reduction.list = list(ifelse("harmony"%in%names(Data@reductions),yes="harmony",no="pca"), 
  ifelse("wharmony"%in%names(Data@reductions),yes="wharmony",no="wpca")
  ), 
  dims.list = list(1:max(pc.changes+1), 1:max(wpc.changes+1)), modality.weight.name = "multi_modal_weight"
)
Data <- RunUMAP(Data, nn.name = "weighted.nn")
Data <- RunTSNE(Data, check_duplicates = FALSE, nn.name = "weighted.nn")

} else {
	Data=FindNeighbors(Data,dims = 1:ifelse(params$npcs==0,yes=max(pc.changes)+1,no=params$npcs),reduction = ifelse("harmony"%in%names(Data@reductions),
                   yes="harmony",
                   no="pca"))

}



```


## Clustering {.tabset}

And then Graph Based Community Detection Algorithm is used to cluster the cells.

In order to select the best clustering resolution, the sc3 stability index is calculate for each resolution. The resolution with the highest mean sc3 stability index (marked by red line in the figure below).

```{r Clustering, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}

if ("wsnn"%in%names(Data@graphs)) {
	Data <- FindClusters(Data, graph.name = "wsnn",algorithm = params$algorithm,resolution = seq(params$maxRes,params$minRes,-0.1))

} else {
	Data=FindClusters(Data,algorithm = params$algorithm,resolution = seq(params$maxRes,params$minRes,-0.1))

}


if ("wsnn"%in%names(Data@graphs)) {
names=paste0("wsnn","_res.")

} else {
names=paste0(DefaultAssay(Data),"_snn_res.")

}

SC3_Stability=clustree(Data,prefix = names)
SC3_Stability.results=SC3_Stability$data
SC3_Stability.results=SC3_Stability.results[,c(names,"sc3_stability")]
colnames(SC3_Stability.results)[1]="resolution"
SC3_Stability.results.mean=aggregate(sc3_stability~resolution,SC3_Stability.results,mean)
colnames(SC3_Stability.results.mean)[2]="sc3_stability_mean"
if ("wsnn"%in%names(Data@graphs)) {
Idents(Data)=paste0("wsnn","_res.",max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean==max(SC3_Stability.results.mean$sc3_stability_mean)]))

} else {
Idents(Data)=paste0(DefaultAssay(Data),"_snn_res.",max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean==max(SC3_Stability.results.mean$sc3_stability_mean)]))
}
Data$seurat_clusters=Idents(Data)

Cluster.distribution=data.frame(table(Data$seurat_clusters,Data$sample))
colnames(Cluster.distribution)=c("Cluster","Sample","CellNumber")



```

### Cluster stability assessment

```{r Cluster stability assessment, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}
ggplot(SC3_Stability.results,aes(x=resolution,y=sc3_stability))+geom_boxplot()+geom_line(data=SC3_Stability.results.mean,aes(x=resolution,y=sc3_stability_mean,group=1))+geom_vline(xintercept = SC3_Stability.results.mean$resolution[SC3_Stability.results.mean$sc3_stability_mean==max(SC3_Stability.results.mean$sc3_stability_mean)],color='red')
```

### Clustree assessment

This figure shows how the cells are assigned as the resolution changes. The color of the arrow shows the amount of cells going into the cluster in the next level and the direction of the arrow shows the identity of cluster that the cells are going to.

As the resolution increases, the arrows will start to appear "messy". This means that the clustering algorithm is having trouble assigning cells.

```{r Clustree assessment, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}
SC3_Stability
```

### CellFindR assessment

This is an algorithm developed by Kevin Yu in the Tward lab at UCSF. The algorithm tries to find the optimal clustering results from the single cell RNA-Seq data.

```{r CellFindR, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}

if (as.logical("!{runCellFindR}")) {
	if ("wsnn"%in%names(Data@graphs)) {
	#CellFindR.resolution=find_res(Data)

	#CellFindR=FindClusters(Data,resolution=CellFindR.resolution,graph.name="wsnn")
	CellFindR=Data
} else {
	CellFindR.resolution=find_res(Data)

	CellFindR=FindClusters(Data,resolution=CellFindR.resolution)

}

Data$CellFindR.clustering.results=CellFindR$seurat_clusters

}


```








### Sample Distribution over Cluster
```{r Sample distribution, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}



ggplot(Cluster.distribution,aes(y=Cluster,x=CellNumber,fill=Sample))+geom_bar(stat="identity",position="fill")

```

### Cluster Distribution over Sample
```{r Cluster distribution, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}

Cluster.distribution=data.frame(table(Data$seurat_clusters,Data$sample))
colnames(Cluster.distribution)=c("Cluster","Sample","CellNumber")

ggplot(Cluster.distribution,aes(y=Sample,x=CellNumber,fill=Cluster))+geom_bar(stat="identity",position="fill")

```

### tSNE Visualization of the Cluster

```{r tSNE Visualization of the Cluster, fig.width=10,fig.height=10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}
DimPlot(Data,reduction = "tsne",label = T)

```

### tSNE distributed by sample

```{r tSNE distributed by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "tsne")+NoLegend()
} else {
  DimPlot(Data,reduction = "tsne",split.by="sample",ncol=2)+NoLegend()
}
```

### UMAP Visualization of the Cluster

```{r UMAP Visualization of the Cluster, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,results=FALSE,echo=FALSE}
DimPlot(Data,reduction = "umap",label = T)

```

### UMAP distributed by sample

```{r UMAP distributed by sample, fig.width=10,fig.height =10,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

if (length(unique(Data$sample))==1){
DimPlot(Data,reduction = "umap")+NoLegend()
} else {
  DimPlot(Data,reduction = "umap",split.by="sample",ncol=2)+NoLegend()
}
```



## Cluster Markers {.tabset}

Use differential expression analysis to find Gene markers for each cluster. These gene markers are very helpful in identifying Cell types.

```{r Find Cluster markers, fig.width=15,fig.height=20,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}

Data=NormalizeData(Data,assay = "RNA")


Data.markers=FindAllMarkers(Data,only.pos = T,assay = "RNA")
write.table(Data.markers,"Cluster.Markers.tsv",quote=F,sep="\t")
if (as.logical("!{findClusterforallResolution}")) {
	for (i in colnames(Data@meta.data)[grepl("snn_res.",colnames(Data@meta.data))]) {
	temp=Data
	Idents(temp)=i
	temp=FindAllMarkers(temp,only.pos = T)
	write.table(temp,paste0(i,".Cluster.Markers.tsv"),quote=F,sep="\t")

}
}



if ("cluster" %in% colnames(Data.markers)){
top10 = Data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
} else {
	top10=NULL
}
saveRDS(Data,"Final_Analysis.rds",compress = FALSE)

```

### Top gene markers for clusters

```{r Top gene markers for clusters, fig.width=15,fig.height=20,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}
#if (!is.null(top10)){
DT::datatable(top10,filter = "top",options = list(autoWidth = TRUE))
#}
```

### Heatmap of top gene markers for clusters

```{r heatmap of Top gene markers for clusters, fig.width=15,fig.height=20,error=FALSE,warning=FALSE,message=FALSE,echo=FALSE}
if (!is.null(top10)){

Vis=sub_sampling(Data)


Vis=ScaleData(Vis,features = top10$gene,assay = "RNA")
DoHeatmap(Vis, features = top10$gene,assay = "RNA") + NoLegend()
}
```


## Cluster Results Quality Control {.tabset}

In this section the number of genes, UMIs, mitochondrial percentages and ribosomal percentages will be plotted by cluster. This step is to check whether the clustering results are significantly biased by the sequencing depth, sequencing coverage and cell viability.

However, researches have shown (insert reference here later) that number of genes, UMIs and mitochondrial contents will vary between cell types and sub-populations. 

### Violin Plots {.tabset}

#### number of genes per cluster

```{r number of genes per cluster v, fig.width=10,fig.height =10}

VlnPlot(Data,"nFeature_RNA",group.by="seurat_clusters",pt.size = 0)+NoLegend()

```

#### number of UMIs per cluster

```{r number of UMIs per cluster v, fig.width=10,fig.height =10}

VlnPlot(Data,"nCount_RNA",group.by="seurat_clusters",pt.size = 0)+NoLegend()

```

#### mitochondrial percentages per cluster

```{r mitochondrial percentages per cluster v, fig.width=10,fig.height =10}

VlnPlot(Data,"percent.mt",group.by="seurat_clusters",pt.size = 0)+NoLegend()

```

#### ribosomal percentages per cluster

```{r ribosomal percentages per cluster v, fig.width=10,fig.height =10}

VlnPlot(Data,"percent.ribo",group.by="seurat_clusters",pt.size = 0)+NoLegend()

```

### Ridge Plots {.tabset}

#### number of genes per cluster

```{r number of genes per cluster r, fig.width=10,fig.height =10}

RidgePlot(Data,"nFeature_RNA",group.by="seurat_clusters")+NoLegend()

```

#### number of UMIs per cluster

```{r number of UMIs per cluster r, fig.width=10,fig.height =10}

RidgePlot(Data,"nCount_RNA",group.by="seurat_clusters")+NoLegend()

```

#### mitochondrial percentages per cluster

```{r mitochondrial percentages per cluster r, fig.width=10,fig.height =10}

RidgePlot(Data,"percent.mt",group.by="seurat_clusters")+NoLegend()

```

#### ribosomal percentages per cluster

```{r ribosomal percentages per cluster r, fig.width=10,fig.height =10}

RidgePlot(Data,"percent.ribo",group.by="seurat_clusters")+NoLegend()

```




EOF

open OUT, "> ClusterandMarker.rmd";
print OUT $script;
close OUT;

runCommand("Rscript -e 'rmarkdown::render(\\"ClusterandMarker.rmd\\",\\"html_document\\", output_file = \\"Final_Report.html\\",
params = list(SamplePath=\\"!{seurat_object}\\",minRes=as.numeric(!{minRes}),
maxRes=as.numeric(!{maxRes}),npcs=as.numeric(!{npcs})))'");

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
 path "*.h5ad"  ,emit:g36_22_h5_file00 

label 'scrna_seurat'

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

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30
}
//* platform
//* platform
//* autofill

process scRNA_Analysis_Module_seurat_to_sce {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /sce_obj.rds$/) "Analysis_Apps/$filename"}
input:
 path seurat_object

output:
 path "sce_obj.rds"  ,emit:g36_25_rdsFile00_g36_27 

label 'scrna_seurat'

script:
"""
#!/usr/bin/env Rscript

library(Seurat)
library(SingleCellExperiment) 

seurat = readRDS("${seurat_object}")

sce = as.SingleCellExperiment(seurat)
saveRDS(sce, file='sce_obj.rds')
"""
}


process scRNA_Analysis_Module_launch_isee {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /launch_iSEE.R$/) "Analysis_Apps/$filename"}
input:
 path sce_object

output:
 path 'launch_iSEE.R'  ,emit:g36_27_outputFileOut00 


shell:

'''
#!/usr/bin/env perl

my $script = <<'EOF';

library(iSEE)

getPath <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
        path <- dirname(normalizePath(sub(needle, "", cmdArgs[match]))[1])
        return(path)
    } else {
        return(normalizePath(getwd()))
    }
}

path = normalizePath(paste0(getPath()))
sce_file = "!{sce_object}"

sce = readRDS(normalizePath(paste0(path, '/', sce_file)))
app = iSEE(sce)
options(shiny.host = "0.0.0.0", shiny.port = 8789)
shiny::runApp(app)

EOF

open OUT, ">launch_iSEE.R";
print OUT $script;
close OUT;

'''
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
 path "Data.loom" ,optional:true  ,emit:g36_30_outputFileOut00 

label 'scrna_seurat'


when:
(params.run_scRNA_Analysis && (params.run_scRNA_Analysis == "yes")) || !params.run_scRNA_Analysis

shell:
Generate_loom_file = params.scRNA_Analysis_Module_SCEtoLOOM.Generate_loom_file

'''
#!/usr/bin/env Rscript

#library
library(Seurat)

if (as.logical("!{Generate_loom_file}")) {
	
Data=readRDS("!{seurat_object}")



annotation=read.csv("https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/raw/main/example_input_files/gene_info_table.csv",header = T,row.names = 1)
annotation=annotation[annotation$gene_name%in%names(table(annotation$gene_name))[table(annotation$gene_name)==1],]
rownames(annotation)=annotation$gene_name
metadata=Data@meta.data
matrix=Data@assays$RNA@counts

matrix=matrix[rowSums(matrix)>0,]

matrix=matrix[intersect(rownames(matrix),rownames(annotation)),]

annotation=annotation[intersect(rownames(matrix),rownames(annotation)),]

rownames(matrix)=annotation$ensembl_id
matrix=matrix[rownames(matrix)[order(rownames(matrix),decreasing = F)],]

NewData=CreateSeuratObject(matrix,meta.data = metadata)

NewData.loom <- SeuratDisk::as.loom(NewData, filename = "Data.loom", verbose = FALSE)


}
'''


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
g17_54_genome10_g_19 = Check_and_Build_Module_check_files.out.g17_54_genome10_g_19
g17_54_genomeSizes22 = Check_and_Build_Module_check_files.out.g17_54_genomeSizes22
g17_54_bed33 = Check_and_Build_Module_check_files.out.g17_54_bed33


cell_ranger_mkref(g17_54_genome10_g_19,g17_54_gtfFile01_g_19)
g_19_reference00_g_20 = cell_ranger_mkref.out.g_19_reference00_g_20

g_19_reference00_g_20= g_19_reference00_g_20.ifEmpty(ch_empty_file_1) 


cellranger_ref_checker(g_19_reference00_g_20)
g_20_reference02_g_5 = cellranger_ref_checker.out.g_20_reference02_g_5


Count(g_51_reads00_g_5,g_2_1_g_5,g_20_reference02_g_5)
g_5_outputDir00 = Count.out.g_5_outputDir00
g_5_outputHTML11 = Count.out.g_5_outputHTML11
g_5_outputDir22 = Count.out.g_5_outputDir22
g_5_outputDir33 = Count.out.g_5_outputDir33
g_5_h5_file40_g_49 = Count.out.g_5_h5_file40_g_49
g_5_h5_file51_g_49 = Count.out.g_5_h5_file51_g_49


Ambient_RNA_QC(g_5_h5_file40_g_49,g_5_h5_file51_g_49)
g_49_h5_file00_g36_0 = Ambient_RNA_QC.out.g_49_h5_file00_g36_0
g_49_outputFileHTML11 = Ambient_RNA_QC.out.g_49_outputFileHTML11



scRNA_Analysis_Module_Load_Data_h5(g_49_h5_file00_g36_0,g_41_1_g36_0)
g36_0_rdsFile00_g36_14 = scRNA_Analysis_Module_Load_Data_h5.out.g36_0_rdsFile00_g36_14
g36_0_outputFileHTML11 = scRNA_Analysis_Module_Load_Data_h5.out.g36_0_outputFileHTML11


scRNA_Analysis_Module_Merge_Seurat_Objects(g36_0_rdsFile00_g36_14.collect())
g36_14_rdsFile00_g36_17 = scRNA_Analysis_Module_Merge_Seurat_Objects.out.g36_14_rdsFile00_g36_17


scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction(g36_14_rdsFile00_g36_17)
g36_17_rdsFile00_g36_19 = scRNA_Analysis_Module_PCA_and_Batch_Effect_Correction.out.g36_17_rdsFile00_g36_19


scRNA_Analysis_Module_Clustering_and_Find_Markers(g36_17_rdsFile00_g36_19)
g36_19_outputHTML00 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g36_19_outputHTML00
g36_19_rdsFile10_g36_22 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g36_19_rdsFile10_g36_22
(g36_19_rdsFile10_g36_25,g36_19_rdsFile10_g36_30) = [g36_19_rdsFile10_g36_22,g36_19_rdsFile10_g36_22]
g36_19_outFileTSV22 = scRNA_Analysis_Module_Clustering_and_Find_Markers.out.g36_19_outFileTSV22


scRNA_Analysis_Module_Create_h5ad(g36_19_rdsFile10_g36_22)
g36_22_h5_file00 = scRNA_Analysis_Module_Create_h5ad.out.g36_22_h5_file00


scRNA_Analysis_Module_seurat_to_sce(g36_19_rdsFile10_g36_25)
g36_25_rdsFile00_g36_27 = scRNA_Analysis_Module_seurat_to_sce.out.g36_25_rdsFile00_g36_27


scRNA_Analysis_Module_launch_isee(g36_25_rdsFile00_g36_27)
g36_27_outputFileOut00 = scRNA_Analysis_Module_launch_isee.out.g36_27_outputFileOut00


scRNA_Analysis_Module_SCEtoLOOM(g36_19_rdsFile10_g36_30)
g36_30_outputFileOut00 = scRNA_Analysis_Module_SCEtoLOOM.out.g36_30_outputFileOut00


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
