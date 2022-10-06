shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc; ")

configfile: "config.yaml"
PREFIX = config['raw_reads']
SUFFIX = config['r1_ext']
IDS, = glob_wildcards(f'{PREFIX}{{id}}{SUFFIX}')

bwa_threads = config['bwa_threads']
cdhit_mem = config['cdhit_mem']
cdhit_threads = config['cdhit_threads']
cov_split = config['cov_split']
db = config['db']
dia_threads = config['dia_threads']
diamond_db = config['diamond_db']
kofamscan_split = config['kofamscan_split']
kofamscan_threads = config['kofamscan_threads']
max_qc = config['max_qc']
mega_kmers = config['mega_kmers']
mega_len = config['mega_len']
mega_mem = config['mega_mem']
mega_threads = config['mega_threads']
min_pid = config['min_pid']
min_qc = config['min_qc']
minlen = config['minlen']
mmseq_db = config['mmseq_db']
mmseq_keggs_threads = config['mmseq_keggs_threads']
mmseq_nohit_threads = config['mmseq_nohit_threads']
outdir = config['outdir']
overlap = config['overlap']
combined_bdg = config['combined_bdg']
polyg = config['polyg']
qual = config['qual']
r1_ext = config['r1_ext']
r1ad = config['r1ad']
r2_ext = config['r2_ext']
r2ad = config['r2ad']
raw_reads = config['raw_reads']
split_num = config['split_num']
tmp_cat = config['tmp_cat']
trim_threads = config['trim_threads']
unmapped_threads = config['unmapped_threads']
kallisto = config['kallisto']
kallisto_threads = config['kallisto_threads']
barrnap = config['barrnap']
nohit = config['nohit']
splitting_multiples = config['splitting_multiples']
grouping_multiples = config['grouping_multiples']

#default
rule all:
  input:
    expand(f"{outdir}Touch/trim_{{id}}", id=IDS),
    expand(f"{outdir}assemblies/{{id}}/{{id}}.contigs.fa", id=IDS),
    expand(f"{outdir}Prodigal_Outputs/filtered_aminos/{{id}}/{{id}}.faa", id=IDS),
    expand(f"{outdir}Abundance/{{id}}/{{id}}_bwa", id=IDS),
    "scripts/BamDeal-0.24/bin/BamDeal_Linux",
    expand(f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage", id=IDS),
    f"{db}ko_list",
    expand(f"{outdir}Kofamscan/{{id}}/{{id}}", id=IDS),
    f"{outdir}Results/KOunts_Kofamscan.csv", 
    f"{outdir}Touch/cdhit", 
    f"{outdir}Touch/split_kegg",
    expand("{outdir}Touch/mmseq_kegg_{split}", outdir=config["outdir"],split=config["split"]),
    f"{outdir}Touch/mmseq_nohit", 
    f"{outdir}Results/Number_of_clusters.csv", 
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_KOunt", id=IDS),
    f"{outdir}Touch/kegg_fasta",
    expand(f"{outdir}Touch/nohit_fa_{{id}}", id=IDS),
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohitp_KOunt", id=IDS),   
    expand(f"{outdir}barrnap/nohit_{{id}}_all_output", id=IDS),
    expand(f"{outdir}trnascan/{{id}}", id=IDS),
    expand(f"{outdir}Touch/{{id}}_nohit_fastq", id=IDS),
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_KOunts", id=IDS),
    expand(f"{outdir}kallisto/{{id}}_kallisto_still_missing_KOunt", id=IDS),
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_unmapped_KOunt", id=IDS),
    expand(f"{outdir}kallisto/{{id}}_unmapped_kallisto_still_missing_KOunt", id=IDS),
    expand(f"{outdir}Results/All_KOunts_nohit_unmapped_default.csv", id=IDS)

#not performing protein clustering
rule all_without_clustering:
  input:
    expand(f"{outdir}Touch/trim_{{id}}", id=IDS),
    expand(f"{outdir}assemblies/{{id}}/{{id}}.contigs.fa", id=IDS),
    expand(f"{outdir}Prodigal_Outputs/filtered_aminos/{{id}}/{{id}}.faa", id=IDS),
    expand(f"{outdir}Abundance/{{id}}/{{id}}_bwa", id=IDS),
    "scripts/BamDeal-0.24/bin/BamDeal_Linux",
    expand(f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage", id=IDS),
    f"{db}ko_list",
    expand(f"{outdir}Kofamscan/{{id}}/{{id}}", id=IDS),
    f"{outdir}Results/KOunts_Kofamscan_without_clustering.csv",
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noc_KOunt", id=IDS),
    expand(f"{outdir}Touch/nohit_fa_{{id}}", id=IDS),    
    expand(f"{outdir}barrnap/nohit_{{id}}_all_output", id=IDS),
    expand(f"{outdir}trnascan/{{id}}", id=IDS),
    expand(f"{outdir}Touch/{{id}}_nohit_fastq", id=IDS),
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_KOunts", id=IDS),
    expand(f"{outdir}kallisto/{{id}}_kallisto_still_missing_KOunt", id=IDS),
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_unmapped_KOunt", id=IDS),
    expand(f"{outdir}kallisto/{{id}}_unmapped_kallisto_still_missing_KOunt", id=IDS),
    expand(f"{outdir}Results/All_KOunts_nohit_unmapped_no_clustering.csv", id=IDS)

#without protein clustering and the reference database
rule all_without_reference:
  input:
    expand(f"{outdir}Touch/trim_{{id}}", id=IDS),
    expand(f"{outdir}assemblies/{{id}}/{{id}}.contigs.fa", id=IDS),
    expand(f"{outdir}Prodigal_Outputs/filtered_aminos/{{id}}/{{id}}.faa", id=IDS),
    expand(f"{outdir}Abundance/{{id}}/{{id}}_bwa", id=IDS),
    "scripts/BamDeal-0.24/bin/BamDeal_Linux",
    expand(f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage_noref", id=IDS),
    f"{db}ko_list",
    expand(f"{outdir}Kofamscan/{{id}}/{{id}}", id=IDS),
    f"{outdir}Results/All_KOunts_without_reference.csv",
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noref_KOunt", id=IDS)

#not performing protein clustering or calculating RNA abundance
rule all_without_RNA:
  input:
    expand(f"{outdir}Touch/trim_{{id}}", id=IDS),
    expand(f"{outdir}assemblies/{{id}}/{{id}}.contigs.fa", id=IDS),
    expand(f"{outdir}Prodigal_Outputs/filtered_aminos/{{id}}/{{id}}.faa", id=IDS),
    expand(f"{outdir}Abundance/{{id}}/{{id}}_bwa", id=IDS),
    "scripts/BamDeal-0.24/bin/BamDeal_Linux",
    expand(f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage", id=IDS),
    f"{db}ko_list",
    expand(f"{outdir}Kofamscan/{{id}}/{{id}}", id=IDS),
    f"{outdir}Results/KOunts_Kofamscan_without_clustering.csv",    
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noc_KOunt", id=IDS),
    expand(f"{outdir}Touch/nohit_fa_without_RNA_{{id}}", id=IDS),
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_without_RNA_KOunts", id=IDS),
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohitp_KOunt_norna", id=IDS),
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_unmapped_KOunt_no_RNA", id=IDS),
    expand(f"{outdir}Results/All_KOunts_without_RNA.csv", id=IDS),
    expand(f"{outdir}Touch/{{id}}_nohit_fastq_without_RNA", id=IDS)

rule trim:
  input:
    r1 = f'{raw_reads}{{id}}{r1_ext}',
    r2 = f'{raw_reads}{{id}}{r2_ext}'
  output: f'{outdir}Touch/trim_{{id}}'
  params:
    trim = f'{outdir}trimmed',
    trim_id = f'{outdir}trimmed/{{id}}',
    threads=f'{trim_threads}',
    r1_adapter = f'{r1ad}',
    r2_adapter = f'{r2ad}',
    polyG = f'{polyg}',
    quality = f'{qual}',
    min = f'{minlen}',
    over = f'{overlap}',
    r1 = f'{outdir}trimmed/{{id}}/{{id}}.trimmed.R1.fastq.gz',
    r2 = f'{outdir}trimmed/{{id}}/{{id}}.trimmed.R2.fastq.gz'
  conda: "envs/fastp.yaml"
  shell:
    '''
    mkdir -p {params.trim}
    mkdir -p {params.trim_id}
    fastp -i {input.r1} -o {params.r1} -I {input.r2} -O {params.r2} --adapter_sequence {params.r1_adapter} --adapter_sequence_r2 {params.r2_adapter} {params.quality} -l {params.min} -w {params.threads} --overlap_len_require {params.over} {params.polyG} --json /dev/null --html /dev/null
    touch {output}
    '''

rule megahit:
  input: f'{outdir}Touch/trim_{{id}}'
  params:
    r1 = f'{outdir}trimmed/{{id}}/{{id}}.trimmed.R1.fastq.gz',
    r2 = f'{outdir}trimmed/{{id}}/{{id}}.trimmed.R2.fastq.gz',
    outdir=f'{outdir}megahit/{{id}}/',
    mega=f'{outdir}megahit',
    memory=f'{mega_mem}',
    kmers=f'{mega_kmers}',
    len=f'{mega_len}',
    threads=f'{mega_threads}'
  output:
    final=f'{outdir}assemblies/{{id}}/{{id}}.contigs.fa'
  conda: "envs/Megahit.yaml"
  shell:
    '''
    mkdir -p {params.mega}
    megahit -m {params.memory} --kmin-1pass --k-list {params.kmers} --min-contig-len {params.len} -t {params.threads} -1 {params.r1} -2 {params.r2} -o {params.outdir}
    mv {params.outdir}final.contigs.fa {output.final}
    ./scripts/add_sample_name.sh {output.final} {wildcards.id}
    rm -r {params.outdir}
    '''

rule prodigal:
  input: f'{outdir}assemblies/{{id}}/{{id}}.contigs.fa'
  params: 
    output=f'{outdir}Prodigal_Outputs/{{id}}/'
  output:
    filtered_nuc=f'{outdir}Prodigal_Outputs/filtered_nucleotides/{{id}}/{{id}}.fa',
    filtered_aa=f'{outdir}Prodigal_Outputs/filtered_aminos/{{id}}/{{id}}.faa',
    filtered_gff=f'{outdir}Prodigal_Outputs/filtered_gffs/{{id}}/{{id}}.gff'
  conda: "envs/prodigal.yaml"
  shell:
    '''
    mkdir -p {params.output}
    prodigal -i {input} -p meta -o {params.output}{wildcards.id}.gff -d {params.output}{wildcards.id}.fa -a {params.output}{wildcards.id}.faa -f gff
    awk '/^>/ {{P=index($0,"partial=00")}} {{if(P) print}} ' {params.output}{wildcards.id}.fa > {output.filtered_nuc}
    awk '/^>/ {{P=index($0,"partial=00")}} {{if(P) print}} ' {params.output}{wildcards.id}.faa > {output.filtered_aa}
    sed -i 's/ .*$//' {output.filtered_nuc}
    sed -i 's/ .*$//' {output.filtered_aa}
    awk '!/^>/ {{ printf "%s", $0; n="\\n" }} /^>/ {{ print n $0; n="" }} END {{ printf "%s", n }}' {output.filtered_nuc} > {output.filtered_nuc}_tmp && mv {output.filtered_nuc}_tmp {output.filtered_nuc}
    awk '!/^>/ {{ printf "%s", $0; n="\\n" }} /^>/ {{ print n $0; n="" }} END {{ printf "%s", n }}' {output.filtered_aa} > {output.filtered_aa}_tmp && mv {output.filtered_aa}_tmp {output.filtered_aa}
    awk '{{P=index($0,"partial=00")}} {{if(P) print}}' {params.output}{wildcards.id}.gff | sed 's/ID=[^_]*_/ID=/g' > {output.filtered_gff}
    rm -r {params.output}
    '''

rule bwa:
  input: 
    fasta=f"{outdir}assemblies/{{id}}/{{id}}.contigs.fa",
    touch=f"{outdir}Touch/trim_{{id}}"
  output: 
    touchfile=f"{outdir}Abundance/{{id}}/{{id}}_bwa"
  params: 
    r1 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R1.fastq.gz",
    r2 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R2.fastq.gz",
    bam = f"{outdir}Abundance/{{id}}/{{id}}.bam",
    threads=f"{bwa_threads}"
  conda: "envs/bwa_map.yaml"
  shell:
    '''
    bwa index {input.fasta}
    bwa mem -t {params.threads} {input.fasta} {params.r1} {params.r2} | samtools view -bS - > {params.bam}
    rm {input.fasta}.*
    touch {output.touchfile}
    '''

rule bamdeal:
  output:
    bamdeal="scripts/BamDeal-0.24/bin/BamDeal_Linux"
  params:
    download="https://github.com/BGI-shenzhen/BamDeal/archive/v0.24.tar.gz",
    bamdeal_location="scripts/"
  shell:
    '''
    wget -P {params.bamdeal_location} {params.download}
    tar -C {params.bamdeal_location} -zxvf {params.bamdeal_location}v0.24.tar.gz
    chmod 755 {output.bamdeal}
    '''

rule coverage:
  input:
    filtered_gff=f"{outdir}Prodigal_Outputs/filtered_gffs/{{id}}/{{id}}.gff",
    bamdeal="scripts/BamDeal-0.24/bin/BamDeal_Linux",
    touchfile=f"{outdir}Abundance/{{id}}/{{id}}_bwa"
  output:
    coverage=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage"
  params:
    bam=f"{outdir}Abundance/{{id}}/{{id}}.bam",
    wd=f"{outdir}Abundance/{{id}}/{{id}}",
    abun=f"{outdir}Abundance/",
    out=f"{outdir}Abundance_Outputs/{{id}}/",
    stats=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_stats",
    gff=f"{outdir}Prodigal_Outputs/filtered_gffs/",
    split=f"{cov_split}",
    mapped=f"{outdir}Abundance/{{id}}/{{id}}_mapped_reads"
  threads:1
  conda: "envs/bedtools_samtools.yaml"
  shell:
    '''
    mkdir -p {params.wd}
    mkdir -p {params.out}
    awk '{{print $1}}' {input.filtered_gff} | uniq > {params.wd}_contigs
    split -l {params.split} {params.wd}_contigs {params.wd}
    for file in {params.wd}??; do awk -v file=$file '{{print $1 "\t" file}}' $file | awk -F'/|\t' '{{print $1 "\t" $NF}}'; done > {params.wd}_bedtools_list
    {input.bamdeal} modify bamAssign -i {params.bam} -a {params.wd}_bedtools_list -q 0 -o {params.wd}
    for file in {params.wd}??; do ./scripts/bedtools.sh $file {params.gff} {params.wd}; done > {output.coverage}
    bedtools bamtobed -i {params.bam} > {params.bam}.bed
    rm {params.bam}
    sort -k1,1 -k2,2n {params.bam}.bed > {params.bam}.bed_sorted
    sort -k1,1 -k4,4n {input.filtered_gff} > {input.filtered_gff}_sorted
    bedtools intersect -a {input.filtered_gff}_sorted -b {params.bam}.bed_sorted -wb -sorted | awk '{{print $13}}' | sort -S 50% | uniq > {params.mapped}
    rm {input.filtered_gff}_sorted
    rm {params.bam}.bed_sorted
    rm {params.wd}/NaAss.bam
    rm {params.wd}/UnMap.bam
    rm {params.wd}_contigs
    rm {params.wd}_bedtools_list
    rm {params.wd}*coverage
    rm {params.wd}*gff
    rm {params.wd}??
    '''

rule coverage_without_reference:
  input:
    filtered_gff=f"{outdir}Prodigal_Outputs/filtered_gffs/{{id}}/{{id}}.gff",
    bamdeal="scripts/BamDeal-0.24/bin/BamDeal_Linux",
    touchfile=f"{outdir}Abundance/{{id}}/{{id}}_bwa"
  output:
    coverage=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage_noref"
  params:
    bam=f"{outdir}Abundance/{{id}}/{{id}}.bam",
    wd=f"{outdir}Abundance/{{id}}/{{id}}",
    abun=f"{outdir}Abundance/",
    out=f"{outdir}Abundance_Outputs/{{id}}/",
    stats=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_stats",
    gff=f"{outdir}Prodigal_Outputs/filtered_gffs/",
    split=f"{cov_split}"
  threads:1
  conda: "envs/bedtools_samtools.yaml"
  shell:
    '''
    mkdir -p {params.wd}
    mkdir -p {params.out}
    awk '{{print $1}}' {input.filtered_gff} | uniq > {params.wd}_contigs
    split -l {params.split} {params.wd}_contigs {params.wd}
    for file in {params.wd}??; do awk -v file=$file '{{print $1 "\t" file}}' $file | awk -F'/|\t' '{{print $1 "\t" $NF}}'; done > {params.wd}_bedtools_list
    {input.bamdeal} modify bamAssign -i {params.bam} -a {params.wd}_bedtools_list -q 0 -o {params.wd}
    for file in {params.wd}??; do ./scripts/bedtools.sh $file {params.gff} {params.wd}; done > {output.coverage}
    rm {params.bam}
    rm -r {params.wd}
    rm {params.wd}_contigs
    rm {params.wd}_bedtools_list
    rm {params.wd}*coverage
    rm {params.wd}*gff
    rm {params.wd}??
    '''

rule kegg_db:
  output:
    profiles=directory(f"{db}profiles/"),
    list=f"{db}ko_list",
    touch=f"{db}touch"
  params:
    profiles="ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz",
    ko_list="ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz",
    location=f"{db}"
  shell:
    '''
    touch {output.touch}
    wget -P {params.location} {params.profiles}
    wget -P {params.location} {params.ko_list}
    tar -C {params.location} -xf {params.location}/profiles.tar.gz
    gunzip {params.location}/ko_list.gz
    '''

rule kofamscan:
  input:
    filtered_aa=f"{outdir}Prodigal_Outputs/filtered_aminos/{{id}}/{{id}}.faa"
  params:
    profiles=f"{db}profiles/",
    list=f"{db}ko_list",
    tmp=f"{outdir}Kofamscan/tmp_{{id}}",
    threads=f"{kofamscan_threads}"
  output:
    ko=f"{outdir}Kofamscan/{{id}}/{{id}}"
  conda:"envs/kofamscan_hmmr.yaml"
  shell:
    '''
    exec_annotation --cpu {params.threads} -p {params.profiles} -k {params.list} -o {output.ko} --tmp-dir {params.tmp} {input.filtered_aa}
    rm -r {params.tmp}
    '''

rule kofamscan_results:
  input:
    ko=f"{outdir}Kofamscan/{{id}}/{{id}}",
    coverage=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage",
    filtered_aa=f"{outdir}Prodigal_Outputs/filtered_aminos/{{id}}/{{id}}.faa"
  params:
    split_keggs=f"{outdir}CDHit/",
    splitting_multiples=f"{splitting_multiples}",
    grouping_multiples=f"{grouping_multiples}"
  output:
    out=f"{outdir}Kofamscan/Results/{{id}}/{{id}}",
    KOunt=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_KOunt",
    nohit=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit"
  shell:
    '''
    grep '^\*' {input.ko} | awk '{{print $2 "\t" $3}}' > {output.out}
    {params.splitting_multiples}./scripts/splitting_multiples.sh {input.coverage} {output.out} {output.out} #include if splitting multiples
    {params.grouping_multiples}awk '{{print $1}}' {output.out} | sort | uniq -d > {output.out}_duplicates #include if grouping multiples
    {params.grouping_multiples}cat {output.out} | parallel --colsep '\t' ./scripts/multiples.sh {{1}} {output.out} {output.out}_duplicates | sort | uniq > {output.out}_multiples #include if grouping multiples
    {params.grouping_multiples}awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' {output.out}_multiples {input.coverage} | awk 'BEGIN{{OFS=FS}} $3 == "" {{$3 = "NoHit"}} 1' > {output.out}_coverage #include if grouping multiples
    {params.grouping_multiples}rm {output.out}_duplicates #include if grouping multiples
    {params.grouping_multiples}rm {output.out}_multiples #include if grouping multiples
    {params.grouping_multiples}awk '{{print $2 "\t" $3}}' {output.out}_coverage | awk '{{a[$2]+=$1}}END{{for(i in a) print i,a[i]}}' > {output.KOunt} #include if grouping multiples
    mkdir -p {params.split_keggs}
    awk '{{print $3}}' {output.out}_coverage > {output.out}_coverage_tmp && while read l; do mkdir -p {params.split_keggs}/"$l"; done < {output.out}_coverage_tmp && rm {output.out}_coverage_tmp
    awk '{{print $1>f"{params.split_keggs}" $3 "/{wildcards.id}_" $3}}' {output.out}_coverage
    for file in {params.split_keggs}*/{wildcards.id}_*; do ./scripts/extracting_genes.sh $file {input.filtered_aa}; done
    grep NoHit {output.out}_coverage | awk '{{print $1}}' > {output.nohit}
    rm {output.out}_coverage
    '''

rule kofamscan_results_no_clustering:
  input:
    ko=f"{outdir}Kofamscan/{{id}}/{{id}}",
    coverage=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage"
  params:
    splitting_multiples=f"{splitting_multiples}",
    grouping_multiples=f"{grouping_multiples}"
  output:
    out=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noc",
    KOunt=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noc_KOunt",
    nohit=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_noc"
  shell:
    '''
    grep '^\*' {input.ko} | awk '{{print $2 "\t" $3}}' > {output.out}
    {params.splitting_multiples}./scripts/splitting_multiples.sh {input.coverage} {output.out} {output.out} #include if splitting multiples
    {params.grouping_multiples}awk '{{print $1}}' {output.out} | sort | uniq -d > {output.out}_duplicates #include if grouping multiples
    {params.grouping_multiples}cat {output.out} | parallel --colsep '\t' ./scripts/multiples.sh {{1}} {output.out} {output.out}_duplicates | sort | uniq > {output.out}_multiples #include if grouping multiples
    {params.grouping_multiples}awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' {output.out}_multiples {input.coverage} | awk 'BEGIN{{OFS=FS}} $3 == "" {{$3 = "NoHit"}} 1' > {output.out}_coverage #include if grouping multiples
    {params.grouping_multiples}rm {output.out}_duplicates #include if grouping multiples
    {params.grouping_multiples}rm {output.out}_multiples #include if grouping multiples
    {params.grouping_multiples}awk '{{print $2 "\t" $3}}' {output.out}_coverage | awk '{{a[$2]+=$1}}END{{for(i in a) print i,a[i]}}' > {output.KOunt} #include if grouping multiples
    grep NoHit {output.out}_coverage | awk '{{print $1}}' > {output.nohit}
    rm {output.out}_coverage
    '''

rule kofamscan_results_noref:
  input:
    ko=f"{outdir}Kofamscan/{{id}}/{{id}}",
    coverage=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage_noref"
  params:
    splitting_multiples=f"{splitting_multiples}",
    grouping_multiples=f"{grouping_multiples}"
  output:
    out=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noref",
    KOunt=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noref_KOunt",
    nohit=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_noref"
  shell:
    '''
    grep '^\*' {input.ko} | awk '{{print $2 "\t" $3}}' > {output.out}
    {params.splitting_multiples}./scripts/splitting_multiples.sh {input.coverage} {output.out} {output.out} #include if splitting multiples
    {params.grouping_multiples}awk '{{print $1}}' {output.out} | sort | uniq -d > {output.out}_duplicates #include if grouping multiples
    {params.grouping_multiples}cat {output.out} | parallel --colsep '\t' ./scripts/multiples.sh {{1}} {output.out} {output.out}_duplicates | sort | uniq > {output.out}_multiples #include if grouping multiples
    {params.grouping_multiples}awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' {output.out}_multiples {input.coverage} | awk 'BEGIN{{OFS=FS}} $3 == "" {{$3 = "NoHit"}} 1' > {output.out}_coverage #include if grouping multiples
    {params.grouping_multiples}rm {output.out}_duplicates #include if grouping multiples
    {params.grouping_multiples}rm {output.out}_multiples #include if grouping multiples
    {params.grouping_multiples}awk '{{print $2 "\t" $3}}' {output.out}_coverage | awk '{{a[$2]+=$1}}END{{for(i in a) print i,a[i]}}' > {output.KOunt} #include if grouping multiples
    grep NoHit {output.out}_coverage | awk '{{print $1}}' > {output.nohit}
    rm {output.out}_coverage
    '''

rule abundance_results:
  input:
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_KOunt", id=IDS)
  params:
    tmp=f"{outdir}Kofamscan/tmp/",
    list=f"{db}ko_list",
    kegglist=f"{outdir}kegglist"
  output:
    results=f"{outdir}Results/KOunts_Kofamscan.csv"
  shell:
    '''
    mkdir -p {params.tmp}
    if [ ! -e {params.kegglist} ]; then cat {input} | awk '{{print $1}}' | sort | uniq > {params.kegglist}; fi
    for i in {input}; do cp $i {params.tmp}; done
    ./scripts/abundance_matrix.sh {params.tmp} {params.list} {params.kegglist} > {output.results}
    rm -r {params.tmp}
    if [ -e {outdir}Touch/kegg_fasta ]; then rm {params.kegglist}; fi
    '''

rule abundance_results_noc:
  input:
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noc_KOunt", id=IDS)
  params:
    tmp=f"{outdir}Kofamscan/tmp/",
    list=f"{db}ko_list",
    kegglist=f"{outdir}kegglist"
  output:
    results=f"{outdir}Results/KOunts_Kofamscan_without_clustering.csv"
  shell:
    '''
    mkdir -p {params.tmp}
    if [ ! -e {params.kegglist} ]; then cat {input} | awk '{{print $1}}' | sort | uniq > {params.kegglist}; fi
    for i in {input}; do cp $i {params.tmp}; done
    ./scripts/abundance_matrix_noc.sh {params.tmp} {params.list} {params.kegglist} > {output.results}
    rm -r {params.tmp}
    if [ -e {outdir}Touch/kegg_fasta ]; then rm {params.kegglist}; fi
    '''

rule abundance_results_noref:
  input:
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noref_KOunt", id=IDS)
  params:
    tmp=f"{outdir}Kofamscan/tmp/",
    list=f"{db}ko_list",
    kegglist=f"{outdir}kegglist"
  output:
    results=f"{outdir}Results/All_KOunts_without_reference.csv"
  shell:
    '''
    mkdir -p {params.tmp}
    if [ ! -e {params.kegglist} ]; then cat {input} | awk '{{print $1}}' | sort | uniq > {params.kegglist}; fi
    for i in {input}; do cp $i {params.tmp}; done
    ./scripts/abundance_matrix_noref.sh {params.tmp} {params.list} {params.kegglist} > {output.results}
    rm -r {params.tmp}
    if [ -e {outdir}Touch/kegg_fasta ]; then rm {params.kegglist}; fi
    '''

rule kegg_fastas:
  input:
    expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_KOunt", id=IDS)
  params:
    split_dir=f"{outdir}CDHit/",
    cd=f"{outdir}CD-Hit_Outputs",
    cdhit=f"{outdir}CD-Hit_Outputs/By_Sample/",
    mega=f"{outdir}megahit",
    results=f"{outdir}Kofamscan/Results/all_KOunts.csv",
    kegglist=f"{outdir}kegglist"
  output: f"{outdir}Touch/kegg_fasta"
  shell:
    '''
    mkdir -p {params.cd}
    mkdir -p {params.cdhit}
    if [ ! -e {params.kegglist} ]; then cat {input} | awk '{{print $1}}' | sort | uniq > {params.kegglist}; fi
    cat {params.kegglist} | parallel "cat {params.split_dir}{{1}}/*faa > {params.cdhit}{{1}}.faa"
    rm -r {params.split_dir}
    if [ -e {params.mega} ]; then rm -r {params.mega}; fi
    touch {output}
    if [ -e {params.results} ]; then rm {params.kegglist}; fi
    '''

rule cdhit:
  input:
    f"{outdir}Touch/kegg_fasta"
  params:
    cdhit=f"{outdir}CD-Hit_Outputs/By_Sample/",
    memory=f"{cdhit_mem}",
    threads=f"{cdhit_threads}"
  output: f"{outdir}Touch/cdhit"
  conda:"envs/CDHit.yaml"
  
  shell:
    '''
    for file in {params.cdhit}*faa; do cd-hit -i $file -o "$file"_1.0 -c 1.0 -d 0 -M {params.memory} -T {params.threads}; done
    for file in {params.cdhit}*clstr; do clstr2txt.pl $file | awk '{{print $2 "\t" $1}}' > "$file".txt; done
    for file in {params.cdhit}*clstr; do clstr_rep.pl $file > "$file"_reps; done
    touch {output}
    '''

rule split_keggs:
  input: f"{outdir}Touch/cdhit"
  output: f"{outdir}Touch/split_kegg"
  params:
    splitting_multiples=f"{splitting_multiples}",
    grouping_multiples=f"{grouping_multiples}",
    sample=f"{outdir}CD-Hit_Outputs/By_Sample/",
    num=f"{split_num}",
    lists=f"{outdir}lists"
  shell: 
    '''
    mkdir -p {params.lists}
    {params.splitting_multiples}(cd {params.sample} && ls K*_1.0) > keggslist #include if splitting multiples
    {params.grouping_multiples}(cd {params.sample} && ls Multiple*1.0) >> keggslist #include if grouping multiples
    split -d --number=l/{params.num} keggslist {params.lists}/keggslist
    rm keggslist
    touch {output}
    '''

rule mmseq_keggs:
  input:
    keggs=f"{outdir}Touch/split_kegg",
    cdhit=f"{outdir}Touch/cdhit"
  params:
    cdhit=f"{outdir}CD-Hit_Outputs/By_Sample/",
    list="{outdir}lists/keggslist{split}",
    mm=f"{outdir}MMSeq_Outputs/",
    mmseq=f"{outdir}MMSeq_Outputs/By_Sample/",
    threads=f"{mmseq_keggs_threads}"
  output: "{outdir}Touch/mmseq_kegg_{split}"
  conda:"envs/mmseq.yaml"
  shell:
    '''
    mkdir -p {params.mm}
    mkdir -p {params.mmseq}
    while read file; do ./scripts/mmseq_keggs.sh $file {params.cdhit} {params.mmseq}; done < {params.list}
    while read file; do ./scripts/count_mmseq_clusters.sh {params.mmseq}"$file"_0.9_all; done < {params.list}
    while read file; do ./scripts/count_mmseq_clusters.sh {params.mmseq}"$file"_0.5_all; done < {params.list}
    rm {params.list}
    touch {output}
    '''

rule mmseq_nohit:
  input:
    cdhit=f"{outdir}Touch/cdhit"
  params:
    cdhit=f"{outdir}CD-Hit_Outputs/By_Sample/",
    mm=f"{outdir}MMSeq_Outputs/",
    mmseq=f"{outdir}MMSeq_Outputs/By_Sample/",
    threads=mmseq_nohit_threads
  output: f"{outdir}Touch/mmseq_nohit"
  conda:"envs/mmseq.yaml"
  shell:
    '''
    mkdir -p {params.mm}
    mkdir -p {params.mmseq}
    ./scripts/mmseq_keggs.sh NoHit.faa_1.0 {params.cdhit} {params.mmseq}
    ./scripts/count_mmseq_clusters.sh {params.mmseq}NoHit.faa_1.0_0.9_all
    ./scripts/count_mmseq_clusters.sh {params.mmseq}NoHit.faa_1.0_0.5_all
    touch {output}
    '''

rule mmseq_cluster_count:
  input:
    keggs=expand("{outdir}Touch/mmseq_kegg_{split}", outdir=config["outdir"], split=config['split']),
    nohit=f"{outdir}Touch/mmseq_nohit"
  params:
    mmseq=f"{outdir}MMSeq_Outputs/By_Sample/"
  output:f"{outdir}Results/Number_of_clusters.csv"
  shell:
    '''
    echo -e "KEGG,Number of clusters,Number of singletons,Number with multiple genes,Number of proteins\n$(cat {params.mmseq}*txt)" > {output}
    rm {params.mmseq}*txt
    '''

rule nohit_fastas:
  input:
    nohit=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit",
    faa=f"{outdir}Prodigal_Outputs/filtered_aminos/{{id}}/{{id}}.faa",
    fa=f"{outdir}Prodigal_Outputs/filtered_nucleotides/{{id}}/{{id}}.fa"
  params:
    faa=f"{outdir}diamond/NoHit_{{id}}.faa",
    fa=f"{outdir}diamond/NoHit_{{id}}.fa",
    wd=f"{outdir}diamond/"
  output:f"{outdir}Touch/nohit_fa_{{id}}"
  shell:
    '''
    mkdir -p {params.wd}
    ./scripts/extracting_fastas.sh {input.nohit} {input.faa} > {params.faa}
    ./scripts/extracting_fastas.sh {input.nohit} {input.fa} > {params.fa}
    touch {output}
    '''

rule nohit_fastas_without_RNA:
  input:
    nohit=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_noc",
    faa=f"{outdir}Prodigal_Outputs/filtered_aminos/{{id}}/{{id}}.faa"
  params:
    faa=f"{outdir}diamond/NoHit_{{id}}.faa",
    fa=f"{outdir}diamond/NoHit_{{id}}.fa",
    wd=f"{outdir}diamond/"
  output:f"{outdir}Touch/nohit_fa_without_RNA_{{id}}"
  shell:
    '''
    mkdir -p {params.wd}
    ./scripts/extracting_fastas.sh {input.nohit} {input.faa} > {params.faa}
    touch {output}
    '''

rule diamond_search:
  input:
    nohit=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit",
    cov=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage",
    touch=f"{outdir}Touch/nohit_fa_{{id}}"
  output: 
    raw=f"{outdir}diamond/NoHit_{{id}}",
    fil=f"{outdir}diamond/{{id}}/{{id}}_filtered",
    best=f"{outdir}diamond/NoHit_{{id}}_besthits",
    ko=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohitp_ko",
    KOunt=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohitp_KOunt",
    no=f"{outdir}diamond/still_nohit_{{id}}",
    mmseq=f"{outdir}mmseqs/nohitp_{{id}}"
  params:
    faa=f"{outdir}diamond/NoHit_{{id}}.faa",
    fa=f"{outdir}diamond/NoHit_{{id}}.fa",
    rna=f"{mmseq_db}",
    of="6 qseqid sseqid pident qlen slen qstart qend sstart send evalue",
    ofm="query,target,pident,qlen,tlen,qstart,qend,tstart,tend,evalue",
    db=f"{diamond_db}",
    tmp=f"{outdir}diamond/{{id}}_tmp",
    cat=f"{outdir}diamond/NoHit_{{id}}_cat",
    min_qc=f"{min_qc}",
    max_qc=f"{max_qc}",
    min_pid=f"{min_pid}",
    threads=f"{dia_threads}"
  conda: "envs/diamond_bedtools_mmseq2.yaml"
  shell:
    '''
    diamond blastp --query {params.faa} --threads {params.threads} --outfmt {params.of} -d {params.db} > {output.raw}
    mmseqs easy-search {params.fa} {params.rna} {output.mmseq} {params.tmp} --threads {params.threads} --search-type 3 --format-output {params.ofm}
    cat {output.raw} {output.mmseq} > {params.cat}
    awk -F '\t' '{{ print $0 "\t" ($4/$5) * 100 }}' {params.cat} | gawk -F '\t' '$NF>{params.min_qc} {{print $0}}' | gawk -F '\t' '$NF<{params.max_qc} {{print $0}}' | gawk -F '\t' '$3>{params.min_pid} {{print $0}}' > {output.fil}
    awk '! a[$1]++' {output.fil} > {output.fil}_tmp
    awk '{{print $1 "," $3 "," $10 "\t" $2}}' {output.fil}_tmp > {output.fil}_tmp.1
    awk '{{print $1 "," $3 "," $10 "\t" $2}}' {output.fil} > {output.fil}.1
    awk 'NR==FNR{{a[$1];next}}$1 in a' {output.fil}_tmp.1 {output.fil}.1 | tr ',' '\t' > {output.best}
    rm {output.fil}_tmp.1
    rm {output.fil}.1
    awk -F '_' '{{print $0 "\t" $NF}}' {output.best} | awk '{{print $1 "\t" $NF}}' | sort | uniq > {output.ko}
    awk '{{print $1}}' {output.ko} | uniq -c | awk '{{print $2 "\t" $1}}' > {output.ko}_count
    awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' {input.cov} {output.ko} | awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' {output.ko}_count - | awk '{{print $1 "\t" $2 "\t" ($3/$4)}}' > {output.ko}_coverage
    awk '{{a[$2]+=$3}}END{{for(i in a) print i,a[i]}}' {output.ko}_coverage | awk -F'-' '{{print $0 "\t" NF-1}}' | awk -v OFMT='%.5f' '{{print $1"\t",($2/$3)}}' | cut -d- -f2- | awk '{{gsub(/-/,"\t"$2",")}}1' | tr ',' '\n' | awk -v OFMT='%.5f' '{{a[$1]+=$2}}END{{for(i in a) print i,a[i]}}' > {output.KOunt}
    awk '{{print $1}}' {output.ko}_coverage > {output.ko}_coverage_tmp
    cat {input.nohit} {output.ko}_coverage_tmp | sort | uniq -u > {output.no}
    rm {output.ko}_coverage_tmp
    rm {output.ko}_coverage
    rm {output.ko}_count
    '''

rule diamond_search_without_RNA:
  input:
    nohit=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_noc",
    cov=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage",
    touch=f"{outdir}Touch/nohit_fa_without_RNA_{{id}}"
  output: 
    raw=f"{outdir}diamond/NoHit_{{id}}_norna",
    fil=f"{outdir}diamond/{{id}}/{{id}}_filtered_norna",
    best=f"{outdir}diamond/NoHit_{{id}}_besthits_norna",
    ko=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohitp_ko_norna",
    KOunt=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohitp_KOunt_norna",
    no=f"{outdir}diamond/still_nohit_{{id}}_norna"
  params:
    faa=f"{outdir}diamond/NoHit_{{id}}.faa",
    rna=f"{mmseq_db}",
    of="6 qseqid sseqid pident qlen slen qstart qend sstart send evalue",
    db=f"{diamond_db}",
    tmp=f"{outdir}diamond/{{id}}_tmp",
    cat=f"{outdir}diamond/NoHit_{{id}}_cat",
    min_qc=f"{min_qc}",
    max_qc=f"{max_qc}",
    min_pid=f"{min_pid}",
    threads=f"{dia_threads}"
  conda: "envs/diamond_bedtools_mmseq2.yaml"
  shell:
    '''
    diamond blastp --query {params.faa} --threads {params.threads} --outfmt {params.of} -d {params.db} > {output.raw}
    awk -F '\t' '{{ print $0 "\t" ($4/$5) * 100 }}' {output.raw} | gawk -F '\t' '$NF>{params.min_qc} {{print $0}}' | gawk -F '\t' '$NF<{params.max_qc} {{print $0}}' | gawk -F '\t' '$3>{params.min_pid} {{print $0}}' > {output.fil}
    awk '! a[$1]++' {output.fil} > {output.fil}_tmp
    awk '{{print $1 "," $3 "," $10 "\t" $2}}' {output.fil}_tmp > {output.fil}_tmp.1
    awk '{{print $1 "," $3 "," $10 "\t" $2}}' {output.fil} > {output.fil}.1
    awk 'NR==FNR{{a[$1];next}}$1 in a' {output.fil}_tmp.1 {output.fil}.1 | tr ',' '\t' > {output.best}
    rm {output.fil}_tmp.1
    rm {output.fil}.1
    awk -F '_' '{{print $0 "\t" $NF}}' {output.best} | awk '{{print $1 "\t" $NF}}' | sort | uniq > {output.ko}
    awk '{{print $1}}' {output.ko} | uniq -c | awk '{{print $2 "\t" $1}}' > {output.ko}_count
    awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' {input.cov} {output.ko} | awk 'NR==FNR{{a[$1]=$2;next}}{{print $0,a[$1]}}' {output.ko}_count - | awk '{{print $1 "\t" $2 "\t" ($3/$4)}}' > {output.ko}_coverage
    awk '{{a[$2]+=$3}}END{{for(i in a) print i,a[i]}}' {output.ko}_coverage | awk -F'-' '{{print $0 "\t" NF-1}}' | awk -v OFMT='%.5f' '{{print $1"\t",($2/$3)}}' | cut -d- -f2- | awk '{{gsub(/-/,"\t"$2",")}}1' | tr ',' '\n' | awk -v OFMT='%.5f' '{{a[$1]+=$2}}END{{for(i in a) print i,a[i]}}' > {output.KOunt}
    awk '{{print $1}}' {output.ko}_coverage > {output.ko}_coverage_tmp
    cat {input.nohit} {output.ko}_coverage_tmp | sort | uniq -u > {output.no}
    rm {output.ko}_coverage_tmp
    rm {output.ko}_coverage
    rm {output.ko}_count
    '''

rule barrnap:
  input:
    ann=f"{outdir}diamond/still_nohit_{{id}}",
    fa=f"{outdir}Prodigal_Outputs/filtered_nucleotides/{{id}}/{{id}}.fa"
  params:
    nohit=f"{outdir}diamond/NoHit_{{id}}.fa",
    threads=f"{barrnap}"
  conda:"envs/barrnap.yaml"
  output:
    mito=f"{outdir}barrnap/nohit_{{id}}_mito_barrnap.gff",
    bac=f"{outdir}barrnap/nohit_{{id}}_bac_barrnap.gff",
    arc=f"{outdir}barrnap/nohit_{{id}}_arc_barrnap.gff",
    euk=f"{outdir}barrnap/nohit_{{id}}_euk_barrnap.gff",
    all=f"{outdir}barrnap/nohit_{{id}}_all_output"
  shell:
    '''
    ./scripts/extracting_fastas.sh {input.ann} {input.fa} > {params.nohit}
    barrnap --threads {params.threads} --kingdom mito {params.nohit} > {output.mito}
    barrnap --threads {params.threads} --kingdom bac {params.nohit} > {output.bac}
    barrnap --threads {params.threads} --kingdom arc {params.nohit} > {output.arc}
    barrnap --threads {params.threads} --kingdom euk {params.nohit} > {output.euk}
    cat {output.mito} {output.bac} {output.arc} {output.euk} | grep -v gff | awk -F'\t|;' 'BEGIN{{OFS="\t"}} NR>1{{print $1,$4,$5,$9}}' > {output.all}
    '''

rule trnascan:
  input:f"{outdir}barrnap/nohit_{{id}}_all_output"
  params:f"{outdir}diamond/NoHit_{{id}}.fa"
  output:f"{outdir}trnascan/{{id}}"
  conda:"envs/trnascan.yaml"
  shell:
    '''
    tRNAscan-SE -G -o {output} {params}
    '''

rule rna_abundance:
  input:
    trna=f"{outdir}trnascan/{{id}}",
    rrna=f"{outdir}barrnap/nohit_{{id}}_all_output",
    gff=f"{outdir}Prodigal_Outputs/filtered_gffs/{{id}}/{{id}}.gff",
    nohit=f"{outdir}diamond/still_nohit_{{id}}"
  conda: "envs/diamond_bedtools_mmseq2.yaml"
  params:
    bed=f"{outdir}barrnap/nohit_{{id}}.bed",
    bedko=f"{outdir}barrnap/nohit_{{id}}_all_output_catKO",
    bedfull=f"{outdir}barrnap/nohit_{{id}}_full.bed",
    bedco=f"{outdir}barrnap/nohit_{{id}}_all_output_catKO_ID",
    bam=f"{outdir}Abundance/{{id}}/{{id}}.bam.bed",
    rnahits=f"{outdir}barrnap/RNA_{{id}}_hits",
    index="scripts/index",
    cov=f"{outdir}barrnap/nohit_{{id}}_coverage",
    miss=f"{outdir}barrnap/RNA_{{id}}_missing"
  output:
    raw=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_rna_coverage",
    miss=f"{outdir}barrnap/RNA_{{id}}_missing",
    kos=f"{outdir}barrnap/{{id}}_rna_kos",
    KOunt=f"{outdir}barrnap/{{id}}_rna_KOunt"
  shell:
    '''
    ./scripts/rrna_trna_cov.sh {input.rrna} {input.trna} {params.index} > {params.bed}
    while read file; do grep -w $(echo "$file" | awk '{{print $1}}' | rev | cut -d_ -f2- | rev) {input.gff} | grep $(echo "$file" | awk '{{print $1}}' | awk -F'_' '{{print "ID="$NF";"}}') | awk -F '[\t|=|;]' -v "file=$file" '{{print $1 "\t" $4 "\t" file}}' | awk '{{print $1"\t"($2+$4)"\t"($2+$5)"\t"$7"\t"$3}}'; done < {params.bed} > {params.bedfull}
    while read file; do grep -w $(echo "$file" | awk '{{print $1}}' | rev | cut -d_ -f2- | rev) {input.gff} | grep $(echo "$file" | awk '{{print $1}}' | awk -F'_' '{{print "ID="$NF";"}}') | awk -F '[\t|=|;]' -v "file=$file" '{{print $1 "\t" $4 "\t" file}}' | awk '{{print $1"\t"($2+$4)"\t"($2+$5)"\t"$6}}'; done < {params.bedko} > {params.bedco}
    bedtools coverage -a {params.bedfull} -b {params.bam} -d > {output.raw}
    bedtools intersect -a {output.raw} -b {params.bedco} > {params.cov}
    awk '{{print $5}}' {params.cov} | sort | uniq > {params.rnahits}
    cat {params.rnahits} {input.nohit} | sort | uniq -u > {params.miss}
    awk -F',' '{{print $0"\t"NF}}' {params.cov} | awk '{{OFS="\t"}}{{print $1"-"$2"-"$3,($7/$8),$4}}' | awk '{{gsub(/,/,","$1"\t"$2"\t")}}1' | tr ',' '\n' | awk '{{print $1"-"$3"\t"$2}}' | awk -v OFMT='%.5f' '{{OFS = "\t"}}{{sum[$1]+=$2;count[$1]++}}END{{for (i in sum) print i,(sum[i]/count[i])}}' > {output.kos}
    rev {output.kos} | cut -d'-' -f1 | rev | awk '{{a[$1]+=$2}}END{{for(i in a) print i,a[i]}}' > {output.KOunt}
    rm {params.bed}
    rm {params.bedko}
    rm {params.bedfull}
    rm {params.bedco}
    rm {params.rnahits}
    rm {params.cov}
    '''

rule nohit_fastq:
  input:
    touchfile=f"{outdir}Abundance/{{id}}/{{id}}_bwa",
    nohit=f"{outdir}barrnap/RNA_{{id}}_missing",
    gff=f"{outdir}Prodigal_Outputs/filtered_gffs/{{id}}/{{id}}.gff"
  params:
    bed=f"{outdir}diamond/{{id}}.bed",
    nobed=f"{outdir}barrnap/RNA_{{id}}_missing.bed",
    bam=f"{outdir}Abundance/{{id}}/{{id}}.bam.bed",
    fastq=f"{outdir}barrnap/RNA_{{id}}_missing.fastq.gz",
    reads=f"{outdir}barrnap/RNA_{{id}}_missing_reads",
    r1 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R1.fastq.gz",
    r2 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R2.fastq.gz"
  conda: "envs/bedtools_seqtk.yaml"
  output: f"{outdir}Touch/{{id}}_nohit_fastq"
  shell:
    '''
    awk -F '=|\t|;' '{{print $1 "_" $10 "\t" $4 "\t" $5}}' {input.gff} > {params.bed}
    awk 'NR==FNR{{a[$1]=$2;b[$1]=$3;next}}{{print $0 "\t" a[$1] "\t" b[$1]}}' {params.bed} {input.nohit} > {params.nobed}_tmp
    rm {params.bed}
    awk '{{print $1}}' {params.nobed}_tmp | rev | cut -d "_" -f 2- | rev > {params.nobed}_tmp.1
    awk '{{print $2 "\t" $3}}' {params.nobed}_tmp > {params.nobed}_tmp.2
    paste {params.nobed}_tmp.1 {params.nobed}_tmp.2 > {params.nobed}
    rm {params.nobed}_tmp
    rm {params.nobed}_tmp.1
    rm {params.nobed}_tmp.2
    bedtools intersect -a {params.nobed} -b {params.bam} -bed -wb | awk '{{print $7}}' | sort | uniq > {params.reads}
    rm {params.nobed}
    gawk -F '/' '$NF==1' {params.reads} > {params.reads}_R1
    gawk -F '/' '$NF==2' {params.reads} > {params.reads}_R2
    zcat {params.r1} | seqtk subseq - {params.reads}_R1 > {params.fastq}_R1
    zcat {params.r2} | seqtk subseq - {params.reads}_R2 > {params.fastq}_R2
    gzip {params.fastq}_R1
    gzip {params.fastq}_R2
    rm {params.reads}_R1
    rm {params.reads}_R2
    cat {params.fastq}_R1.gz {params.fastq}_R2.gz > {params.fastq}
    touch {output}
    '''

rule nohit_fastq_without_RNA:
  input:
    touchfile=f"{outdir}Abundance/{{id}}/{{id}}_bwa",
    nohit=f"{outdir}diamond/still_nohit_{{id}}_norna",
    gff=f"{outdir}Prodigal_Outputs/filtered_gffs/{{id}}/{{id}}.gff"
  params:
    bed=f"{outdir}diamond/{{id}}.bed",
    nobed=f"{outdir}diamond/no_RNA_{{id}}_missing.bed",
    bam=f"{outdir}Abundance/{{id}}/{{id}}.bam.bed",
    fastq=f"{outdir}diamond/no_RNA_{{id}}_missing.fastq.gz",
    reads=f"{outdir}diamond/no_RNA_{{id}}_missing_reads",
    r1 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R1.fastq.gz",
    r2 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R2.fastq.gz"
  conda: "envs/bedtools_seqtk.yaml"
  output: f"{outdir}Touch/{{id}}_nohit_fastq_without_RNA"
  shell:
    '''
    awk -F '=|\t|;' '{{print $1 "_" $10 "\t" $4 "\t" $5}}' {input.gff} > {params.bed}
    awk 'NR==FNR{{a[$1]=$2;b[$1]=$3;next}}{{print $0 "\t" a[$1] "\t" b[$1]}}' {params.bed} {input.nohit} > {params.nobed}_tmp
    rm {params.bed}
    awk '{{print $1}}' {params.nobed}_tmp | rev | cut -d "_" -f 2- | rev > {params.nobed}_tmp.1
    awk '{{print $2 "\t" $3}}' {params.nobed}_tmp > {params.nobed}_tmp.2
    paste {params.nobed}_tmp.1 {params.nobed}_tmp.2 > {params.nobed}
    rm {params.nobed}_tmp
    rm {params.nobed}_tmp.1
    rm {params.nobed}_tmp.2
    bedtools intersect -a {params.nobed} -b {params.bam} -bed -wb | awk '{{print $7}}' | sort | uniq > {params.reads}
    rm {params.nobed}
    gawk -F '/' '$NF==1' {params.reads} > {params.reads}_R1
    gawk -F '/' '$NF==2' {params.reads} > {params.reads}_R2
    zcat {params.r1} | seqtk subseq - {params.reads}_R1 > {params.fastq}_R1
    zcat {params.r2} | seqtk subseq - {params.reads}_R2 > {params.fastq}_R2
    gzip {params.fastq}_R1
    gzip {params.fastq}_R2
    rm {params.reads}_R1
    rm {params.reads}_R2
    cat {params.fastq}_R1.gz {params.fastq}_R2.gz > {params.fastq}
    touch {output}
    '''

rule nohit_annotate_reads:
  input:
    touch=f"{outdir}Touch/{{id}}_nohit_fastq",
    db=f"{diamond_db}",
    rna=f"{mmseq_db}",
    bg=f"{combined_bdg}"
  params:
    fq=f"{outdir}barrnap/RNA_{{id}}_missing.fastq.gz",
    of="6 qseqid sseqid pident qlen slen qstart qend sstart send evalue",
    ofm="query,target,pident,qlen,tlen,qstart,qend,tstart,tend,evalue",
    wd=f"{outdir}Abundance/{{id}}/nohit/",
    tmp=f"{outdir}diamond/{{id}}_reads_tmp",
    cat=f"{outdir}diamond/NoHit_{{id}}_cat",
    nohit=f"{outdir}diamond/{{id}}/{{id}}_hits_nohit",
    mmseq=f"{outdir}mmseqs/nohit_reads_{{id}}",
    reads=f"{outdir}barrnap/RNA_{{id}}_missing_reads",
    missing=f"{outdir}barrnap/RNA_{{id}}_still_missing_reads",
    threads=f"{nohit}"
  conda: "envs/diamond_bedtools_mmseq2.yaml"
  output:
    besthits=f"{outdir}diamond/{{id}}/{{id}}_nohit_besthits",
    ko=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_KOunts"
  shell:
    '''
    mkdir -p {params.wd}
    diamond blastx --query {params.fq} --threads {params.threads} --outfmt {params.of} -d {input.db} > {params.nohit}
    mmseqs easy-search {params.fq} {input.rna} {params.mmseq} {params.tmp} --threads {params.threads} --search-type 3 --format-output {params.ofm}
    cat {params.nohit} {params.mmseq} > {params.cat}
    rm -r {params.tmp}
    ./scripts/unmapped_reads_ko_100.sh {params.cat} {params.wd} {output.besthits} {input.bg} {output.ko}
    awk '{{print $1}}' {output.besthits} | sort | uniq | cat - {params.reads} | sort | uniq -u > {params.missing}
    rm {params.cat}
    '''

rule nohit_annotate_reads_without_RNA:
  input:
    touch=f"{outdir}Touch/{{id}}_nohit_fastq_without_RNA",
    db=f"{diamond_db}",
    rna=f"{mmseq_db}",
    bg=f"{combined_bdg}"
  params:
    fq=f"{outdir}diamond/no_RNA_{{id}}_missing.fastq.gz",
    of="6 qseqid sseqid pident qlen slen qstart qend sstart send evalue",
    wd=f"{outdir}Abundance/{{id}}/nohit/",
    nohit=f"{outdir}diamond/{{id}}/{{id}}_hits_nohit",
    reads=f"{outdir}diamond/no_RNA_{{id}}_missing_reads",
    missing=f"{outdir}diamond/no_RNA_{{id}}_still_missing_reads",
    threads=f"{nohit}"
  conda: "envs/diamond_bedtools_mmseq2.yaml"
  output:
    besthits=f"{outdir}diamond/{{id}}/{{id}}_nohit_without_RNA_besthits",
    ko=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_without_RNA_KOunts"
  shell:
    '''
    mkdir -p {params.wd}
    diamond blastx --query {params.fq} --threads {params.threads} --outfmt {params.of} -d {input.db} > {params.nohit}
    ./scripts/unmapped_reads_ko_100.sh {params.nohit} {params.wd} {output.besthits} {input.bg} {output.ko}
    awk '{{print $1}}' {output.besthits} | sort | uniq | cat - {params.reads} | sort | uniq -u > {params.missing}
    '''

rule kallisto:
  input:
    besthits=f"{outdir}diamond/{{id}}/{{id}}_nohit_besthits",
    ref=f"{kallisto}"
  params:
    missing=f"{outdir}barrnap/RNA_{{id}}_still_missing_reads",
    fq=f"{outdir}barrnap/RNA_{{id}}_missing.fastq.gz",
    r1=f"{outdir}barrnap/RNA_{{id}}_still_missing_R1.fastq.gz",
    r2=f"{outdir}barrnap/RNA_{{id}}_still_missing_R2.fastq.gz",
    threads=f"{kallisto_threads}",
    kal=directory(f"{outdir}kallisto/{{id}}_kallisto_still_missing")
  conda: "envs/kallisto.yaml"
  output:
    ko=f"{outdir}kallisto/{{id}}_kallisto_still_missing_KOunt"
  shell:
    '''
    gawk -F '/' '$NF==1' {params.missing} > {params.missing}_R1
    gawk -F '/' '$NF==2' {params.missing} > {params.missing}_R2
    zcat {params.fq} | seqtk subseq - {params.missing}_R1 | gzip > {params.r1}
    zcat {params.fq} | seqtk subseq - {params.missing}_R2 | gzip > {params.r2}
    rm {params.fq}
    rm {params.missing}_R1
    rm {params.missing}_R2
    kallisto quant -i {input.ref} -o {params.kal} -t {params.threads} -c {input.ref}_len {params.r1} {params.r2}
    awk -v len=$(cat {params.r1} {params.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum/(NR/4)}}') 'NR>1{{print $1"\t"(len*$4)*$2}}' {params.kal}/abundance.tsv > {params.kal}/tmp
    awk '{{print $1}}' {params.kal}/tmp | awk -F'-' '{{print NF-1}}' | paste {params.kal}/tmp - | awk -v OFMT='%.10f' -v OFS='\t' '{{print $1,($2/$3)}}' | awk '{{gsub(/-/,"\t"$2",")}}1' | cut -d, -f2- | sed 's/,$//g' | tr ',' '\n' | awk -v OFMT='%.10f' '{{a[$1]+=$2}}END{{for(i in a) print i,a[i]}}' > {output.ko}
    rm {params.kal}/tmp
    rm {params.r1}
    rm {params.r2}
    '''

rule unmapped_reads:
  input: 
    db=f"{diamond_db}",
    rna=f"{mmseq_db}",
    bg=f"{combined_bdg}",
    bam=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage"
  params:
    r1 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R1.fastq.gz",
    r2 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R2.fastq.gz",
    mapped=f"{outdir}Abundance/{{id}}/{{id}}_mapped_reads",
    unmapped=f"{outdir}Abundance/{{id}}/{{id}}_unmapped_reads",
    fastq=f"{outdir}Abundance/{{id}}/{{id}}_unmapped_reads.fastq.gz",
    bam=f"{outdir}Abundance/{{id}}/{{id}}/",
    of="6 qseqid sseqid pident qlen slen qstart qend sstart send evalue",
    ofm="query,target,pident,qlen,tlen,qstart,qend,tstart,tend,evalue",
    cat=f"{outdir}Abundance/{{id}}/{{id}}_catted",
    wd=f"{outdir}Abundance/{{id}}/unmapped/",
    sorted=f"{outdir}Abundance/{{id}}/unmapped/sorted",
    threads=f"{unmapped_threads}",
    tmp=f"{outdir}MMSeq_Outputs/{{id}}_tmp",
    unmap=f"{outdir}diamond/{{id}}/{{id}}_hits_unmapped",
    mm=f"{outdir}MMSeq_Outputs",
    mm2=f"{outdir}MMSeq_Outputs/{{id}}"
  conda: "envs/bedtools_seqtk.yaml"
  output:
    unmap=f"{outdir}diamond/{{id}}/{{id}}_hits_unmapped",
    mmseq=f"{outdir}MMSeq_Outputs/{{id}}/{{id}}_hits_unmapped",
    besthits=f"{outdir}diamond/{{id}}/{{id}}_unmapped_besthits",
    ko=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_unmapped_KOunt"
  shell:
    '''
    mkdir -p {params.wd}
    mkdir -p {params.mm}
    mkdir -p {params.mm2}
    mkdir -p {params.tmp}
    zcat {params.r1} {params.r2} | awk 'NR%4==1 {{print substr($1,2)}}' | cat - {params.mapped} | sort --parallel=8 | uniq -u > {params.unmapped}
    zcat {params.r1} {params.r2} | seqtk subseq - {params.unmapped} | gzip > {params.fastq}
    diamond blastx --query {params.fastq} --threads {params.threads} --outfmt 6 qseqid sseqid pident qlen slen qstart qend sstart send evalue -d {input.db} > {output.unmap}
    mmseqs easy-search {params.fastq} {input.rna} {output.mmseq} {params.tmp} --threads {params.threads} --search-type 3 --format-output {params.ofm}
    rm -r {params.tmp}
    cat {output.unmap} {output.mmseq} > {params.cat}
    ./scripts/unmapped_reads_ko_100.sh {params.cat} {params.wd} {output.besthits} {input.bg} {output.ko}
    rm {params.cat}
    '''

rule unmapped_reads_no_RNA:
  input: 
    db=f"{diamond_db}",
    bg=f"{combined_bdg}",
    cov=f"{outdir}Abundance_Outputs/{{id}}/{{id}}_coverage"
  params:
    r1 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R1.fastq.gz",
    r2 = f"{outdir}trimmed/{{id}}/{{id}}.trimmed.R2.fastq.gz",
    mapped=f"{outdir}Abundance/{{id}}/{{id}}_mapped_reads",
    unmapped=f"{outdir}Abundance/{{id}}/{{id}}_unmapped_reads",
    fastq=f"{outdir}Abundance/{{id}}/{{id}}_unmapped_reads.fastq.gz",
    of="6 qseqid sseqid pident qlen slen qstart qend sstart send evalue",
    cat=f"{outdir}Abundance/{{id}}/{{id}}_catted",
    wd=f"{outdir}Abundance/{{id}}/unmapped/",
    sorted=f"{outdir}Abundance/{{id}}/unmapped/sorted",
    threads=f"{unmapped_threads}"
  conda: "envs/bedtools_seqtk.yaml"
  output:
    unmap=f"{outdir}diamond/{{id}}/{{id}}_hits_unmapped_norna",
    besthits=f"{outdir}diamond/{{id}}/{{id}}_unmapped_besthits_norna",
    ko=f"{outdir}Kofamscan/Results/{{id}}/{{id}}_unmapped_KOunt_no_RNA"
  shell:
    '''
    mkdir -p {params.wd}
    zcat {params.r1} {params.r2} | awk 'NR%4==1 {{print substr($1,2)}}' | cat - {params.mapped} | sort | uniq -u > {params.unmapped}
    zcat {params.r1} {params.r2} | seqtk subseq - {params.unmapped} | gzip > {params.fastq}
    diamond blastx --query {params.fastq} --threads {params.threads} --outfmt 6 qseqid sseqid pident qlen slen qstart qend sstart send evalue -d {input.db} > {output.unmap}
    ./scripts/unmapped_reads_ko_100.sh {output.unmap} {params.wd} {output.besthits} {input.bg} {output.ko}
    '''

rule kallisto_unmapped:
  input:
    ref=f"{kallisto}",
    besthits=f"{outdir}diamond/{{id}}/{{id}}_unmapped_besthits"
  params:
    unmapped=f"{outdir}Abundance/{{id}}/{{id}}_unmapped_reads",
    still_missing=f"{outdir}Abundance/{{id}}/{{id}}_unmapped_still_missing",
    fq=f"{outdir}Abundance/{{id}}/{{id}}_unmapped_reads.fastq.gz",
    r1=f"{outdir}barrnap/RNA_{{id}}_unmap_still_missing_R1.fastq.gz",
    r2=f"{outdir}barrnap/RNA_{{id}}_unmap_still_missing_R2.fastq.gz",
    threads=f"{kallisto_threads}",
    kal=directory(f"{outdir}kallisto/{{id}}_unmapped_kallisto_still_missing")
  conda: "envs/kallisto.yaml"
  output:
    ko=f"{outdir}kallisto/{{id}}_unmapped_kallisto_still_missing_KOunt"
  shell:
    '''
    awk '{{print $1}}' {input.besthits} | sort | uniq | cat - {params.unmapped} | sort | uniq -u > {params.still_missing}
    gawk -F '/' '$NF==1' {params.still_missing} > {params.still_missing}_R1
    gawk -F '/' '$NF==2' {params.still_missing} > {params.still_missing}_R2
    zcat {params.fq} | seqtk subseq - {params.still_missing}_R1 | gzip > {params.r1}
    zcat {params.fq} | seqtk subseq - {params.still_missing}_R2 | gzip > {params.r2}
    rm {params.still_missing}_R1
    rm {params.still_missing}_R2
    kallisto quant -i {input.ref} -o {params.kal} -t {params.threads} -c {input.ref}_len {params.r1} {params.r2}
    awk -v len=$(zcat {params.r1} {params.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum/(NR/4)}}') 'NR>1{{print $1"\t"(len*$4)*$2}}' {params.kal}/abundance.tsv > {params.kal}/tmp
    awk '{{print $1}}' {params.kal}/tmp | awk -F'-' '{{print NF-1}}' | paste {params.kal}/tmp - | awk -v OFMT='%.10f' -v OFS='\t' '{{print $1,($2/$3)}}' | awk '{{gsub(/-/,"\t"$2",")}}1' | cut -d, -f2- | sed 's/,$//g' | tr ',' '\n' | awk -v OFMT='%.10f' '{{a[$1]+=$2}}END{{for(i in a) print i,a[i]}}' > {output.ko}
    rm {params.kal}/tmp
    rm {params.r1}
    rm {params.r2}
    '''

rule abundance_matrix_annotate_with_mapping:
  input:
    kofam=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_KOunt", id=IDS),
    nohit_p=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohitp_KOunt", id=IDS),
    rna=expand(f"{outdir}barrnap/{{id}}_rna_KOunt", id=IDS),
    kal=expand(f"{outdir}kallisto/{{id}}_kallisto_still_missing_KOunt", id=IDS),
    unmap=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_unmapped_KOunt", id=IDS),
    kal_unmap=expand(f"{outdir}kallisto/{{id}}_unmapped_kallisto_still_missing_KOunt", id=IDS)
  params:
    tmp=f"{outdir}tmp/",
    list=f"{db}ko_list",
    kegglist=f"{outdir}all_kos_found",
    id=expand(f"{{id}}", id=IDS)
  output:
    ann_results=f"{outdir}Results/All_KOunts_nohit_unmapped_default.csv",
  shell:
    '''
    mkdir {params.tmp}
    for i in {input.kofam}; do cp $i {params.tmp}; done
    for i in {input.nohit_p}; do cp $i {params.tmp}; done
    for i in {input.rna}; do cp $i {params.tmp}; done
    for i in {input.kal}; do cp $i {params.tmp}; done
    for i in {input.unmap}; do cp $i {params.tmp}; done
    for i in {input.kal_unmap}; do cp $i {params.tmp}; done
    sed -i '/NoHit/d' {params.tmp}/*
    if [ ! -e {params.kegglist} ]; then cat {params.tmp}* | awk '{{print $1}}' | sort | uniq > {params.kegglist}; fi
    for i in {params.id}; do echo $i; done > {params.tmp}samples
    while read file; do cat {params.tmp}"$file"_* | awk '{{a[$1]+=$2}}END{{for(i in a) print i,a[i]}}' > {params.tmp}"$file"_catted; done < {params.tmp}samples
    ./scripts/abundance_matrix_catted.sh {params.tmp} {params.list} {params.kegglist} > {output.ann_results}
    rm -r {params.tmp}
    '''

rule abundance_matrix_annotate_with_mapping_noclust:
  input:
    kofam=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noc_KOunt", id=IDS),
    nohit_p=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohitp_KOunt", id=IDS),
    rna=expand(f"{outdir}barrnap/{{id}}_rna_KOunt", id=IDS),
    kal=expand(f"{outdir}kallisto/{{id}}_kallisto_still_missing_KOunt", id=IDS),
    unmap=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_unmapped_KOunt", id=IDS),
    kal_unmap=expand(f"{outdir}kallisto/{{id}}_unmapped_kallisto_still_missing_KOunt", id=IDS)
  params:
    tmp=f"{outdir}tmp/",
    list=f"{db}ko_list",
    kegglist=f"{outdir}all_kos_found",
    id=expand(f"{{id}}", id=IDS)
  output:
    ann_results=f"{outdir}Results/All_KOunts_nohit_unmapped_no_clustering.csv",
  shell:
    '''
    mkdir {params.tmp}
    for i in {input.kofam}; do cp $i {params.tmp}; done
    for i in {input.nohit_p}; do cp $i {params.tmp}; done
    for i in {input.rna}; do cp $i {params.tmp}; done
    for i in {input.kal}; do cp $i {params.tmp}; done
    for i in {input.unmap}; do cp $i {params.tmp}; done
    for i in {input.kal_unmap}; do cp $i {params.tmp}; done
    sed -i '/NoHit/d' {params.tmp}/*
    if [ ! -e {params.kegglist} ]; then cat {params.tmp}* | awk '{{print $1}}' | sort | uniq > {params.kegglist}; fi
    for i in {params.id}; do echo $i; done > {params.tmp}samples
    while read file; do cat {params.tmp}"$file"_* | awk '{{a[$1]+=$2}}END{{for(i in a) print i,a[i]}}' > {params.tmp}"$file"_catted; done < {params.tmp}samples
    ./scripts/abundance_matrix_catted.sh {params.tmp} {params.list} {params.kegglist} > {output.ann_results}
    rm -r {params.tmp}
    '''

rule abundance_matrix_annotate_with_mapping_noRNA:
  input:
    kofam=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_noc_KOunt", id=IDS),
    nohit_p=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohitp_KOunt_norna", id=IDS),
    nohit_reads=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_nohit_without_RNA_KOunts", id=IDS),
    unmap=expand(f"{outdir}Kofamscan/Results/{{id}}/{{id}}_unmapped_KOunt_no_RNA", id=IDS),
  params:
    tmp=f"{outdir}tmp/",
    list=f"{db}ko_list",
    kegglist=f"{outdir}all_kos_found",
    id=expand(f"{{id}}", id=IDS)
  output:
    ann_results=f"{outdir}Results/All_KOunts_without_RNA.csv",
  shell:
    '''
    mkdir {params.tmp}
    for i in {input.kofam}; do cp $i {params.tmp}; done
    for i in {input.nohit_p}; do cp $i {params.tmp}; done
    for i in {input.nohit_reads}; do cp $i {params.tmp}; done
    for i in {input.unmap}; do cp $i {params.tmp}; done
    sed -i '/NoHit/d' {params.tmp}/*
    if [ ! -e {params.kegglist} ]; then cat {params.tmp}* | awk '{{print $1}}' | sort | uniq > {params.kegglist}; fi
    for i in {params.id}; do echo $i; done > {params.tmp}samples
    while read file; do cat {params.tmp}"$file"_* | awk '{{a[$1]+=$2}}END{{for(i in a) print i,a[i]}}' > {params.tmp}"$file"_catted; done < {params.tmp}samples
    ./scripts/abundance_matrix_catted.sh {params.tmp} {params.list} {params.kegglist} > {output.ann_results}
    rm -r {params.tmp}
    '''
