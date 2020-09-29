#!/usr/bin/env nextflow


// Define output directory
def exp_dir = null
// Check if an output dir name is provided in the config, otherwise use the selex_name instead.
if (params.output.out_dir == null) {
    today = new Date().format("yyyy-MM-dd")
    output_dir = "./output/$today/"
    exp_dir = new File("${output_dir}/${params.selex_name}/")

    // If the directory already exists a number in the range (0,9999) is appended
    for (int i = 0; exp_dir.exists(); i++) {
        nmbr = String.format("%04d", i)
        exp_dir = new File("${output_dir}/${params.selex_name}_${nmbr}/")
    }
} else {
    exp_dir = new File(params.output.out_dir)
    if (exp_dir.exists()) {
        throw new Exception("Specified output dir already exists.")
    }
}
def exp_data_dir = new File(exp_dir.getPath() + "/data/")
def exp_results_dir = new File(exp_dir.getPath() + "/results/")

def dir_trim = "1_trim_primers"
def dir_filter = "2_filtered"
def dir_merge = "3_merged"
def dir_unique = "4_unique"
def dir_analysis = "5_analysis"


// Create defined output dirs
exp_dir.mkdir()
exp_data_dir.mkdir()
exp_results_dir.mkdir()

// ===============================================================
// ===========  FILE READING  ====================================
// ===============================================================

// Reading file pairs (forward and reverse) in fastq format.
// Output is a flattended list of arrays: [[round_name, fwd_file, rev_file],...]
Channel
    .fromFilePairs(params.input.raw_reads, maxDepth:0, flat:true, checkIfExists:true, followLinks: true)
    { file -> file.getName().substring(0, file.getName().indexOf(params.input.round_delimiter)) }
    .filter { file -> params.rounds.contains(file[0]) }
    .set{ ch_raw_reads }

// ===============================================================
// =========== PRIMER CUTTING ====================================
// ===============================================================

def cutadapt_saveAs(f, rid) {
    disc = "discarded"
    switch(f) {
        case "no_adapters_found.fwd.fastq":
        case "no_adapters_found.rev.fastq":
        case "too_short.fwd.fastq":
        case "too_short.rev.fastq":
        case "too_long.fwd.fastq":
        case "too_long.rev.fastq":
            return "${disc}/${rid}/${rid}.${f}"

        case "fwd.fastq":
            return "${rid}.fwd.fastq"
        case "rev.fastq":
            return "${rid}.rev.fastq"
    }
}

process trim_selex_primers {
    publishDir exp_data_dir.getPath() + "/" + dir_trim, \
        mode: 'symlink', \
        overwrite: true, \
        pattern: "*.fastq", \
        saveAs: { fn -> cutadapt_saveAs(fn, read_id) }

    publishDir exp_results_dir.getPath() + "/" + dir_trim, \
        mode: 'copy', \
        overwrite: true, \
        pattern: "*.log"

    input:
        tuple read_id, file("${read_id}.fwd.fastq"), file("${read_id}.fastq") from ch_raw_reads

    output:
        tuple read_id, file("fwd.fastq"), file("rev.fastq") into ch_cutadapt_trimmed
        file("${read_id}.log") into ch_multiqc_cutadapt
        tuple read_id, file("no_adapters_found.fwd.fastq"),\
                    file("no_adapters_found.rev.fastq"),\
                    file("too_short.fwd.fastq"),\
                    file("too_short.rev.fastq"),\
                    file("too_long.fwd.fastq"),\
                    file("too_long.rev.fastq") into ch_cutadapt_discarded
        /* Quality Meta data */
        val(read_id) into ch_multiqc_sample_name
        tuple val(read_id), file("${read_id}.fwd.fastq"), file("fwd.fastq") into trace_read_processing_cutadapt

    script:
        def max_untrimmed = params.cutadapt.N + params.cutadapt.max_deviation 
        def min_untrimmed = params.cutadapt.N - params.cutadapt.max_deviation 
        def reads_raw = 0
        def reads_cut = 0        

        """        
        cutadapt \
            -e $params.cutadapt.max_error \
            --action=trim \
            --report=full \
            \
            --minimum-length $min_untrimmed \
            --maximum-length $max_untrimmed \
            -g ^${params.primers.p5_f}...${params.primers.p3_f} \
            -G ^${params.primers.p5_r}...${params.primers.p3_r}  \
            \
            --untrimmed-output no_adapters_found.fwd.fastq \
            --untrimmed-paired-output no_adapters_found.rev.fastq \
            --too-short-output too_short.fwd.fastq \
            --too-short-paired-output too_short.rev.fastq \
            --too-long-output too_long.fwd.fastq \
            --too-long-paired-output too_long.rev.fastq \
            \
            --output fwd.fastq \
            --paired-output rev.fastq \
            \
            ${read_id}.fwd.fastq \
            ${read_id}.fastq \
            > ${read_id}.log
        """
}

process trim_multiqc {
    echo true
    publishDir exp_results_dir.getPath() + "/" + dir_trim,
        mode: 'copy', \
        overwrite: true
    
    input:
        file("*") from ch_multiqc_cutadapt.collect()
        
    output:
        file("trimming_report*") into ch_multiqc_cutadapt_fin
        
    script:
    """
          multiqc --title 'Preprocessing of ${params.input.data_name}: Primer trimming' \
            -n trimming_report.html \
            .
    """
}


// ===============================================================
// =========== READ FILTERING ====================================
// ===============================================================

process filter_reads_by_read_quality {    
    publishDir exp_data_dir.getPath()  + "/" + dir_filter, \
        mode: 'symlink', \
        overwrite: true, \
        pattern: "filtered.*.fastq", \
        saveAs: { fn -> "${read_id}.${fn}" }
        
    publishDir exp_data_dir.getPath()  + "/" + dir_filter + "/discarded/", \
        mode: 'symlink', \
        overwrite: true, \
        pattern: "tfm_bad_reads.fastq", \
        saveAs: { fn -> "${read_id}.fwd.rev.bad_quality.fastq" }

    input:
        tuple read_id, file("${read_id}.fastq"), file("${read_id}.rev.fastq") from ch_cutadapt_trimmed
    output:
        tuple read_id, file("filtered.fwd.fastq"), file("filtered.rev.fastq") into ch_fastp_filtered
        tuple read_id, file("tfm_bad_reads.fastq") into ch_fastp_filter_discarded
        file("${read_id}.filter.fastp.json") into ch_fastp_filter_log
        tuple val(read_id), file("filtered.fwd.fastq") into trace_read_processing_filter
    script:
        def reads_filtered = 0
        """
        fastp -i ${read_id}.fastq -I ${read_id}.rev.fastq \
            -o filtered.fwd.fastq -O filtered.rev.fastq \
            --failed_out=tfm_bad_reads.fastq \
            --disable_adapter_trimming \
            --average_qual $params.fastp.filter_min_phred \
            --json=${read_id}.filter.fastp.json
        """
}

process filter_multiqc {
    echo true
    publishDir exp_results_dir.getPath() + "/" + dir_filter,
        mode: 'copy', \
        overwrite: true
    
    input:
        file("*") from ch_fastp_filter_log.collect()
        
    output:
        file("filter_report*") into ch_multiqc_filter_fin
        
    script:
    """
          multiqc --title 'Preprocessing of ${params.input.data_name}: Quality Filtering' \
            -n filter_report.html \
            .
    """
}


// ===============================================================
// ===========  READ MERGING  ====================================
// ===============================================================
process merge_reads {    
    publishDir exp_data_dir.getPath() + "/" + dir_merge, \
        mode: 'symlink', \
        overwrite: true, \
        pattern: "merged.fastq", \
        saveAs: { fn -> "${read_id}.fastq" }
    
    publishDir exp_data_dir.getPath()  + "/" + dir_merge + "/discarded/", \
        mode: 'symlink', \
        overwrite: true, \
        pattern: "*disc*", \
        saveAs: { fn -> "${read_id}.${fn}" }

    input:
        tuple read_id, file("${read_id}.fastq"), file("${read_id}.rev.fastq") from ch_fastp_filtered
    output:
        tuple read_id, file("merged.fastq") into ch_merged
        tuple read_id, file("discarded.fwd.fastq"), file("discarded.rev.fastq") into ch_fastp_merging_discarded
        tuple read_id, file("${read_id}.merging.fastp.json") into ch_fastp_merging_log
        tuple val(read_id), file("merged.fastq") into trace_read_processing_merge
    script:
        """
        mkdir ./fastp_meta/
        fastp -i ${read_id}.fastq -I ${read_id}.rev.fastq \
            -o discarded.fwd.fastq -O discarded.rev.fastq \
            --merge \
            --merged_out=merged.fastq \
            --disable_adapter_trimming \
            --json=${read_id}.merging.fastp.json

        lines_merged=\$(cat merged.fastq | wc -l)
        reads_merged=\$((lines_merged/4))
        """
}

process merge_multiqc {
    echo true
    publishDir exp_results_dir.getPath()  + "/" + dir_merge,
        mode: 'copy', \
        overwrite: true
    
    input:
        file("*") from ch_fastp_merging_log.collect()
        
    output:
        file("merge_report*") into ch_multiqc_merge_fin
        
    script:
    """
          multiqc --title 'Preprocessing of ${params.input.data_name}: Read Merging' \
            -n merge_report.html \
            .
    """
}


trace_read_processing_combined = trace_read_processing_cutadapt
    .join(trace_read_processing_filter)
    .join(trace_read_processing_merge)


process trace_read_processing_count_reads {
    input:
        tuple val(round_id), file("raw.fastq"), file("cut.fastq"), file("filtered.fastq"), file("merged.fastq") from trace_read_processing_combined
    output:
        tuple val(round_id), file("${round_id}.counts.csv") into trace_read_processing_counts
    script:
    """
        raw_reads=\$(cat raw.fastq | wc -l)
        raw_reads=\$((raw_reads/4))
    
        cut_reads=\$(cat cut.fastq | wc -l)
        cut_reads=\$((cut_reads/4))

        filtered_reads=\$(cat filtered.fastq | wc -l)
        filtered_reads=\$((filtered_reads/4))

        merged_reads=\$(cat merged.fastq | wc -l)
        merged_reads=\$((merged_reads/4))


        touch ${round_id}.counts.csv
        echo -n $round_id >> ${round_id}.counts.csv
        echo -n '\t' >> ${round_id}.counts.csv
        echo -n \$raw_reads >> ${round_id}.counts.csv
        echo -n '\t' >> ${round_id}.counts.csv
        echo -n \$cut_reads >> ${round_id}.counts.csv
        echo -n '\t' >> ${round_id}.counts.csv
        echo -n \$filtered_reads >> ${round_id}.counts.csv
        echo -n '\t' >> ${round_id}.counts.csv
        echo -n \$merged_reads >> ${round_id}.counts.csv
        echo >> ${round_id}.counts.csv # new line
    """

}

trace_read_processing_counts_sorted = trace_read_processing_counts
    .toSortedList( { a -> a[0] } )
    .transpose()
    .last()
    .collect()

process trace_read_processing_combine_to_csv {
    publishDir exp_results_dir.getPath(), \
        mode: 'copy', \
        overwrite: true 

    input:
        file(f) from trace_read_processing_counts_sorted
    output:
        file("reads_filtering.csv") into trace_read_counts
    script:
    """
        touch reads_filtering.csv
        echo 'round_id\traw\tcutadapt\tfilter\tmerge' >> reads_filtering.csv
        cat $f >> reads_filtering.csv
    """
}

// ===============================================================
// =========== Analysis of SELEX Data ============================
// ===============================================================

process merged_fastq_to_fasta {
    input:
        tuple val(round_id), file("selex_round.fastq") from ch_merged
    output:
        tuple val(round_id), file("${round_id}.fasta") into ch_merged_fasta

    script:
    """
        sed -n 'p;n;p;n;n' selex_round.fastq | sed 's/@M/>M/g' > ${round_id}.fasta
    """
}

ch_merged_fasta.into { ch_merged_fasta_1; ch_merged_fasta_2 }

ch_merged_fasta_flat = ch_merged_fasta_1
    .toSortedList { read_id -> read_id[0] }
    .transpose()
    .last()
    .collect()


process derep_fasta_files {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis, \
        mode: 'copy', \
        overwrite: true, \
        pattern: "*.csv"
    publishDir exp_data_dir.getPath()  + "/" + dir_analysis, \
        mode: 'copy', \
        overwrite: true, \
        pattern: "*.fasta"

    input:
        file(f) from ch_merged_fasta_flat
    output:
        file("selex_derep.csv") into ch_derep_csv
        file("selex_derep.fasta") into ch_derep_fasta

    script:
    """
        SELEXderep.py -o selex_derep.fasta -c selex_derep.csv -t ${params.specs.cpus} --add-primers --primer-forward ${params.primers.p5_f} --primer-reverse ${params.primers.p3_f} ${f}
    """

}

/*process fold_sequences {
    storeDir "store/"

    input:
        file("selex_derep.fasta") from ch_derep_fasta
    output:
        file("selex_seqs.folded") into ch_vienna_rna_fold

    script:
    """
        head -n 100 selex_derep.fasta > selex_derep100.fasta
        RNAfold -T ${params.vienna_rna.T} --paramFile=${params.vienna_rna.mathews2004} \
            --noPS --gquad --noconv --jobs=${params.specs.cpus} \
            -i selex_derep100.fasta \
            --outfile=selex_seqs.folded
    """
}

awk = 'BEGIN { RS = ">"; FS = "\\n|( \\\\( *)|)\\n" } { print $1, "\\t", $2, "\\t", $3, "\\t", $4 }'
process folded_sequences_to_csv {
    storeDir "store/"

    input:
        file("selex_seqs.folded") from ch_vienna_rna_fold

    output:
        file("selex_seqs.csv") into ch_vienna_rna_csv

    script:
    """
        awk '$awk' selex_seqs.folded > selex_seqs.csv
    """
}*/


process assess_rpm {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis, \
        mode: 'copy', \
        overwrite: true, \
        pattern: "*"

    input:
        file("selex_derep.csv") from ch_derep_csv
    output:
        file("rpm.csv") into ch_rpm

    script:
    """
        SELEXrpm.r -i selex_derep.csv -o rpm.csv
    """
}

process assess_top1000_rpm {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis,\
        mode: 'copy', \
        overwrite: true, \
        pattern: "*"

    input:
        file("rpm.csv") from ch_rpm
    output:
        file("top1000.xlsx") into ch_rpm_top1000

    script:
    """
        SELEXrpm_top1000.R -i rpm.csv -o top1000.csv -n 1000
    """
}

logs = [2,10]
process assess_rpm_log_duplicates {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis, \
        mode: 'copy', \
        overwrite: true, \
        pattern: "*"

    input:
        file("rpm.csv") from ch_rpm
        each logX from logs
    output:
        file("$logX") into ch_rpm_csv
        tuple val(logX), file("$logX/rpm_percent.log_duplicates.csv") into ch_rpm_plot

    script:
    """
        mkdir $logX
        SELEXanalyse_log_duplicates.r --in-rpm-csv rpm.csv --bin-base $logX --out-dir $logX
        plot_rpm.r $logX/rpm_percent.log_duplicates.csv $logX $logX
    """
}

/*process plot_rpm_log_duplicates {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis + "/" + logX, \
        mode: 'copy', \
        overwrite: true, \
        pattern: "*"

    input:    
        tuple val(logX), file("rpm_stp.csv") from ch_rpm_plot
    output:
        file("*.tiff") into ch_rpm_tiff

    script:
    """
        # mkdir $logX
        # plot_rpm.r rpm_stp.csv $logX $logX
        plot_rpm.r rpm_stp.csv $logX .
    """
}*/


process assess_selex_composition {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis + "/composition/", \
        mode: 'copy', \
        overwrite: true

    input:
        file(f) from ch_merged_fasta_flat
    output:
        file("selex.acgt.csv") into ch_selex_composition_plot

    script:
    """
        analyse_selex_composition.r ${params.n} ${f}
    """
}

process plot_selex_composition {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis + "/composition/", \
        mode: 'copy', \
        overwrite: true

    input:
        file("selex.acgt.csv") from ch_selex_composition_plot
    output:
        file("*.tiff") into ch_selex_composition_tiff


    script:
    """
        plot_selex_composition.r selex.acgt.csv
    """
}


process assess_round_composition {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis + "/composition/", \
        mode: 'copy', \
        overwrite: true

    input:
        tuple val(round_id), file(f) from ch_merged_fasta_2
    output:
        tuple val(round_id), file("${round_id}.acgt.csv") into ch_acgt_compositions_csv

    script:
    """
        SELEXntcomposition.py -i ${f} -o ${round_id}.acgt.csv -n ${params.n} --DNA -p
    """
}

process plot_round_composition {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis + "/composition/", \
        mode: 'copy', \
        overwrite: true

    input:
        tuple val(round_id), file("${round_id}.acgt.csv") from ch_acgt_compositions_csv
    output:
        file("${round_id}.acgt.tiff") into ch_acgt_compositions_tiff

    script:
    """
        plot_round_composition.r ${round_id}.acgt.csv ${round_id}.acgt.tiff
    """
}


process assess_duplicates {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis + "/composition/", \
        mode: 'copy', \
        overwrite: true

    input:
        file("selex_derep.csv") from ch_derep_csv
    output:
        file("duplicates.csv") into ch_duplicates_csv

    script:
    """
        SELEXanalyse_duplicates.r -o duplicates.csv -i selex_derep.csv
    """
}

process plot_duplicates {
    publishDir exp_results_dir.getPath()  + "/" + dir_analysis + "/composition/", \
        mode: 'copy', \
        overwrite: true

    input:
        file("duplicates.csv") from ch_duplicates_csv
    output:
        file("duplicates.tiff") into ch_duplicates_tiff

    script:
    """
        plot_duplicates.r duplicates.csv duplicates.tiff
    """
}


