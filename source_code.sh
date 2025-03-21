# ==============================================================================
# TỔNG HỢP QUY TRÌNH PHÂN TÍCH RNA-SEQ
# ==============================================================================

# Cấu hình các biến
BASE_DIR="$HOME/Downloads/upstr_lung_hisat"
RAW_DATA_DIR="${BASE_DIR}/raw_data"
HISAT2_INDEX_DIR="${BASE_DIR}/ref_hisat2index/hg38/"
HISAT2_OUTPUT_DIR="${BASE_DIR}/hisat2_output"
CONVERT_SAM_TO_BAM_DIR="${BASE_DIR}/convert_sam_to_bam"
GENE_ANNOTATION_DIR="${BASE_DIR}/gene_annotation"
STRINGTIE_MERGE_DIR="${BASE_DIR}/stringtie_merge"
GTF_FILE_DIR="${BASE_DIR}/gtf_file"

SAMPLES=("8755796" "8755797" "8755798" "8755799" "8755800" "8755801" "8755802" "8755803" "8755804")

# ==============================================================================
# 1. MAP THE READS FOR EACH SAMPLE TO THE REFERENCE GENOME
# ==============================================================================

for i in "${SAMPLES[@]}"; do
    SAMPLE_NAME="SRR${i}"
    FASTQ_FILE="${RAW_DATA_DIR}/${SAMPLE_NAME}.fastq"
    OUTPUT_SAM="${HISAT2_OUTPUT_DIR}/${SAMPLE_NAME}.sam"
    
    if [ ! -f "$FASTQ_FILE" ]; then
        continue
    fi
    
    hisat2 -p 8 --dta -x "${HISAT2_INDEX_DIR}genome" -U "$FASTQ_FILE" -S "$OUTPUT_SAM"
done

# ==============================================================================
# 2. SORT AND CONVERT THE SAM FILES TO BAM
# ==============================================================================

for i in "${SAMPLES[@]}"; do
    SAMPLE_NAME="SRR${i}"
    SAM_FILE="${HISAT2_OUTPUT_DIR}/${SAMPLE_NAME}.sam"
    BAM_FILE="${CONVERT_SAM_TO_BAM_DIR}/${SAMPLE_NAME}.bam"
    
    samtools sort -@ 8 -o "$BAM_FILE" "$SAM_FILE"
done

# ==============================================================================
# 3. ASSEMBLE TRANSCRIPTS FOR EACH SAMPLE
# ==============================================================================

GTF_REF="${GENE_ANNOTATION_DIR}/hg38.ncbiRefSeq.gtf"

for i in "${SAMPLES[@]}"; do
    SAMPLE_NAME="SRR${i}"
    BAM_FILE="${CONVERT_SAM_TO_BAM_DIR}/${SAMPLE_NAME}.bam"
    GTF_OUTPUT="${GENE_ANNOTATION_DIR}/${SAMPLE_NAME}.gtf"
    
    stringtie -p 8 -G "$GTF_REF" -o "$GTF_OUTPUT" "$BAM_FILE"
done

# ==============================================================================
# 4. MERGE TRANSCRIPTS FROM ALL SAMPLES
# ==============================================================================

GTF_REF="${GENE_ANNOTATION_DIR}/hg38.ncbiRefSeq.gtf"
GTF_MERGED="${STRINGTIE_MERGE_DIR}/stringtie_merge.gtf"
MERGE_LIST="${GTF_FILE_DIR}/mergelist.txt"

stringtie --merge -p 8 -G "$GTF_REF" -o "$GTF_MERGED" "$MERGE_LIST"

# ==============================================================================
# 5. EXAMINE HOW THE TRANSCRIPTS COMPARE WITH THE REFERENCE ANNOTATION (OPTIONAL)
# ==============================================================================

REFERENCE_GTF="chrX_data/genes/chrX.gtf" # Thay thế bằng đường dẫn GTF tham chiếu của bạn
GFFCOMPARE_OUTPUT="${STRINGTIE_MERGE_DIR}/merged"

gffcompare -r "$REFERENCE_GTF" -G -o "$GFFCOMPARE_OUTPUT" "${STRINGTIE_MERGE_DIR}/stringtie_merge.gtf"

# ==============================================================================
# 6. ESTIMATE TRANSCRIPT ABUNDANCES AND CREATE TABLE COUNTS FOR BALLGOWN
# ==============================================================================

mkdir -p "${BASE_DIR}/ballgown"

for i in "${SAMPLES[@]}"; do
    SAMPLE_NAME="SRR${i}"
    BAM_FILE="${CONVERT_SAM_TO_BAM_DIR}/${SAMPLE_NAME}.bam"
    BALLGOWN_OUTPUT="${BASE_DIR}/ballgown/${SAMPLE_NAME}"

    mkdir -p "$BALLGOWN_OUTPUT"

    stringtie -e -B -p 8 -G "${STRINGTIE_MERGE_DIR}/stringtie_merge.gtf" \
              -o "${BALLGOWN_OUTPUT}/${SAMPLE_NAME}.gtf" "$BAM_FILE"
done

# ==============================================================================
# 7. COUNT READS WITH FEATURECOUNTS
# ==============================================================================

mkdir -p "${BASE_DIR}/counts"

for i in "${SAMPLES[@]}"; do
    SAMPLE_NAME="SRR${i}"
    BAM_FILE="${CONVERT_SAM_TO_BAM_DIR}/${SAMPLE_NAME}.bam"
    
    featureCounts -p -t exon -g gene_id -a "${GENE_ANNOTATION_DIR}/hg38.ncbiRefSeq.gtf" \
                  -o "${BASE_DIR}/counts/${SAMPLE_NAME}.counts" "$BAM_FILE"
done
