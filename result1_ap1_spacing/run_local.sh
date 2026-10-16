#!/bin/bash
# Result 1: local steps after the Kaggle notebooks (run from code/ on Linux/WSL).
# Needs: python (numpy, pandas, scipy), UCSC twoBitToFa and faToTwoBit, MEME suite FIMO 5.5.x.
#
# Before running, extract the notebook outputs (enhancers_results.tgz -> enhancers/) into
#   data/result1/zebrafish_heart/  data/result1/killifish_fin/  data/result1/zebrafish_fin/
#   data/result1/mouse_p1/         data/result1/mouse_p8/
set -e
S=result1_ap1_spacing; D=data/result1; G=data/genomes; O=output/result1
mkdir -p $G $O

# ---- reference files ----
U=https://hgdownload.soe.ucsc.edu/goldenPath/danRer10/bigZips
[ -s $G/danRer10.2bit ] || wget -q $U/danRer10.2bit -O $G/danRer10.2bit
[ -s $G/danRer10.chrom.sizes ] || wget -q $U/danRer10.chrom.sizes -O $G/danRer10.chrom.sizes
if [ ! -s $G/killifish.2bit ]; then   # Nothobranchius furzeri, NCBI RefSeq GCF_043380555.1 (sequence names = accessions)
  A=GCF_043380555.1
  wget -q "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$A/download?include_annotation_type=GENOME_FASTA" -O $G/kf.zip
  python -c "import zipfile; zipfile.ZipFile('$G/kf.zip').extractall('$G/kf')"
  awk '/^>/{print $1; next}{print}' $G/kf/ncbi_dataset/data/$A/*.fna > $G/killifish.fa
  faToTwoBit $G/killifish.fa $G/killifish.2bit
  awk '/^>/{if (n) print n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)} END{print n"\t"l}' $G/killifish.fa > $G/killifish.sizes
  rm -rf $G/kf.zip $G/kf $G/killifish.fa
fi
B=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/835/505/GCF_016835505.1_ASM1683550v1/GCF_016835505.1_ASM1683550v1
[ -s $G/polypterus.fna.gz ] || wget -q ${B}_genomic.fna.gz -O $G/polypterus.fna.gz
[ -s $D/polypterus/GSE315999_Ps_ATAC_seq_diff_peaks.txt.gz ] || { mkdir -p $D/polypterus; wget -q -P $D/polypterus \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE315nnn/GSE315999/suppl/GSE315999_Ps_ATAC_seq_diff_peaks.txt.gz; }
[ -s $G/jaspar2024_vert_nr.meme ] || wget -q \
  https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt -O $G/jaspar2024_vert_nr.meme

ZH="$D/zebrafish_heart/rre.bed $D/zebrafish_heart/control_pool.bed $G/danRer10.2bit $G/danRer10.chrom.sizes"
KF="$D/killifish_fin/rre.bed $D/killifish_fin/control_pool.bed $G/killifish.2bit $G/killifish.sizes"
ZF="$D/zebrafish_fin/rre.bed $D/zebrafish_fin/control_pool.bed $G/danRer10.2bit $G/danRer10.chrom.sizes"

# step 02: primary spacing test for the two confirmatory sets (Figure 1A rows 1-2)
python $S/02_ap1_spacing_test.py --rre $D/zebrafish_heart/rre.bed --ctrl $D/zebrafish_heart/control_pool.bed \
  --genome-2bit $G/danRer10.2bit --genome-sizes $G/danRer10.chrom.sizes --out $O/zebrafish_heart_spacing.json
python $S/02_ap1_spacing_test.py --rre $D/killifish_fin/rre.bed --ctrl $D/killifish_fin/control_pool.bed \
  --genome-2bit $G/killifish.2bit --genome-sizes $G/killifish.sizes --out $O/killifish_fin_spacing.json

# step 03: repeat-masked test (two confirmatory sets) and the zebrafish fin replication set (Figure 1A row 3)
python $S/03_spacing_sets_and_repeat_masking.py --set zebrafish_heart $ZH --set killifish_fin $KF --out $O/repeat_masked_zebrafish_heart_killifish_fin.json
python $S/03_spacing_sets_and_repeat_masking.py --set zebrafish_fin $ZF --out $O/zebrafish_fin_spacing.json

# step 04: local composition null (pooled over the two confirmatory sets)
python $S/04_local_composition_null.py --set zebrafish_heart $ZH --set killifish_fin $KF --set zebrafish_fin $ZF --primary zebrafish_heart killifish_fin \
  --out $O/local_composition_null.json

# step 05: JASPAR motif panel (needs FIMO; slow)
python $S/05_jaspar_motif_panel.py --meme $G/jaspar2024_vert_nr.meme --workdir $O/jaspar_work --jobs 4 --chunks 24 \
  --set zebrafish_heart $ZH --set killifish_fin $KF --set zebrafish_fin $ZF --out $O/jaspar_motif_panel.json

# step 06 + 03: Polypterus (Figure 1A row 4); primary and two sensitivity definitions
python $S/06_polypterus_prepare_units.py --table $D/polypterus/GSE315999_Ps_ATAC_seq_diff_peaks.txt.gz \
  --genome $G/polypterus.fna.gz --out $D/polypterus/units
for s in rre rre_median rre_fold1; do
  n=polypterus; [ $s = rre ] || n=polypterus_$s
  python $S/03_spacing_sets_and_repeat_masking.py --set $n $D/polypterus/units/$s.bed $D/polypterus/units/control_pool.bed \
    $D/polypterus/units/units.2bit $D/polypterus/units/units.sizes --out $O/polypterus_spacing_$s.json
done

# step 07: mouse heart P1 / P8 (Figure 1A rows 5-6)
python $S/07_mouse_p1_p8_summary.py --p1 $D/mouse_p1 --p8 $D/mouse_p8 --out $O/mouse_heart_p1_p8_spacing.json

# step 08: spacing vs strength of induction (pooled over the four fish sets)
python $S/08_spacing_vs_induction.py \
  --set zebrafish_heart $D/zebrafish_heart/rre.bed $G/danRer10.2bit $G/danRer10.chrom.sizes \
  --set killifish_fin $D/killifish_fin/rre.bed $G/killifish.2bit $G/killifish.sizes \
  --set zebrafish_fin $D/zebrafish_fin/rre.bed $G/danRer10.2bit $G/danRer10.chrom.sizes \
  --set polypterus $D/polypterus/units/rre.bed $D/polypterus/units/units.2bit $D/polypterus/units/units.sizes \
  --out $O/spacing_vs_induction.json

# step 09 (ChromBPNet) runs on a Kaggle GPU: notebooks/chrombpnet_ap1_spacing.ipynb
echo "Result 1 local steps done -> $O"
