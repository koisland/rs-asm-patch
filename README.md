# rs-asm-patch
Patch two assemblies using misassembly calls and alignment between two assemblies.

> WIP

Align two assemblies with minimap2 on asm5 preset.
```bash
minimap2 -t 12 -ax asm5 --eqx --cs -s 25000 -K 8G mPanPan1_merged_dedup_asm.fa.gz mPanPan1_merged_dedup_asm_query.fa.gz > mPanPan1.paf
```

Then trimmed with rustybam trim.
```bash
# Not included.
rb trim-paf mPanPan1.paf | sed 's/\\tid:Z://g' > mPanPan1_trim.paf
```

Generate misassembly calls with nucflag.
```bash
data/mPanPan1_cen_misassemblies.bed
data/mPanPan1_cen_misassemblies_query.bed
```

### Test
Run script on mPanPan1.
* Verkko reference
* hifiasm query

```bash
cargo run -- \
-i data/mPanPan1_trim.paf \
-r data/mPanPan1_merged_dedup_asm.fa.gz \
-q data/mPanPan1_merged_dedup_asm_query.fa.gz \
--ref-misasm-bed <(grep -v "HET" data/mPanPan1_cen_misassemblies.bed) \
--qry-misasm-bed <(grep -v "HET" data/mPanPan1_cen_misassemblies_query.bed) \
--log-level Debug \
-b test.bed \
-o test.fa 2> test.log
```
