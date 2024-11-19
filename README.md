# rs-asm-patch
Patch two assemblies using misassembly calls and alignment between two assemblies.

> [!NOTE]
> WIP

Align two assemblies with minimap2 on asm5 preset.
```bash
minimap2 -t 12 -ax asm5 --eqx -K 8G mPanPan1_merged_dedup_asm.fa.gz mPanPan1_merged_dedup_asm_query.fa.gz > mPanPan1.paf
```

Align reads to reference assembly and generate misassembly calls with [NucFlag](https://github.com/logsdon-lab/NucFlag).
* For workflow, see [`Snakemake-NucFlag`](https://github.com/logsdon-lab/Snakemake-NucFlag)

Use [`impg`](https://github.com/pangenome/impg) to liftover query coordinates to reference coordinates and substitute query sequence.
* [Fork](https://github.com/koisland/impg/tree/dev) allowing reverse intervals used.

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
