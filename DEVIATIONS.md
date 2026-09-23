# Deviations from the pre-specified analysis plan

1. **23 September 2026, before results.** P1 gene set added: `cheon2014_collagen` (10 genes; Cheon et al. 2014, Clin Cancer Res 20(3):711-723, doi 10.1158/1078-0432.CCR-13-1256). The Cheon signature named in plan v1.0 was not previously in the project gene-set JSON.
2. **23 September 2026, before results.** P2 gene set changed from `resistance_up_core` to `resistance_up_intrinsic` (14 genes: `resistance_up_core` without FOS, JUN, NR4A1, COL11A1, LOX, FABP4, IL6, STAT3, STAB1), to avoid overlap with the stress positive control (`stress_IEG_state_GSE165897`) and with stromal genes.
3. **23 September 2026, before results.** Gene sets with fewer than 5 genes (`matrix_CAF`, `sensitivity_up`, `partial_EMT_Xu2022`) moved from set-level testing to gene-level reporting. The 5-gene minimum is unchanged.
4. **23 September 2026, before results.** liana-py (rank_aggregate by sample) chosen as the cell-cell communication tool for secondary analysis 7.
